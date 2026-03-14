FROM ubuntu:22.04

# Avoid interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install build dependencies
# - gfortran        : Fortran compiler (all modes)
# - openmpi-bin     : mpirun/mpiexec (MPI, Hybrid)
# - libopenmpi-dev  : mpif90 wrapper and MPI headers (MPI, Hybrid)
# - make            : build tool
RUN apt-get update && apt-get install -y \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \
    make \
    && rm -rf /var/lib/apt/lists/*

# Build arguments
# PARA: (empty, default = serial) | OpenMP | MPI | Hybrid
# MODE: release (default) | debug
ARG PARA=
ARG MODE=release

# Set working directory
WORKDIR /work

# Copy source code and makefiles
COPY src/ ./src/
COPY make/ ./make/

# Build
# For MPI and Hybrid, FC must be mpif90 and FFLAGS must be set explicitly,
# because makefile's ifeq($(FC),gfortran) block is skipped when FC=mpif90.
# The makefile outputs the executable to ../msg123 relative to the make directory,
# which resolves to /work/msg123.
WORKDIR /work/make
RUN if [ "${PARA}" = "MPI" ]; then \
        if [ "${MODE}" = "debug" ]; then \
            FFLAGS="-cpp -Wall -g -fbacktrace -fbounds-check -Wuninitialized -pedantic -ffpe-trap=invalid,zero,overflow"; \
        else \
            FFLAGS="-cpp -Wall -O3 -fmax-stack-var-size=65535"; \
        fi; \
        make FC=mpif90 PARA=${PARA} MODE=${MODE} \
             FFLAGS="${FFLAGS}" \
             CPPFLAGS="-DMPI_MSG"; \
    elif [ "${PARA}" = "Hybrid" ]; then \
        if [ "${MODE}" = "debug" ]; then \
            FFLAGS="-cpp -Wall -g -fbacktrace -fbounds-check -Wuninitialized -pedantic -ffpe-trap=invalid,zero,overflow -fopenmp"; \
        else \
            FFLAGS="-cpp -Wall -O3 -fmax-stack-var-size=65535 -fopenmp"; \
        fi; \
        make FC=mpif90 PARA=${PARA} MODE=${MODE} \
             FFLAGS="${FFLAGS}" \
             CPPFLAGS="-DMPI_MSG"; \
    else \
        make FC=gfortran PARA=${PARA} MODE=${MODE}; \
    fi

# Install the executable to PATH so it can be called from any working directory
RUN ls -la /work/msg123 && mv /work/msg123 /usr/local/bin/msg123 && chmod +x /usr/local/bin/msg123

# Default working directory for the user
WORKDIR /work

# Default command
CMD ["msg123"]
