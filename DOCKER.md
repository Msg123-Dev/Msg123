# Running Msg123 with Docker

## Requirements

- [Docker](https://docs.docker.com/get-docker/) installed on your system

## Build

Build the Docker image from the repository root. Use `--build-arg` to select the
parallelization mode. Each mode produces a separate image.

### Serial (default)

```bash
docker build -t msg123 .
```

### OpenMP (shared-memory parallel)

```bash
docker build --build-arg PARA=OpenMP -t msg123-omp .
```

### MPI (distributed-memory parallel)

```bash
docker build --build-arg PARA=MPI -t msg123-mpi .
```

### Hybrid (MPI + OpenMP)

```bash
docker build --build-arg PARA=Hybrid -t msg123-hybrid .
```

### Debug build

Combine `MODE=debug` with any `PARA` option:

```bash
docker build --build-arg PARA=OpenMP --build-arg MODE=debug -t msg123-omp-debug .
```

## Run

Run the container from the directory that contains your input files (`msg123.main` and others).
Use `-it` to enter the container interactively.

> **Note:** Always use `"$PWD"` (absolute path) with `-v`. Relative paths such as `./` do not work.

### Serial

```bash
cd /path/to/your/data
docker run --rm -it -v "$PWD":/work msg123 bash
```

Inside the container:

```bash
msg123
exit
```

### OpenMP

Set the number of threads with `OMP_NUM_THREADS`:

```bash
cd /path/to/your/data
docker run --rm -it -e OMP_NUM_THREADS=4 -v "$PWD":/work msg123-omp bash
```

Inside the container:

```bash
msg123
exit
```

### MPI

```bash
cd /path/to/your/data
docker run --rm -it -v "$PWD":/work msg123-mpi bash
```

Inside the container:

```bash
mpirun --allow-run-as-root -np 4 msg123
exit
```

### Hybrid (MPI + OpenMP)

```bash
cd /path/to/your/data
docker run --rm -it -e OMP_NUM_THREADS=2 -v "$PWD":/work msg123-hybrid bash
```

Inside the container:

```bash
mpirun --allow-run-as-root -np 4 msg123
exit
```

## Notes

- Inside the container, the current directory is `/work`, which maps to your data directory on the host.
- Output files written to `/work` inside the container will appear in your data directory on the host.
- To confirm the executable is available: `which msg123` → `/usr/local/bin/msg123`

## Build argument summary

| Argument | Options | Default |
|----------|---------|---------|
| `PARA` | *(empty)*, `OpenMP`, `MPI`, `Hybrid` | *(empty, serial)* |
| `MODE` | `release`, `debug` | `release` |
