# Msg123 : Multi-scale groundwater model 1phase 2resolution 3dimension
Msg123 is a variably saturated flow model for large-scale problems.

# Quick Start on Unix/Linux
This quickstart shows how to build and run the model with a minimal example.

### Requirements
- Linux or macOS
- Fortran compiler (gfortran >= 9.1.0 or interl fortran >=19.0.3.199)

### Clone
```bash
git clone git@github.com:Msg123-Dev/Msg123.git
cd Msg123
```

## Build
```bash
cd make
vim makefile
```

Define the Fortran compile flags (gfortran or ifort)  
Define the mode (release or debug)  
Define the parallel (MPI, OpenMP, Hybrid)

```bash
make
```

## Usage
Run a minimal example (Coming soon)

## Documentation
Coming soon

## Citation
If you use this model, please cite:  
See CITATION.cff.

## License
This project is licensed under the Apache-2.0 License.
