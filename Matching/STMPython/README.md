# Matching with Python

## Install

Note: a C++ compiler is required!

Install [PDM](https://pdm-project.org) (for example with pipx as explained 
[here](https://fluidhowto.readthedocs.io/en/latest/setup/setup-apps.html)) and run (deactivate any external Python virtual environments) in the folder of the repository :

```sh
make
```

A virtual env `./.venv` should have been created. One can activate it with

```
. .venv/bin/activate
```
 Another way to use PDM is to add at the beginning of each python command : 

```
pdm run [the usual command]
```

## Run tests

```sh
make test
```

## Run the bench

```sh
make bench
make bench_unoptimized
make bench_cpp
```

