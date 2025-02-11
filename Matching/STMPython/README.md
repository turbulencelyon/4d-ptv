# Matching with Python

## Install

Note: a C++ compiler is required!

Install [PDM](https://pdm-project.org) and run (deactivate any external Python virtual environments):

```sh
make
```

A virtual env `./.venv` should have been created. One can activate it with

```
. .venv/bin/activate
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
