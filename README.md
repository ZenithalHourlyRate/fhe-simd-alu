Installation Instructions
=====================================

We enable OpenMP, NTL, tcmalloc and HEXL in OpenFHE.

The needed packages in Debian/Ubuntu:

```bash
sudo apt install build-essential git libntl-dev libgmp-dev cmake autoconf libtool clang libomp5 libomp-dev
```

Now we can build and run an example

```bash
mkdir build; cd build
CC=clang CXX=clang++ cmake -DWITH_INTEL_HEXL=ON -DMATHBACKEND=6 -DWITH_NTL=ON -DWITH_TCM=ON ..
make -j tcm
make -j z-example
./bin/examples/pke/z-example
```

In the example, we demonstrate the results of arithmetic (add/multiplication) and boolean operations (bitwise logic, shift/rotation, compare). We also demonstrate bootstrappings. The noise growth after each operation is also printed.

If you want to benchmark, use

```bash
make -j z-benchmark-full
OMP_NUM_THREADS=1 ./bin/examples/pke/z-benchmark-full 64 bench
```
