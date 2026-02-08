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
make -j example
./bin/examples/pke/example
```

In the example, we demonstrate the results of arithmetic (add/multiplication) and boolean operations (bitwise logic, shift/rotation, compare). We also demonstrate bootstrappings. The noise growth after each operation is also printed. An example log can be checked in [example.log](/log/example.log). Note that for demonstration purpose we use ring dimension of 1024.

If you want to benchmark, use the following program whose RLWE parameter satisfies 128-bit security.

```bash
make -j benchmark-full
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 64 bench
```

Please note that this repository is a prototype and there is no warranty on its reliability.
