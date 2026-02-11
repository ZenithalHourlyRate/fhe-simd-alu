#!/bin/bash
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 64 bench | tee bench-64.log
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 32 bench | tee bench-32.log
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 16 bench | tee bench-16.log
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 8 bench | tee bench-8.log
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 128 bench | tee bench-128.log
OMP_NUM_THREADS=1 ./bin/examples/pke/benchmark-full 256 bench | tee bench-256.log
