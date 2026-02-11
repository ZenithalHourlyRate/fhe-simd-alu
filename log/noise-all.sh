#!/bin/bash
OMP_NUM_THREADS=16 ./bin/examples/pke/noise-test 64 bench | tee noise-64.log
OMP_NUM_THREADS=16 ./bin/examples/pke/noise-test 32 bench | tee noise-32.log
OMP_NUM_THREADS=16 ./bin/examples/pke/noise-test 16 bench | tee noise-16.log
OMP_NUM_THREADS=16 ./bin/examples/pke/noise-test 8 bench | tee noise-8.log
OMP_NUM_THREADS=16 ./bin/examples/pke/noise-test 128 bench | tee noise-128.log
OMP_NUM_THREADS=16 ./bin/examples/pke/noise-test 256 bench | tee noise-256.log
