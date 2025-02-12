#! /bin/bash

make clean
make nosimd || exit 1
PROGRAM="taskset -c 0 ./main_nosimd"

# Iteration range
N_SEQ=$(seq 3 10)

N_=""
for N in $N_SEQ; do
	N_+=" $((2**N))"
done

# Variants
VARIANT_="20 21 22 23 24 25"

source benchmark_base.sh
