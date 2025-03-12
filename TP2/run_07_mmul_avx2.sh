#! /bin/bash

make clean
make avx2 || exit 1
PROGRAM="./main_avx2"

# Iteration range
N_SEQ=$(seq 4 10)

N_=""
for N in $N_SEQ; do
	N_+=" $((2**N))"
	#N_+=" $((512*N))"
done

# Variants
VARIANT_="35"

source benchmark_base.sh
