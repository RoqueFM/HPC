#! /bin/bash

make clean
make nosimd || exit 1
PROGRAM="taskset -c 0 ./main_nosimd"

# Iteration range
N_SEQ=$(seq 4 11)

N_=""
for N in $N_SEQ; do
	N_+=" $((2**N))"
	#N_+=" $((512*N))"
done

# Variants
VARIANT_="33 32 31 22"

#DEBUG_EXEC=false
#DEBUG_EXEC=true

source benchmark_base.sh

