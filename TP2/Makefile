
CXXFLAGS=

#
# Disable loop blocking optimizations by compiler
# We do this since we like to learn about optimization
#
#CXXFLAGS+=-fno-loop-block -fno-loop-strip-mine -fno-loop-interchange

# Generate code specifically for this architecture
CXXFLAGS+=-march=native

# Use newer C++ version
CXXFLAGS+=-std=c++11
CXXFLAGS_SIMD+=-fopenmp-simd
CXXFLAGS_OMP+=-fopenmp
LDFLAGS_OMP+=-fopenmp

# MKL flags
# LDFLAGS+=-lmkl_intel_ilp64 -lmkl_core -lm -lc
# LDFLAGS_NOOMP+=-lmkl_sequential
# LDFLAGS_OMP+=-lmkl_gnu_thread

# AVX2
CXXFLAGS+=-ftree-vectorize -fstrict-aliasing -O3
# CXXFLAGS+= -fno-tree-vectorize -fstrict-aliasing 

all: nosimd simd openmp debug
	@echo "Compiling all targets"


# default is to use NO vectorization
nosimd: main.cpp
	# Try to disable vectorization.
	# Still uses SSE, but not in a vectorized way
	$(CXX) -O3 $(CXXFLAGS) -fno-tree-vectorize -mno-avx main.cpp -c -o main_nosimd.o
	$(CXX) main_nosimd.o $(LDFLAGS) $(LDFLAGS_NOOMP) -o main_nosimd

simd: main.cpp
	$(CXX) -O3 $(CXXFLAGS) $(CXXFLAGS_SIMD) main.cpp -c -o main_simd.o
	$(CXX) main_simd.o $(LDFLAGS) $(LDFLAGS_NOOMP) -o main_simd

openmp: main.cpp
	$(CXX) -O3 $(CXXFLAGS) $(CXXFLAGS_SIMD) $(CXXFLAGS_OMP) main.cpp -c -o main_openmp.o
	$(CXX) main_openmp.o $(LDFLAGS) $(LDFLAGS_OMP) -o main_openmp

avx2: main.cpp
	$(CXX) -O3 -mavx2 $(CXXFLAGS) $(CXXFLAGS_SIMD) $(CXXFLAGS_OMP) main.cpp -c -o main_avx2.o
	$(CXX) main_avx2.o $(LDFLAGS) $(LDFLAGS_OMP) -o main_avx2

opti: main.cpp
	$(CXX) -O3 -mavx2 $(CXXFLAGS) $(CXXFLAGS_SIMD) $(CXXFLAGS_OMP) main.cpp -c -o main_opti.o
	$(CXX) main_opti.o $(LDFLAGS) $(LDFLAGS_OMP) -o main_opti


debug:
	$(CXX) -O0 $(CXXFLAGS) main.cpp -c -o main_debug.o
	$(CXX) main_debug.o $(LDFLAGS) $(LDFLAGS_NOOMP) -o main_debug

clean:
	rm -f main_debug.o
	rm -f main_debug

	rm -f main_nosimd.o
	rm -f main_nosimd

	rm -f main_simd.o
	rm -f main_simd

	rm -f main_openmp.o
	rm -f main_openmp

	rm -f main_avx2.o
	rm -f main_avx2

	rm -f main_opti.o
	rm -f main_opti

