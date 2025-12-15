include ../conf/make.COMPILERS.CPU-GPU
include ../conf/make.in

# Override for debugging
#COMPOPT=-g -pg -O0 -std=c++11 -DADD_ -DDEBUG
#CUDAOPT=-g -pg -O0 -std=c++11 -DADD_ -gencode arch=compute_20,code=sm_20 -DDEBUG

sfb: file_reader.o matrix_operations.o ztocart.o
	$(NVCC) $(CUDAOPT) file_reader.o matrix_operations.o ztocart.o $(LL) $(II) -lcublas -lcusparse -lcudart -lEulerQuaternion -lmagma -llapack -lgomp -lzmat -o z2cart


file_reader.o: file_reader.cpp file_reader.h
	$(CC) $(COMPOPT) $(II) -c file_reader.cpp

matrix_operations.o: matrix_operations.cpp matrix_operations.h
	$(NVCC) $(CUDAOPT) $(CUDA_ARCH) $(II) -c matrix_operations.cpp

ztocart.o: ztocart.cu
	$(NVCC) $(CUDAOPT) $(CUDA_ARCH) $(II) -c ztocart.cu

clean:
	rm -rf ./z2cart ./*.o

