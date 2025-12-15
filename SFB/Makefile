include ../conf/make.COMPILERS.CPU-GPU
include ../conf/make.in

# Override for debugging
#COMPOPT=-g -pg -O0 -std=c++11 -DADD_ -DDEBUG
#CUDAOPT=-g -pg -O0 -std=c++11 -DADD_ -gencode arch=compute_20,code=sm_20 -DDEBUG

all: sfb sfb_RK

sfb: input.o sfb.o
	$(NVCC) $(CUDAOPT) input.o sfb.o -DADD_ $(LL) $(II) -lcublas -lcudart -lcurand -lcusparse -lEulerQuaternion -lmagma -llapack -lblas -lcblas_LINUX -lgomp -o sfb -lopenblas

sfb_RK: input_RK.o sfb_RK.o
	$(NVCC) $(CUDAOPT) input_RK.o sfb_RK.o -DADD_ $(LL) $(II) -lcublas -lcudart -lcurand -lcusparse -lEulerQuaternion -lmagma -llapack -lblas -lcblas_LINUX -lgomp -o sfb_RK -lopenblas

input.o: input.cu input.h
	$(NVCC) $(CUDAOPT) $(CUDA_ARCH) $(II) -c input.cu

sfb.o: sfb.cu GPUkernels.cu sfb.h
	$(NVCC) $(CUDAOPT) $(CUDA_ARCH) $(II) -c sfb.cu

input_RK.o: input_RK.cu input_RK.h
	$(NVCC) $(CUDAOPT) $(CUDA_ARCH) $(II) -c input_RK.cu

sfb_RK.o: sfb_RK.cu GPUkernels_RK.cu sfb_RK.h
	$(NVCC) $(CUDAOPT) $(CUDA_ARCH) $(II) -c sfb_RK.cu

clean:
	rm -rf ./sfb ./sfb_RK ./*.o

