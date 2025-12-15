include ../Makefile.in

# override for debugging
#COMPOPT=-g -O0 -DDEBUG
#CUDAOPT=-g -O0 -DDEBUG 

# CUDA VERSION
# diffusive
salem-srb-diff_GPU: buildAGI.o diften.o inputParser.o srb_diff_CUDA.o nmr.o varcovar.o d2k0.o muEuler.o Cartesian.o d2k0Gradients.o orderParameter.o getNH.o rotomat.o
	$(NVCC) srb_diff_CUDA.o buildAGI.o diften.o inputParser.o nmr.o varcovar.o d2k0.o muEuler.o Cartesian.o d2k0Gradients.o orderParameter.o getNH.o rotomat.o -D__GPU__ -DADD_ $(BRANCH) $(LL) -lmagma -llapack -lblas -lcblas_LINUX -lzmat -lgfortran -lcublas -lcusparse -lcudart -lgomp -lcqp++ -o salem-srb-diff_GPU -Xcompiler -static-libstdc++ -lopenblas 

# CPU-GPU CODE
srb_diff_CUDA.o: srb_diff_CUDA.cu GPUkernels.cu
	$(NVCC) $(CUDAOPT) $(II) $(BRANCH) -D__GPU__ $(CUDA_ARCH) -c srb_diff_CUDA.cu

# CPU ONLY CODE
buildAGI.o: buildAGI.cpp buildAGI.h
	$(CPP) $(COMPOPT) $(II) -c buildAGI.cpp

diften.o: diften.cpp diften.h
	$(CPP) $(COMPOPT) $(II) -c diften.cpp

inputParser.o: inputParser.cpp inputParser.h
	$(CPP) $(COMPOPT) $(II) -c inputParser.cpp

nmr.o: nmr.cpp nmr.h
	$(CPP) $(COMPOPT) $(II) -c nmr.cpp

varcovar.o: varcovar.cpp ./ADEL/Operations/Cartesian.h ./ADEL/Operations/CartesianFirstDerivative.h ./ADEL/Operations.h
	$(CPP) $(COMPOPT) $(II) -c varcovar.cpp

d2k0.o: d2k0.cpp d2k0.h
	$(CPP) $(COMPOPT) $(II) -c d2k0.cpp

muEuler.o: muEuler.cpp muEuler.h
	$(CPP) $(COMPOPT) $(II) -c muEuler.cpp

Cartesian.o: Cartesian.cpp
	$(CPP) $(COMPOPT) $(II) -c Cartesian.cpp

d2k0Gradients.o: d2k0Gradients.cpp d2k0Gradients.h d2k0.h muEuler.h
	$(CPP) $(COMPOPT) $(II) -c d2k0Gradients.cpp

orderParameter.o: orderParameter.cpp orderParameter.h d2k0Gradients.h d2k0.h muEuler.h
	$(CPP) $(COMPOPT) $(II) -c orderParameter.cpp

getNH.o: getNH.cpp getNH.h
	$(CPP) $(COMPOPT) $(II) -c getNH.cpp

rotomat.o: rotomat.cpp rotomat.h
	$(CPP) $(COMPOPT) $(II) -c rotomat.cpp

clean:
	rm -rf ./*.o ./salem-srb-diff_GPU

