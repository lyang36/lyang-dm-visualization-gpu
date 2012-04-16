

all: kernelpp.o skymap.o setparams.o VL2_debug.o main.o
	@nvcc --link kernelpp.o skymap.o setparams.o VL2_debug.o main.o -o VL2_gpu
	@echo success

kernelpp.o:
	@nvcc -c kernelpp.cu 

skymap.o:
	@nvcc -c skymap.cu

setparams.o:
	@nvcc -c setparams.cpp

VL2_debug.o:
	@nvcc -c VL2_debug.cpp

main.o:
	@nvcc -c main.cpp

clean:
	@rm *.o
	@rm VL2_gpu

