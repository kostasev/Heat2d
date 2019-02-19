all: heat2d

heat2d: cuda_Heat2D.o cuda_funcs.o funcs.o
	nvcc -o heat2d cuda_Heat2D.o cuda_funcs.o funcs.o

cuda_Heat2D.o: cuda_Heat2D.c
	nvcc  -c cuda_Heat2D.c

funcs.o: funcs.c
	gcc -Wall -c funcs.c

cuda_funcs.o: cuda_funcs.cu
	nvcc  -c cuda_funcs.cu

clean:
	rm -f *.o heat2d