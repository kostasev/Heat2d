heat2d: cuda_Heat2D.o cuda_funcs.o
	nvcc -o heat2d cuda_Heat2D.o cuda_funcs.o

cuda_Heat2D.o: cuda_Heat2D.c
	gcc -c cuda_Heat2D.c

cuda_funcs.o: cuda_funcs.cu
	nvcc -c cuda_funcs.cu

clean:
	rm -f *.o heat2d