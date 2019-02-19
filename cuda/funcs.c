
#include "funcs.h"

float calculate_Heat(int NX, int NY){
	float	*table_1,*table_2,*table_host;
	int	it,
		size,
		iz = 0;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	/* Calculate size and allocate memory in device for the arrays */
	size = NX*NY*sizeof(float);
	CUDA_SAFE_CALL(cudaMalloc((void**)&table_1,(long) size));
	CUDA_SAFE_CALL(cudaMalloc((void**)&table_2,(long) size));

	/* Allocate memory in host for the array */
	table_host = (float*)malloc(size);
	if (table_host == NULL) {
		printf("Main ERROR: Allocation memory.\n");
		exit(-1);
	}
	
	/* Initialize table_host with zeros and then call inidat*/
	memset(table_host, 0, NX*NY*sizeof(float));
	inidat(NX, NY, table_host);
	
	/* Copy table_1 and table_2 to GPU */
	CUDA_SAFE_CALL(cudaMemcpy(table_1, table_host, size, cudaMemcpyHostToDevice));	
	CUDA_SAFE_CALL(cudaMemcpy(table_2, table_host, size, cudaMemcpyHostToDevice));
	 
	/* Create N blocks of N threads each */
	dim3 NumberOfThreads(NX);			
	dim3 NumberOfBlocks(NY);
	
	/* Start the Clock */
	cudaEventRecord(start);
	/* Go! */
	for (it = 1; it <= STEPS; it++)
	{       
		if ( iz==0 ){
			update<<<NumberOfBlocks,NumberOfThreads>>>(NY, table_1, table_2);
		}
		else {
			update<<<NumberOfBlocks,NumberOfThreads>>>(NY, table_2, table_1);
		}
		
		/* Swap table pointers for next loop */
		iz = 1 - iz;
        	
		/* Sync Cuda Threads */
		CUDA_SAFE_CALL(cudaThreadSynchronize());		
	}
	cudaEventRecord(stop);
	
	/* Copy table with results to table_host from GPU */
	CUDA_SAFE_CALL(cudaMemcpy(table_host, table_2, NX*NY*sizeof(float), cudaMemcpyDeviceToHost));
	
	/* Free Resources */
	CUDA_SAFE_CALL(cudaFree(table_1) );	
	CUDA_SAFE_CALL(cudaFree(table_2) );
	free(table_host);
	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	return milliseconds;
}
