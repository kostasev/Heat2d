/* =============================================
 * Kostas Evangelou
 * cuda_heat2D.c
 * Parallel Systems
 * CUDA Optimization of heat2D problem
 * =============================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "timestamp.h"

#define NXPROB 2000
#define NYPROB 3000
#define STEPS  500
#define CHECK_INTERVAL 20
#define THREADS_PER_ROW 32
#define THREADS_PER_BLOCK (THREADS_PER_ROW*THREADS_PER_ROW)
#define ROW_BLOCKS (((NXPROB-2)/THREADS_PER_ROW)+(((NXPROB-2)%THREADS_PER_ROW==0)?(0):(1)))
#define COL_BLOCKS (((NYPROB-2)/THREADS_PER_ROW)+(((NYPROB-2)%THREADS_PER_ROW==0)?(0):(1)))
#define NUMBER_OF_BLOCKS (ROW_BLOCKS*COL_BLOCKS)
#define PARMS_CX 0.1f
#define PARMS_CY 0.1f

// function to get the index for a 2d array that is saved as an 1d array
#define index2D(i,j,col) ((i)*(col)+(j))

void inidat(int, int, float*), prtdat(int, int, float*, char*);

__global__ void heat(float *old_grid, float *new_grid)
{
        int my_x = (blockIdx.x * blockDim.x) + threadIdx.x;
        int my_y = (blockIdx.y * blockDim.y) + threadIdx.y;
        // Increment each coordinate by 1. This is used in an attempt to save resources: 
        // The first row and column stay always 0, so there is no point in having threads there just sitting.
        my_x += 1;
        my_y += 1;

        // update cell
        if (my_x < NXPROB-1 && my_y < NYPROB-1)
        {
                new_grid[index2D(my_x,my_y,NYPROB)] = old_grid[index2D(my_x,my_y,NYPROB)] +
                                                        PARMS_CX * (old_grid[index2D(my_x+1,my_y,NYPROB)] +
                                                        old_grid[index2D(my_x-1,my_y,NYPROB)] -
                                                        2.0f * old_grid[index2D(my_x,my_y,NYPROB)]) +
                                                        PARMS_CY * (old_grid[index2D(my_x,my_y+1,NYPROB)] +
                                                        old_grid[index2D(my_x,my_y-1,NYPROB)] -
                                                        2.0f * old_grid[index2D(my_x,my_y,NYPROB)]);
        }
}

int main(void)
{
	int i, old = 0;
	float *grid, msecs;
	if ((grid = (float *)malloc(NXPROB*NYPROB*sizeof(float))) == NULL)
	{
		perror("malloc for grid");
		return -1;
	}

	// allocate space in device for the 2 arrays
	float *dvc_grid0, *dvc_grid1, *dvc_grids[2];
	cudaMalloc((void**)&dvc_grid0, NXPROB*NYPROB*sizeof(float));
	cudaMalloc((void**)&dvc_grid1, NXPROB*NYPROB*sizeof(float));

	dvc_grids[0] = dvc_grid0;
	dvc_grids[1] = dvc_grid1;

	// initialize starting grid
	inidat(NXPROB, NYPROB, grid);

	// copy the initialized data to device memory
	cudaMemcpy(dvc_grid0, grid, NXPROB*NYPROB*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dvc_grid1, grid, NXPROB*NYPROB*sizeof(float), cudaMemcpyHostToDevice);

	dim3 numofthreads(THREADS_PER_ROW, THREADS_PER_ROW);
	dim3 numofblocks(ROW_BLOCKS, COL_BLOCKS);

	timestamp t_start = getTimestamp();

	for (i = 0 ; i < STEPS ; ++i)
	{
		// launch kernel on GPU
		heat<<<numofblocks, numofthreads>>>(dvc_grids[old], dvc_grids[1-old]);

		// iteration finished, set current grid as old
		old = 1-old;

	}
	msecs = getElapsedtime(t_start);

	// get the results from GPU ("old" grid has the final values, after the last iteration)
	cudaMemcpy(grid, dvc_grids[old], NXPROB*NYPROB*sizeof(float), cudaMemcpyDeviceToHost);

	// print info
	printf("Elapsed time: %.3f %ssecs\n", (msecs/1000 > 1.0 ? msecs/1000 : msecs), (msecs/1000 > 1.0 ? "" : "m"));

	// free resources and exit
	cudaFree(dvc_grid0);
	cudaFree(dvc_grid1);
	cudaFree(dvc_grids);
	free(grid);

	return 0;
}

/*****************************************************************************
 *  *  subroutine inidat
 *   *****************************************************************************/
void inidat(int nx, int ny, float *u) {
int ix, iy;

for (ix = 0; ix <= nx-1; ix++)
  for (iy = 0; iy <= ny-1; iy++)
     *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
}


