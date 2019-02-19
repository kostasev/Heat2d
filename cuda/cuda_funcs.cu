
#include "cuda_funcs.h"

__global__ void update(int y, float *u1, float *u2)
{	
	int ix = 0;
	int iy = 0;
	
	struct Parms { 
		float cx;
		float cy;
	} parms = {0.1, 0.1};
	
	/* Get coordinators */
	ix = blockIdx.x*blockDim.x+threadIdx.x;
	iy = blockIdx.y*blockDim.y+threadIdx.y;
  	
	*(u2 + ix * y + iy) = *(u1 + ix * y + iy) +
		parms.cx * (*(u1 + (ix + 1) * y + iy) +
		*(u1 + (ix - 1) * y + iy) -
		2.0 * *(u1 + ix * y + iy)) +
		parms.cy * (*(u1 + ix * y + iy + 1) +
		*(u1 + ix * y + iy - 1) -
		2.0 * *(u1 + ix * y + iy));
	
	__syncthreads();
}