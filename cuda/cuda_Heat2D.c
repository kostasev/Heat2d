/* =============================================
 * Kostas Evangelou
 * cuda_heat2D.c
 * Parallel Systems
 * CUDA Optimization of heat2D problem
 * =============================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "funcs.h"

#define  NXPROB        80                /* x dimension of problem grid */
#define  NYPROB        64                /* y dimension of problem grid */
#define  STEPS         500               /* number of time steps */


int main(int argc, char *argv[]) {
	float time_exec;

	time_exec = calculate_Heat(NXPROB,NYPROB);
	printf("Time elapsed: %f seconds\n",time_exec);
	return 0;
}