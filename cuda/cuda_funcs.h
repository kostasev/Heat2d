#ifndef CUDA_FUNCS
#define CUDA_FUNCS
#include <stdlib.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

void update(int y, float *u1, float *u2);

#endif