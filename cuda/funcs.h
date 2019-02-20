
#ifndef FUNCS_H
#define FUNCS_H

#include <cuda.h>

void  inidat0(int nx, int ny, float *u);
void inidat1(int nx, int ny, float *u);
float calculate_Heat(int NX, int NY);

#endif