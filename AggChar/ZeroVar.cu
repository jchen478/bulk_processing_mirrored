#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "ZeroVar.h"

using namespace std;

__global__ void ZeroVar(int *array){

	int tid = threadIdx.x + blockIdx.x*blockDim.x;
	array[tid] = 0; 
}
