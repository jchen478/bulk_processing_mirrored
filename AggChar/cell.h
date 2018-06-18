#ifndef __AggChar__cell__
#define __AggChar__cell__

#include <cuda_runtime_api.h>

__global__ void cell(int *bin, int *list, int *bnum,
	float *rx, float *ry, float *rz,
	int *nxbinpt, int *nybinpt, int *nzbinpt,
	int *maxBinpt, float *sidexpt, float *sideypt, float *sidezpt,
	float *delta_rxpt, float *dxpt, float *dypt, float *dzpt);

#endif /* defined(__AggChar__cell__) */