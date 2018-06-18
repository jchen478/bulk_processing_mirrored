#ifndef __AggChar__lead__
#define __AggChar__lead__

#include <cuda_runtime_api.h>

__global__ void lead(int *ncpf, int *clist, int *status,
	int *nc, int *lead_clist, int *maxConpt, int *maxGrpt);

#endif /* defined(__AggChar__lead__) */