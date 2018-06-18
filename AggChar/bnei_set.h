#ifndef __AggChar__bnei_set__
#define __AggChar__bnei_set__

#include <stdio.h>
#include <cuda_runtime_api.h>
#include "bnei_set.h"

__global__ void bnei_set(int *bnei, int *nxbinpt, int *nybinpt, int *nzbinpt,
	int *bdimxpt, int *bdimypt, int *bdimzpt);

#endif /* defined(__AggChar__bnei_set__) */