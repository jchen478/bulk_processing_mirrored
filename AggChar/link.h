#ifndef __AggChar__link__
#define __AggChar__link__

#include <cuda_runtime_api.h>

__global__ void link(int *bin, int *list, int *bnei, int *bnum,
	int *potCon, int *potConSize, int *npcnpt,
	int *nsegpt, int *nxbinpt, int *nybinpt, int *nzbinpt);

#endif /* defined(__AggChar__link__) */