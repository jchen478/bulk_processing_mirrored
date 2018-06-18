#ifndef __AggChar__contact__
#define __AggChar__contact__

#include <cuda_runtime_api.h>

__global__ void contact(int *total_contact, int *total_overlap,
	int *potCon, int *potConSize, int *npcnpt,
	float *rx, float *ry, float *rz,
	float *px, float *py, float *pz,
	float *sidexpt, float *sideypt, float *sidezpt,
	float *delta_rxpt, float *rppt, float *over_cutpt,
	float *contact_cutoffpt, float *rep_cutoffpt,
	int *maxConpt, int *ncpf, int *clist);

#endif /* defined(__AggChar__contact__) */