#ifndef __AggChar__group__
#define __AggChar__group__

#include <cuda_runtime_api.h>

__global__ void group(int *ifiber, int *ncnt, int *ncpf, int *nsegpt,
	int *clist, int *status, int *lead_clist, int *nc, int *maxConpt,
	int *clist_pos, int *maxGrpt, int *groupId, int *num_groups,
	int *total_contact_no_joints);

#endif /* defined(__AggChar__group__) */