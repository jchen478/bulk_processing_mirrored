#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "group.h"

using namespace std;

__global__ void group(int *ifiber, int *ncnt, int *ncpf, int *nsegpt,
	int *clist, int *status, int *lead_clist, int *nc, int *maxConpt,
	int *clist_pos, int *maxGrpt, int *groupId, int *num_groups,
	int *total_contact_no_joints){

	int  nseg = *nsegpt;
	int maxCon = *maxConpt;
	int maxGr = *maxGrpt;

	int pos, fib, N, ii, jj, ic;
	int segold1, segold2, segnew1, segnew2;
	bool add;
	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	if (ncpf[tid] != 0 && status[tid] == 0){

		pos = atomicAdd(num_groups, 1); 
		groupId[pos] = tid;

		// boss fiber compiles list
		// find total contacts in group and create list of contacts		
		N = ncpf[tid];

		for (ii = 0; ii < N; ii++){
			add = true;
			for (ic = 0; ic < ncnt[tid]; ic++){

				segnew1 = clist[tid*maxCon + ii];
				segold1 = clist[tid*maxCon + ic];

				// if at joints, don't add to group
				if (segnew1 == segold1 + 1 && segnew1%nseg != 0){
					//printf("1 tid segii segic %4d %4d %4d\n", tid, segnew1, segold1);
					add = false;
					break;
				}
				if (segold1 == segnew1 + 1 && segold1%nseg != 0){
					//printf("2 tid segii segic %4d %4d %4d\n", tid, segnew1, segold1);
					add = false;
					break;
				}
			}
			if (add){
				ifiber[tid * 2 * maxGr + 0 * maxGr + ncnt[tid]] = tid;
				ifiber[tid * 2 * maxGr + 1 * maxGr + ncnt[tid]] = clist[tid*maxCon + ii];
				clist_pos[tid*maxGr + ncnt[tid]] = ii;
				ncnt[tid]++;
				atomicAdd(total_contact_no_joints, 1);
			}		
		}
		for (ii = 0; ii < nc[tid]; ii++){
			fib = lead_clist[tid*maxGr + ii];
			segnew1 = fib;
			N = ncpf[fib];
			for (jj = 0; jj < N; jj++){
				add = true;
				if (fib < clist[fib*maxCon + jj]){
					segnew2 = clist[fib*maxCon + jj];
					for (ic = 0; ic < ncnt[tid]; ic++){
						segold1 = ifiber[tid * 2 * maxGr + 0 * maxGr + ic];
						segold2 = ifiber[tid * 2 * maxGr + 1 * maxGr + ic];
						if (segold2 == segnew2){
							if (segnew1 == segold1 + 1 && segnew1%nseg != 0){
								//printf("9 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
							if (segold1 == segnew1 + 1 && segold1%nseg != 0){
								//printf("10 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
						}
						else if (segold1 == segnew1){
							if (segnew2 == segold2 + 1 && segnew2%nseg != 0){
								//printf("3 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2); 
								add = false;
								break;
							}
							if (segold2 == segnew2 + 1 && segold2%nseg != 0){
								//printf("4 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
						}
						else if (segold1 == segnew2){
							if (segnew1 == segold2 + 1 && segnew1%nseg != 0){
								//printf("5 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
							if (segold2 == segnew1 + 1 && segold2%nseg != 0){
								//printf("6 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
						}
						else if (segold2 == segnew1){
							if (segnew2 == segold1 + 1 && segnew2%nseg != 0){
								//printf("7 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
							if (segold1 == segnew2 + 1 && segold1%nseg != 0){
								//printf("8 tid: %4d old (%4d, %4d) new (%4d, %4d) \n", tid, segold1, segold2, segnew1, segnew2);
								add = false;
								break;
							}
						}
					}
					if (add){
						ifiber[tid * 2 * maxGr + 0 * maxGr + ncnt[tid]] = fib;
						ifiber[tid * 2 * maxGr + 1 * maxGr + ncnt[tid]] = clist[fib*maxCon + jj];
						clist_pos[tid*maxGr + ncnt[tid]] = jj;
						ncnt[tid]++;
						atomicAdd(total_contact_no_joints, 1);
					}
				}
			}
		}
		if (ncnt[tid] >= maxGr){
			printf("error: increase maxGr ncnt[%8d] = %3d\n", tid, ncnt[tid]);
		}
	}
}