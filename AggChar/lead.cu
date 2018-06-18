#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "lead.h"

using namespace std;


__global__ void lead(int *ncpf, int *clist, int *status,
	int *nc, int *lead_clist, int *maxConpt, int *maxGrpt){

	int maxCon = *maxConpt;
	int maxGr = *maxGrpt;

	int tid = threadIdx.x + blockIdx.x*blockDim.x;

	int pos, fiber, locfiber, i, j, k; 
	bool add;

	pos = 0;
	add = false;

	for (i = 0; i < ncpf[tid]; i++){
		fiber = clist[tid*maxCon + i];
		if (fiber < tid){
			status[tid] = 1;
			return;
		}
		lead_clist[maxGr*tid + nc[tid]] = clist[tid*maxCon + i];
		nc[tid]++;
	}
	while (pos != nc[tid]){
		fiber = lead_clist[maxGr*tid + pos];
		for (j = 0; j < ncpf[fiber]; j++){
			locfiber = clist[fiber*maxCon + j];
			if (locfiber < tid){
				status[tid] = 1;
				return;
			}
			if (locfiber == tid){
				continue;
			}
			add = true;
			for (k = 0; k < nc[tid]; k++){
				if (locfiber == lead_clist[maxGr*tid + k]){
					add = false;
					break;
				}
			}
			if (add){
				lead_clist[maxGr*tid + nc[tid]] = locfiber;
				nc[tid]++;
				if (nc[tid] >= maxGr){
					printf("error in lead: increase maxGr\n"); 
				}
			}
		}
		pos++;
	}

}