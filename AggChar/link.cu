#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "link.h"

using namespace std;

__global__ void link(int *bin, int *list, int *bnei, int *bnum,
	int *potCon, int *potConSize, int *npcnpt, 
	int *nsegpt, int *nxbinpt, int *nybinpt, int *nzbinpt){

	int npcn = *npcnpt;
	int nseg = *nsegpt;
	int nxbin = *nxbinpt;
	int nybin = *nybinpt;
	int nzbin = *nzbinpt;

	int mi = blockIdx.x; 
	int tid = threadIdx.x; 
	int miBin = bnum[mi]; 
	int nextBin = bnei[miBin + tid*nxbin*nybin*nzbin]; 
	int nextTot = bin[nextBin]; 

	int nj, i, pos; 

	if (tid == 0){
		potConSize[mi] = 0; 
	}
	__syncthreads(); 

	for (i = 0; i < nextTot; i++){
		nj = list[nextBin + i*nxbin*nybin*nzbin];
		// skip when segments are connected
		if (mi >= nj){
			continue;
		}
		if ((mi - (mi / nseg)*nseg) != 0 && nj == mi - 1){
			continue;
		}
		if ((nj - (nj / nseg)*nseg) != 0 && mi == nj - 1){
			continue;
		}
		if ((mi - (mi / nseg)*nseg) != nseg - 1 && nj == mi + 1){
			continue;
		}
		if ((nj - (nj / nseg)*nseg) != nseg - 1 && mi == nj + 1){
			continue;
		}
		//printf("link mi nj %4d %4d\n", mi, nj); 
		pos = atomicAdd(potConSize + mi, 1); 
		if (pos >= npcn - 1) printf("allocate more space for potCon pos %4d\n", pos); 
		potCon[mi * npcn + pos] = nj; 
	}
}
