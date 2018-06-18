#ifndef __AggChar__errorCheck__
#define __AggChar__errorCheck__

#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>

#define gpuErrchk(ans) {gpuAssert((ans),__FILE__,__LINE__);}
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true){

	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) {
			getchar();
			exit(code);
		}
	}
}

#endif /* defined(__AggChar__errorCheck__) */