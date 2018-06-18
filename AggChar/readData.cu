/**
* \brief Read input into arrays
*
* Read data into arrays
* First read time stamp then data
* @param file file handle
* @param array array that holds data
* @param size size of array 
*/

#include "cuda_runtime.h"

#include <stdio.h>
#include <stdlib.h>
#include "readData.h"

using namespace std;

void readData(FILE *file, float *array, int size){

	float dum;
	fread(&dum, sizeof(float), 1, file);
	fread(array, sizeof(float), size, file);

}