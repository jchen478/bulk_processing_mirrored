//  calcI.cpp
//  INSSimple 
//
//  Created by Jing-Yao Chen on 6/27/18.
//  Copyright (c) 2018 Jing-Yao Chen. All rights reserved.
//

#include "calcI.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

using namespace std;

float calcI(float *bin, int nxbin, int nybin, int nzbin, 
		float volfrac, float Vcell, float Vseg){
    
	int ind;
	float I, phi;
	I = 0.0;
	for (int ii = 0; ii < nxbin; ii++){
		for (int jj = 0; jj < nybin; jj++){
			for (int kk = 0; kk < nzbin; kk++){
				ind = ii + jj*nxbin + kk*nxbin*nybin;
				phi = bin[ind] * Vseg / Vcell;
				if (phi > 1.0){
					printf("local concentratione exceeds one\n");
					printf("increase bin dimensions\n");
					return -1;
				}
				I += (phi - volfrac)*(phi - volfrac);
			}
		}
	}
	return I;
}
