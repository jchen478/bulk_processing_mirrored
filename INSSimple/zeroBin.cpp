//  zeroBin.cpp
//  INSSimple 
//
//  Created by Jing-Yao Chen on 6/27/18.
//  Copyright (c) 2018 Jing-Yao Chen. All rights reserved.
//

#include "zeroBin.h"
#include <iostream>
using namespace std;

int zeroBin(float *bin, int nxbin, int nybin, int nzbin){
    
	int ind;
	for (int ii = 0; ii < nxbin; ii++){
		for (int jj = 0; jj < nybin; jj++){
			for (int kk = 0; kk < nzbin; kk++){
				ind = ii + jj*nxbin + kk*nxbin*nybin;
				bin[ind] = 0.0;
			}
		}
	}
	return 0;
}
