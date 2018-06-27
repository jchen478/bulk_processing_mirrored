//  binning.cpp
//  INSSimple 
//
//  Created by Jing-Yao Chen on 6/27/18.
//  Copyright (c) 2018 Jing-Yao Chen. All rights reserved.
//

#include "binning.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

using namespace std;

int binning(float *bin, int nxbin, int nybin, int nzbin, 
		float dx, float dy, float dz,
		float *rx, float *ry, float *rz, 
		int nfib, int nseg, float sidex, 
		float sidey, float sidez){
    
	int mi, xloc, yloc, zloc;
	float rxmi, rymi, rzmi;
	for (int m = 0; m < nfib; m++){
		for (int i = 0; i < nseg; i++){
			mi = m*nseg+i; 
			rxmi = rx[mi];
			rymi = ry[mi];
			rzmi = rz[mi];					
			// shift coordinate of segments to stay in the box
			rxmi -= sidex*roundf(rxmi / sidex);
			rymi -= sidey*roundf(rymi / sidey);
			rzmi -= sidez*roundf(rzmi / sidez);
			// shift coordinates to all be positive
			rxmi += sidex / 2.0;
			rymi += sidey / 2.0;
			rzmi += sidez / 2.0;
			// find associated binning index
			xloc = int(rxmi / dx);
			yloc = int(rymi / dy);
			zloc = int(rzmi / dz);
			// if right at edge of box
			// count as the smaller box index next to it
			if (xloc == nxbin) xloc--; 
			if (yloc == nybin) yloc--;
			if (zloc == nzbin) zloc--;
			// check for out of bound indices
			if (xloc >= nxbin || yloc >= nybin || zloc >= nzbin){
				printf("%8d %8d %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %6d %6d %6d\n",
					m + 1, i + 1, rx[m*nseg + i], ry[m*nseg + i], rz[m*nseg + i], rxmi, rymi, rzmi, xloc, yloc, zloc);						
				getchar(); 
			}
			if (xloc < 0 || yloc < 0 || zloc < 0){
				printf("%8d %8d %13.6f %13.6f %13.6f %13.6f %13.6f %13.6f %6d %6d %6d\n",
					m + 1, i + 1, rx[m*nseg + i], ry[m*nseg + i], rz[m*nseg + i], rxmi, rymi, rzmi, xloc, yloc, zloc);
				getchar();
			}
			// increase bin count
			bin[xloc + yloc*nxbin + zloc*nxbin*nybin] += 1.0; 
		}
	}
	return 0;
}
