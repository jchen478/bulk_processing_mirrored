//
//  binning.h
//  INSSimple 
//
//  Created by Jing-Yao Chen on 6/27/18.
//  Copyright (c) 2018 Jing-Yao Chen. All rights reserved.
//

#ifndef __INSSimple__binning__
#define __INSSimple__binning__

#include <stdio.h>

int binning(float *bin, int nxbin, int nybin, int nzbin, 
		float dx, float dy, float dz,
		float *rx, float *ry, float *rz, 
		int nfib, int nseg, float sidex, 
		float sidey, float sidez);

#endif /* defined(__INSSimple__binning__) */
