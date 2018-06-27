//
//  calcI.h
//  INSSimple 
//
//  Created by Jing-Yao Chen on 6/27/18.
//  Copyright (c) 2018 Jing-Yao Chen. All rights reserved.
//

#ifndef __INSSimple__calcI__
#define __INSSimple__calcI__

#include <stdio.h>

float calcI(float *bin, int nxbin, int nybin, int nzbin, 
		float volfrac, float Vcell, float Vseg);

#endif /* defined(__INSSimple__calcI__) */
