//
//  diffusion.cpp
//
//  Created by Jing-Yao Chen on 10/25/17.
//  Copyright (c) 2017 Jing-Yao Chen. All rights reserved.
//

#include <cstring>
#include <stdio.h>
#include <stdlib.h> /* malloc, calloc, free, exit */
#include <cmath>

using namespace std;

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *center_mass, *MSD, *msdin;

	// open files //
	Parameters = fopen("Parameters.in", "r");
	center_mass = fopen("center_mass.txt", "rb");
	msdin = fopen("diffusion.in", "r");
	MSD = fopen("MSD.txt", "w");

	int nfib, nseg, config_write, idum;
	float dt, strain, sidex, sidey, sidez, start_MSD, dum;

	fscanf(msdin, "%f", &start_MSD);
	fclose(msdin);

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &dum); // rps
	fscanf(Parameters, "%*[^\n]%f", &dum); // kb
	fscanf(Parameters, "%*[^\n]%f", &dum); // mu_stat
	fscanf(Parameters, "%*[^\n]%f", &dum); // mu_kin
	fscanf(Parameters, "%*[^\n]%f", &dum); // contact_cutoff
	fscanf(Parameters, "%*[^\n]%f", &dum); // rep_cutoff
	fscanf(Parameters, "%*[^\n]%f", &dum); // overlap
	fscanf(Parameters, "%*[^\n]%f", &dt); 
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &dum); // fraction_rp
	fscanf(Parameters, "%*[^\n]%d", &config_write);
	fclose(Parameters); 
	
	// memory allocation
	float *rcmx, *rcmy, *rcmz;
	float *rcmx0, *rcmy0, *rcmz0;
	float *rcmxp, *rcmyp, *rcmzp;
	float *MSDxx, *MSDyy, *MSDzz, *r; 
	float diffx, diffy, diffz, sumx, sumy, sumz, sumr; 
	float *crossxx, *crossyy, *crosszz;
	float dx, dy, dz;
	int step, nConfig, start_MSD_step, m;
	
	nConfig = int(strain / dt / float(config_write)) + 1;
	rcmx = (float*)malloc(nfib*sizeof(float));
	rcmy = (float*)malloc(nfib*sizeof(float));
	rcmz = (float*)malloc(nfib*sizeof(float));
	rcmx0 = (float*)malloc(nfib*sizeof(float));
	rcmy0 = (float*)malloc(nfib*sizeof(float));
	rcmz0 = (float*)malloc(nfib*sizeof(float));
	rcmxp = (float*)malloc(nfib*sizeof(float));
	rcmyp = (float*)malloc(nfib*sizeof(float));
	rcmzp = (float*)malloc(nfib*sizeof(float));
	r = (float*)malloc(nConfig*sizeof(float));
	MSDxx = (float*)calloc(nConfig,sizeof(float));
	MSDyy = (float*)calloc(nConfig,sizeof(float));
	MSDzz = (float*)calloc(nConfig,sizeof(float));
	crossxx = (float*)calloc(nfib,sizeof(float));
	crossyy = (float*)calloc(nfib,sizeof(float));
	crosszz = (float*)calloc(nfib,sizeof(float));
	
	// based on writing frequency, find the step to start calculate MSD
	start_MSD_step = int(start_MSD / dt / (float(config_write))); 

	// read through configuration vs. strain
	for (step = 0; step < nConfig; step++){

		// read center of mass at step
		fread(&dum, sizeof(float), 1, center_mass); 
		fread(rcmx, sizeof(float), nfib, center_mass); 
		fread(rcmy, sizeof(float), nfib, center_mass);
		fread(rcmz, sizeof(float), nfib, center_mass);

		// store initial configuration to calculate MSD
		if (step == start_MSD_step){
			// initial configuration
			memcpy(rcmx0,rcmx,nfib*sizeof(float));
			memcpy(rcmy0,rcmy,nfib*sizeof(float));
			memcpy(rcmz0,rcmz,nfib*sizeof(float));
			// configuration at previous step
			memcpy(rcmxp,rcmx,nfib*sizeof(float));
			memcpy(rcmyp,rcmy,nfib*sizeof(float));
			memcpy(rcmzp,rcmz,nfib*sizeof(float));
			continue;
		}
		// calculate difference with previous step 
		// to update number of times 
		// crossing periodic boundary
		for (m = 0; m < nfib; m++){
			diffx = (rcmx[m] - rcmxp[m]);
			diffy = (rcmy[m] - rcmyp[m]);
			diffz = (rcmz[m] - rcmzp[m]);
			crossxx[m] -= roundf(diffx/sidex);		
			crossyy[m] -= roundf(diffy/sidey);		
			crosszz[m] -= roundf(diffz/sidez);	
		}
		// update configuration at previous step 
		memcpy(rcmxp,rcmx,nfib*sizeof(float));
		memcpy(rcmyp,rcmy,nfib*sizeof(float));
		memcpy(rcmzp,rcmz,nfib*sizeof(float));
		// average msd across fibers
		sumr = 0.0; sumx = 0.0; sumy = 0.0; sumz = 0.0; 
		for (m = 0; m < nfib; m++){
			// correct distance due to periodic boundary
			rcmx[m] += sidex*crossxx[m];			
			rcmy[m] += sidey*crossyy[m];			
			rcmz[m] += sidez*crosszz[m];				
			// update sum
			dx = rcmx[m]-rcmx0[m];
			dy = rcmy[m]-rcmy0[m];
			dz = rcmz[m]-rcmz0[m];
			sumx += dx*dx;
			sumy += dy*dy;
			sumz += dz*dz;
			sumr += dx*dx+dy*dy+dz*dz;
		}
		MSDxx[step] = sumx / float(2*nfib);
		MSDyy[step] = sumy / float(2*nfib);
		MSDzz[step] = sumz / float(2*nfib);
		r[step] = sumr / float(2*nfib);
	}
	// printf results to file
	for (step = 0; step < nConfig; step++){
		
		if (step >= start_MSD_step){
			fprintf(MSD,"%6.2f %15.4f %15.4f %15.4f %15.4f\n", 
				float(step*config_write)*dt,
				MSDxx[step], MSDyy[step],MSDzz[step], 
				r[step]);
		}
	}
	return 0;
}
