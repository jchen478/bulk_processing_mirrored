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
	FILE *Parameters, *center_mass, *MSD;

	// open files //
	Parameters = fopen("Parameters.in", "r");
	center_mass = fopen("center_mass.txt", "rb");
	MSD = fopen("MSD.txt", "w");

	int nfib, nseg, config_write, idum;
	float dt, strain, sidex, sidey, sidez, start_MSD, dum;

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
	fscanf(Parameters, "%*[^\n]%d", &idum); // contact_write
	fscanf(Parameters, "%*[^\n]%f", &dum); // fstarr
	fscanf(Parameters, "%*[^\n]%f", &dum); // fact
	fscanf(Parameters, "%*[^\n]%f", &dum); // Astart
	fscanf(Parameters, "%*[^\n]%f", &dum); // decatt
	fscanf(Parameters, "%*[^\n]%f", &dum); // delta_rx
	fscanf(Parameters, "%*[^\n]%f", &dum); // duidj
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum); // fac
	fscanf(Parameters, "%*[^\n]%f", &dum); // elf
	fscanf(Parameters, "%*[^\n]%f", &dum); // di
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum); // bdimi
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum); // fiberPerBlock
	fscanf(Parameters, "%*[^\n]%d", &idum); // maxCon
	fscanf(Parameters, "%*[^\n]%d", &idum); // maxGr
	fscanf(Parameters, "%*[^\n]%d", &idum); // maxBin
	fscanf(Parameters, "%*[^\n]%d", &idum); // nfibGrid
	fscanf(Parameters, "%*[^\n]%d", &idum); // nfibBlock
	fscanf(Parameters, "%*[^\n]%d", &idum); // blasGrid
	fscanf(Parameters, "%*[^\n]%d", &idum); // blasBlock
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_stress
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_binToAscii
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_pair
	fscanf(Parameters, "%*[^\n]%d", &idum); // stress_write
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_contacts
	fscanf(Parameters, "%*[^\n]%f", &dum); // concFac
	fscanf(Parameters, "%*[^\n]%f", &dum); // expFac
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_conc
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_exp
	fscanf(Parameters, "%*[^\n]%f", &dum); // eq_conc
	fscanf(Parameters, "%*[^\n]%f", &dum); // eq_exp
	fscanf(Parameters, "%*[^\n]%d", &idum); // n_conc
	fscanf(Parameters, "%*[^\n]%d", &idum); // n_exp
	fscanf(Parameters, "%*[^\n]%f", &dum); // di_entropy
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum); // nxbin_dist
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum); // concPHase
	fscanf(Parameters, "%*[^\n]%d", &idum); // INSPhase
	fscanf(Parameters, "%*[^\n]%d", &idum); // expPhase
	fscanf(Parameters, "%*[^\n]%f", &dum); // stop_strain
	fscanf(Parameters, "%*[^\n]%d", &idum); // elastic_write
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_INS
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_elastic
	fscanf(Parameters, "%*[^\n]%f", &start_MSD); // start_MSD
	fclose(Parameters); 
	
	// Read in configuration at the last frame //
	float *rcmx, *rcmy, *rcmz;
	float *rcmx0, *rcmy0, *rcmz0;
	float *rcmxp, *rcmyp, *rcmzp;
	float *MSDyy, *MSDzz, diffy, diffz, sumy, sumz; 
	float *crossyy, *crosszz;
	int step, nConfig, start_MSD_step,m;
	
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
	MSDyy = (float*)malloc(nConfig*sizeof(float));
	MSDzz = (float*)malloc(nConfig*sizeof(float));
	crossyy = (float*)malloc(nfib*sizeof(float));
	crosszz = (float*)malloc(nfib*sizeof(float));
	
	start_MSD_step = int(start_MSD / dt / (float(config_write))); 

	for (m = 0; m < nfib; m++){
		crossyy[m] = 0.0;
		crosszz[m] = 0.0;
	}
	
	for (step = 0; step < nConfig; step++){

		fread(&dum, sizeof(float), 1, center_mass); 
		fread(rcmx, sizeof(float), nfib, center_mass); 
		fread(rcmy, sizeof(float), nfib, center_mass);
		fread(rcmz, sizeof(float), nfib, center_mass);
		if (step == start_MSD_step){
			MSDyy[step] = 0.0; 
			MSDzz[step] = 0.0; 
			memcpy(rcmx0,rcmx,nfib*sizeof(float));
			memcpy(rcmy0,rcmy,nfib*sizeof(float));
			memcpy(rcmz0,rcmz,nfib*sizeof(float));
			memcpy(rcmxp,rcmx,nfib*sizeof(float));
			memcpy(rcmyp,rcmy,nfib*sizeof(float));
			memcpy(rcmzp,rcmz,nfib*sizeof(float));
			for (m = 0; m < nfib; m++){
				crossyy[m] = 0.0;
				crosszz[m] = 0.0;
			}
			continue;
		}
		for (m = 0; m < nfib; m++){
			diffy = (rcmy[m] - rcmyp[m]);
			diffz = (rcmz[m] - rcmzp[m]);
			crossyy[m] -= roundf(diffy/sidey);		
			crosszz[m] -= roundf(diffz/sidez);	
		}
		memcpy(rcmxp,rcmx,nfib*sizeof(float));
		memcpy(rcmyp,rcmy,nfib*sizeof(float));
		memcpy(rcmzp,rcmz,nfib*sizeof(float));
		sumy = 0.0; sumz = 0.0; 
		for (m = 0; m < nfib; m++){
			rcmy[m] += sidey*crossyy[m];			
			rcmz[m] += sidez*crosszz[m];				
			sumy += (rcmy[m] - rcmy0[m])*(rcmy[m]-rcmy0[m]);
			sumz += (rcmz[m] - rcmz0[m])*(rcmz[m]-rcmz0[m]);
		}
		MSDyy[step] = sumy / float(2*nfib);
		MSDzz[step] = sumz / float(2*nfib);
	}
	// printf results to file
	for (step = 0; step < nConfig; step++){
		
		if (step >= start_MSD_step){
			printf("%5.2f %10.4f %10.4f\n", float(step*config_write)*dt, MSDyy[step], MSDzz[step]);
			fprintf(MSD,"%5.2f %10.4f %10.4f\n", float(step*config_write)*dt,MSDyy[step],MSDzz[step]);
		}
	}

	fclose(center_mass); fclose(MSD); 
	free(MSDyy); free(MSDzz); 
	free(rcmx); free(rcmy); free(rcmz); 
	free(rcmx0); free(rcmy0); free(rcmz0); 
	free(rcmxp); free(rcmyp); free(rcmzp); 
	free(crossyy); free(crosszz);
	return 0;
}
