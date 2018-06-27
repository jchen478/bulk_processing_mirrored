#include "cuda_runtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "readData.h"
#include "bnei_set.h"
#include "cell.h"
#include "link.h"
#include "contact.h"
#include "ZeroVar.h"
#include "lead.h"
#include "group.h"
#include "errorCheck.h"

using namespace std;

/**
* \file main.cu
* \brief Characterizes aggregates that remain after redispersion
*
* Contain for loop through trajectories and outputs the following
* - Total number of contact (broken + unbroken)
* - Total number of contact, eliminated of joint duplications
*/

int main()
{

	int maxCon = 50;
	int maxGr = 672; 
	int maxBin = 672;
	int npcn = 5000;

	// 1. Read parameters associated with redispersion cycle ...
	//    and other fiber parameters
	int nfib, nseg, config_write, contact_write; 
	int bdimx, bdimy, bdimz;
	float rp, contact_cutoff, rep_cutoff, overlap;
	float dt, strain, sidex, sidey, sidez;
	float fstar, fact, Astar, decatt;
	float dx_ref, dy_ref, dz_ref;

	FILE *aggChar_input;
	aggChar_input= fopen("aggChar_input.txt", "r");

	fscanf(aggChar_input, "%d", &nfib);
	fscanf(aggChar_input, "%*[^\n]%d", &nseg);
	fscanf(aggChar_input, "%*[^\n]%f", &rp);
	fscanf(aggChar_input, "%*[^\n]%f", &contact_cutoff);
	fscanf(aggChar_input, "%*[^\n]%f", &rep_cutoff);
	fscanf(aggChar_input, "%*[^\n]%f", &overlap);
	fscanf(aggChar_input, "%*[^\n]%f", &dt);
	fscanf(aggChar_input, "%*[^\n]%f", &strain);
	fscanf(aggChar_input, "%*[^\n]%f", &sidex);
	fscanf(aggChar_input, " %f", &sidey);
	fscanf(aggChar_input, " %f", &sidez);
	fscanf(aggChar_input, "%*[^\n]%d", &config_write);
	fscanf(aggChar_input, "%*[^\n]%d", &contact_write);
	fscanf(aggChar_input, "%*[^\n]%f", &fstar);
	fscanf(aggChar_input, "%*[^\n]%f", &fact);
	fscanf(aggChar_input, "%*[^\n]%f", &Astar);
	fscanf(aggChar_input, "%*[^\n]%f", &decatt);
	fscanf(aggChar_input, "%*[^\n]%f", &dx_ref);
	fscanf(aggChar_input, "%*[^\n]%f", &dy_ref);
	fscanf(aggChar_input, "%*[^\n]%f", &dz_ref);
	fscanf(aggChar_input, "%*[^\n]%d", &bdimx);
	fscanf(aggChar_input, "%*[^\n]%d", &bdimy);
	fscanf(aggChar_input, "%*[^\n]%d", &bdimz);
	fclose(aggChar_input); 

	// Constant calculations
	int nConfig, nxbinMax, nybinMax, nzbinMax;

	// number of configurations
	nConfig = int(strain / (dt*float(config_write))) + 1;
	
	// maximum bin dimensions
	nxbinMax = int(floorf(sidex / dx_ref));
	nybinMax = int(floorf(sidey / dy_ref));
	nzbinMax = int(floorf(sidez / dz_ref));
	if (nxbinMax % 2 != 0){
		nxbinMax--;
	}
	if (nybinMax % 2 != 0){
		nybinMax--;
	}
	if (nzbinMax % 2 != 0){
		nzbinMax--;
	}

	// cutoffs
	contact_cutoff = powf((contact_cutoff + 2.0), 2.0);
	rep_cutoff = powf((rep_cutoff + 2.0), 2.0);


	// 2. Open trajectory files
	FILE *rxfile, *ryfile, *rzfile; 
	FILE *pxfile, *pyfile, *pzfile;
	rxfile = fopen("rx.txt", "rb");
	ryfile = fopen("ry.txt", "rb");
	rzfile = fopen("rz.txt", "rb");
	pxfile = fopen("px.txt", "rb");
	pyfile = fopen("py.txt", "rb");
	pzfile = fopen("pz.txt", "rb");

	// read simulation case
	FILE *INSInput;
	INSInput = fopen("INSinput.gen","r");
	int simcase;
	fscanf(INSInput, "%d", &simcase);
	fclose(INSInput);	

	// Based on simulation cases, determine simulation
	// box size at every frame
	float *Lx, *Ly, *Lz; 
	Lx = (float*)malloc(nConfig*sizeof(float)); 
	Ly = (float*)malloc(nConfig*sizeof(float));
	Lz = (float*)malloc(nConfig*sizeof(float));
	 
	// case 0 - basis - sidex does not change
	if (simcase == 0){
		for (int f = 0; f < nConfig; f++){
			Lx[f] = sidex; 
			Ly[f] = sidey; 
			Lz[f] = sidez; 
		}
	}
	// case 1 - redispersion - read sidex from Lbox.txt
	else {
		// Read box info
		FILE *BoxFile;
		float dum;
		int nLbox, box_write;
		BoxFile = fopen("box.gen","r");
		fscanf(BoxFile, "%f", &dum);
		fscanf(BoxFile, "%*[^\n]%f", &dum); 
		fscanf(BoxFile, "%*[^\n]%f", &dum); 
		fscanf(BoxFile, "%*[^\n]%f", &dum); 
		fscanf(BoxFile, "%*[^\n]%d", &box_write); 
		fclose(BoxFile); 
		nLbox = strain / (dt *float(box_write)) + 1;

		// Read box dimensions
		FILE *LboxFile;
		float LxTmp, LyTmp, LzTmp;
		LboxFile = fopen("Lbox.txt","r");
		for (int box = 0; box < nLbox; box++){
			fscanf(LboxFile, "%f %f %f %f %f %f %f",
				&dum, &LxTmp, &LyTmp, &LzTmp, &dum, &dum, &dum);
			if ((box*box_write) % config_write == 0){
				Lx[box*box_write/config_write] = LxTmp;
				Ly[box*box_write/config_write] = LyTmp;
				Lz[box*box_write/config_write] = LzTmp;
			}
		}
		fclose(LboxFile); 
	}


	// Open output file
	FILE *ContactFile;
	ContactFile = fopen("ContactStat.txt", "w");

	// 4. Allocate memory for both host and device operations
	int nxbin, nybin, nzbin, num_groups;
	int total_overlap, total_contact, total_contact_no_joints;
	float total_forc, total_dist_inCon;
	float delta_rx, dx, dy, dz;
	float *rx, *ry, *rz, *px, *py, *pz;
	rx = (float *)malloc(nfib*nseg*sizeof(float));
	ry = (float *)malloc(nfib*nseg*sizeof(float));
	rz = (float *)malloc(nfib*nseg*sizeof(float));
	px = (float *)malloc(nfib*nseg*sizeof(float));
	py = (float *)malloc(nfib*nseg*sizeof(float));
	pz = (float *)malloc(nfib*nseg*sizeof(float));

	float *d_rx, *d_ry, *d_rz; 
	float *d_px, *d_py, *d_pz;
	cudaMalloc((void**)&d_rx, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_ry, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_rz, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_px, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_py, nfib*nseg*sizeof(float));
	cudaMalloc((void**)&d_pz, nfib*nseg*sizeof(float));

	int *d_nfib, *d_nseg, *d_maxBin;
	int *d_npcn, *d_maxCon, *d_maxGr;
	int *d_total_overlap, *d_total_contact, *d_total_contact_no_joints;
	float *d_total_forc, *d_total_dist_inCon;
	float *d_dx, *d_dy, *d_dz, *d_delta_rx;
	float *d_sidex, *d_sidey, *d_sidez;
	float *d_fstar, *d_fact, *d_Astar, *d_decatt;
	cudaMalloc((void**)&d_nfib, sizeof(int));
	cudaMalloc((void**)&d_nseg, sizeof(int));
	cudaMalloc((void**)&d_maxBin, sizeof(int));
	cudaMalloc((void**)&d_maxGr, sizeof(int));
	cudaMalloc((void**)&d_npcn, sizeof(int));
	cudaMalloc((void**)&d_maxCon, sizeof(int));
	cudaMalloc((void**)&d_total_overlap, sizeof(int));
	cudaMalloc((void**)&d_total_contact, sizeof(int));
	cudaMalloc((void**)&d_total_contact_no_joints, sizeof(int));
	cudaMalloc((void**)&d_total_forc, sizeof(float));
	cudaMalloc((void**)&d_total_dist_inCon, sizeof(float));
	cudaMalloc((void**)&d_dx, sizeof(float));
	cudaMalloc((void**)&d_dy, sizeof(float));
	cudaMalloc((void**)&d_dz, sizeof(float));
	cudaMalloc((void**)&d_sidex, sizeof(float));
	cudaMalloc((void**)&d_sidey, sizeof(float));
	cudaMalloc((void**)&d_sidez, sizeof(float));
	cudaMalloc((void**)&d_delta_rx, sizeof(float));
	cudaMalloc((void**)&d_fstar, sizeof(float));
	cudaMalloc((void**)&d_fact, sizeof(float));
	cudaMalloc((void**)&d_Astar, sizeof(float));
	cudaMalloc((void**)&d_decatt, sizeof(float));

	int *bnei, *bin, *list, *bnum;
	int *d_nxbin, *d_nybin, *d_nzbin;
	int *d_bdimx, *d_bdimy, *d_bdimz;
	int *potCon, *potConSize, *d_num_groups;
	int *status, *lead_clist, *nc, *ifiber, *ncnt, *clist_pos, *groupId; 
	cudaMalloc((void**)&bnei, bdimx*bdimy*bdimz*nxbinMax*nybinMax*nzbinMax*sizeof(int));
	cudaMalloc((void**)&bin, nxbinMax*nybinMax*nzbinMax*sizeof(int));
	cudaMalloc((void**)&list, maxBin*nxbinMax*nybinMax*nzbinMax*sizeof(int));
	cudaMalloc((void**)&bnum, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&d_nxbin, sizeof(int));
	cudaMalloc((void**)&d_nybin, sizeof(int));
	cudaMalloc((void**)&d_nzbin, sizeof(int));
	cudaMalloc((void**)&d_bdimx, sizeof(int));
	cudaMalloc((void**)&d_bdimy, sizeof(int));
	cudaMalloc((void**)&d_bdimz, sizeof(int));
	cudaMalloc((void**)&potCon, nfib*nseg * npcn * sizeof(int));
	cudaMalloc((void**)&potConSize, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&status, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&lead_clist, nfib*nseg*maxGr*sizeof(int));
	cudaMalloc((void**)&nc, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&clist_pos, nfib*nseg*maxGr*sizeof(int));
	cudaMalloc((void**)&ifiber, nfib*nseg * 2 * maxGr*sizeof(int));
	cudaMalloc((void**)&ncnt, nfib*nseg*sizeof(int));
	cudaMalloc((void**)&d_num_groups, sizeof(int));
	cudaMalloc((void**)&groupId, nfib*nseg*sizeof(int));


	int *ncpf, *clist; 
	cudaMalloc((void**)&ncpf, (nfib*nseg)*sizeof(int));
	cudaMalloc((void**)&clist, nfib*nseg*maxCon*sizeof(int));


	float *d_contact_cutoff, *d_rep_cutoff, *d_over_cut, *d_rp;
	cudaMalloc((void**)&d_contact_cutoff, sizeof(float));
	cudaMalloc((void**)&d_rep_cutoff, sizeof(float));
	cudaMalloc((void**)&d_over_cut, sizeof(float));
	cudaMalloc((void**)&d_rp, sizeof(float));

	// copy memory to device
	cudaMemcpy(d_bdimx, &bdimx, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bdimy, &bdimy, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_bdimz, &bdimz, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_maxBin, &maxBin, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_maxCon, &maxCon, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_maxGr, &maxGr, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_npcn, &npcn, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_nseg, &nseg, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_decatt, &decatt, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fstar, &fstar, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fact, &fact, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Astar, &Astar, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_contact_cutoff, &contact_cutoff, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rep_cutoff, &rep_cutoff, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_over_cut, &overlap, sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_rp, &rp, sizeof(float), cudaMemcpyHostToDevice);


	for (int step = 0; step < nConfig; step++){

		readData(rxfile, rx, nfib*nseg);
		readData(ryfile, ry, nfib*nseg);
		readData(rzfile, rz, nfib*nseg);
		readData(pxfile, px, nfib*nseg);
		readData(pyfile, py, nfib*nseg);
		readData(pzfile, pz, nfib*nseg);

		delta_rx = float(step*config_write) * dt;
		delta_rx -= lroundf(delta_rx / sidex)*sidex;

		sidex = Lx[step]; 
		sidey = Ly[step];
		sidez = Lz[step];

		nxbin = int(floorf(sidex / dx_ref));
		nybin = int(floorf(sidey / dy_ref));
		nzbin = int(floorf(sidez / dz_ref));
		if (nxbin % 2 != 0){
			nxbin--;
		}
		if (nybin % 2 != 0){
			nybin--;
		}
		if (nzbin % 2 != 0){
			nzbin--;
		}
		dx = sidex / float(nxbin);
		dy = sidey / float(nybin);
		dz = sidez / float(nzbin);
		
		total_contact = 0; 
		total_contact_no_joints = 0;
		total_overlap = 0; 
		total_forc = 0.0;
		total_dist_inCon = 0.0; 
		num_groups = 0; 

		// copy memory to device
		cudaMemcpy(d_rx, rx, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_ry, ry, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_rz, rz, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_px, px, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_py, py, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_pz, pz, nfib*nseg*sizeof(float), cudaMemcpyHostToDevice);

		cudaMemcpy(d_nxbin, &nxbin, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_nybin, &nybin, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_nzbin, &nzbin, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_num_groups, &num_groups, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_total_forc, &total_forc, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_total_dist_inCon, &total_dist_inCon, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_total_contact, &total_contact, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_total_contact_no_joints, &total_contact_no_joints, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_total_overlap, &total_overlap, sizeof(int), cudaMemcpyHostToDevice);
		cudaMemcpy(d_dx, &dx, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_dy, &dy, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_dz, &dz, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_sidex, &sidex, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_sidey, &sidey, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_sidez, &sidez, sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_delta_rx, &delta_rx, sizeof(float), cudaMemcpyHostToDevice);

		// begin characterization
		// 0. zero the variables
		ZeroVar << <nzbin*nybin, nxbin >> >(bin);
		ZeroVar << <nfib, nseg >> >(potConSize);
		ZeroVar << <nfib, nseg >> >(ncpf);
		ZeroVar << <nfib, nseg >> >(status);
		ZeroVar << <nfib, nseg >> >(nc);
		ZeroVar << <nfib, nseg >> >(ncnt);

		// 1. set the neighbors of the bins
		bnei_set << < nzbin*nybin, nxbin >> > (bnei, d_nxbin, d_nybin, d_nzbin,
			d_bdimx, d_bdimy, d_bdimz);

		// 2. Put fibers in the bins
		cell << < nfib/32, nseg*32 >> >(bin, list, bnum, d_rx, d_ry, d_rz, d_nxbin, d_nybin, d_nzbin,
			d_maxBin, d_sidex, d_sidey, d_sidez, d_delta_rx, d_dx, d_dy, d_dz);

		// 3. Find possible contacting pairs and eliminate adjacent segments
		link << < nfib*nseg, bdimx*bdimy*bdimz >> > 
			(bin, list, bnei, bnum, potCon, potConSize, d_npcn, d_nseg, d_nxbin, d_nybin, d_nzbin);

		// 4. check for contacts and obtain contact statistics
		contact << < nfib / 32, nseg * 32 >> >(d_total_contact, d_total_overlap, 
			potCon, potConSize, d_npcn, d_rx, d_ry, d_rz,
			d_px, d_py, d_pz, d_sidex, d_sidey, d_sidez,
			d_delta_rx, d_rp, d_over_cut, d_contact_cutoff, 
			d_rep_cutoff, d_maxCon, ncpf, clist, d_Astar, d_decatt,
			d_fact, d_fstar, d_total_forc, d_total_dist_inCon);

		// 5. Find leaders of each contacting group
		lead << < nfib / 32, nseg * 32 >> >(ncpf, clist, status, nc, lead_clist, d_maxCon, d_maxGr);

		// 6. Eliminate contacts at joints
		group << < nfib / 32, nseg * 32 >>> (ifiber, ncnt, ncpf, d_nseg, 
			clist, status, lead_clist, nc, d_maxCon, clist_pos, d_maxGr,
			groupId, d_num_groups, d_total_contact_no_joints);

		// 5. copy results back to host
		cudaMemcpy(&total_contact, d_total_contact, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&total_contact_no_joints, d_total_contact_no_joints, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&total_overlap, d_total_overlap, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&num_groups, d_num_groups, sizeof(int), cudaMemcpyDeviceToHost);
		cudaMemcpy(&total_forc, d_total_forc, sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(&total_dist_inCon, d_total_dist_inCon, sizeof(float), cudaMemcpyDeviceToHost);

		// 6. output to file or console
		fprintf(ContactFile, "%10.4f %8d %4d %6.3f %8d %6.3f %4d %10.6f %10.6f\n",
			float(step*config_write)*dt, num_groups, total_contact, float(total_contact) / float(nfib),
			total_contact_no_joints, float(total_contact_no_joints) / float(nfib), total_overlap, 
			total_forc/float(total_contact),total_dist_inCon/float(total_contact));
	}

	fclose(ContactFile);

	fclose(rxfile); fclose(ryfile); fclose(rzfile);
	fclose(pxfile); fclose(pyfile); fclose(pzfile);

	free(rx); free(ry); free(rz); 
	free(px); free(py); free(pz);

	free(Lx); free(Ly); free(Lz); 
	
	cudaFree(d_rx); cudaFree(d_ry); cudaFree(d_rz);
	cudaFree(d_px); cudaFree(d_py); cudaFree(d_pz);

	cudaFree(d_nfib); cudaFree(d_nseg); cudaFree(d_maxBin); cudaFree(d_npcn);
	cudaFree(d_dx); cudaFree(d_dy); cudaFree(d_dz); cudaFree(d_delta_rx);
	cudaFree(d_sidex); cudaFree(d_sidey); cudaFree(d_sidez); cudaFree(d_rp);
	cudaFree(d_rep_cutoff); cudaFree(d_contact_cutoff); cudaFree(d_over_cut);
	cudaFree(d_decatt); cudaFree(d_fstar); cudaFree(d_Astar); cudaFree(d_fact);
	cudaFree(ncpf); cudaFree(clist);  cudaFree(d_maxCon);
	cudaFree(d_maxGr); 
	cudaFree(bnei); cudaFree(list); cudaFree(bin); cudaFree(bnum);
	cudaFree(d_nxbin); cudaFree(d_nybin); cudaFree(d_nzbin);
	cudaFree(d_bdimx); cudaFree(d_bdimy); cudaFree(d_bdimz);
	cudaFree(potCon);    cudaFree(potConSize);
	cudaFree(d_total_contact); cudaFree(d_total_overlap); 
	cudaFree(d_total_contact_no_joints);
	cudaFree(d_total_forc);  
	cudaFree(d_total_dist_inCon);

	cudaFree(status); cudaFree(lead_clist);  cudaFree(nc); cudaFree(d_num_groups); 
	cudaFree(ifiber); cudaFree(ncnt); cudaFree(clist_pos); cudaFree(groupId);

    return 0;
}

