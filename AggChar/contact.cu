#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime_api.h>
#include <math.h>
#include "contact.h"

using namespace std;

__device__ void parallel_sort_para(int mi, int nj, float sx, float sy, float sz, float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj, float pdotp, float rp,
	float *xmin, float *ymin);

__global__ void contact(int *total_contact, int *total_overlap,
	int *potCon, int *potConSize, int *npcnpt,
	float *rx, float *ry, float *rz,
	float *px, float *py, float *pz,
	float *sidexpt, float *sideypt, float *sidezpt,
	float *delta_rxpt, float *rppt, float *over_cutpt,
	float *contact_cutoffpt, float *rep_cutoffpt,
	int *maxConpt, int *ncpf, int *clist){
 

	int mi = threadIdx.x + blockIdx.x*blockDim.x; 
	int nPair = potConSize[mi];
	if (nPair == 0) return;

	int npcn = *npcnpt;
	int maxCon = *maxConpt; 

	float rp = *rppt;
	float over_cut = *over_cutpt;
	float sidex = *sidexpt;
	float sidey = *sideypt;
	float sidez = *sidezpt;
	float delta_rx = *delta_rxpt;
	float contact_cutoff = *contact_cutoffpt;
	float rep_cutoff = *rep_cutoffpt;

	float rxmi, rymi, rzmi, rxnj, rynj, rznj;
	float pxmi, pymi, pzmi, pxnj, pynj, pznj;
	float sxx, syy, szz, corx, cory, corz;
	float rxmi_shift, rymi_shift, rzmi_shift; 
	float pdotp, xmin, ymin, dx, dy, dz, sep; 
	float xi[9], yj[9], gij, nijx, nijy, nijz, forc;
	float Gijx, Gijy, Gijz, Gjix, Gjiy, Gjiz, sep_tmp;

	int nP, nj, ipos, ith, oldmi, oldnj;

	rxmi = rx[mi]; rymi = ry[mi]; rzmi = rz[mi];
	pxmi = px[mi]; pymi = py[mi]; pzmi = pz[mi];
	
	for (nP = 0; nP < nPair; nP++){

		nj = potCon[mi * npcn + nP]; 
		rxnj = rx[nj]; rynj = ry[nj]; rznj = rz[nj];
		pxnj = px[nj]; pynj = py[nj]; pznj = pz[nj];

		// find minimum image (for shear flow system)
		sxx = rxnj - rxmi;
		syy = rynj - rymi;
		szz = rznj - rzmi;	
		cory = roundf(syy / sidey);
		corz = roundf(szz / sidez);
		sxx = sxx - corz*delta_rx;
		corx = roundf(sxx / sidex);
		sxx = sxx - corx*sidex;
		syy = syy - cory*sidey;
		szz = szz - corz*sidez;
		rxmi_shift = rxnj - sxx;
		rymi_shift = rynj - syy;
		rzmi_shift = rznj - szz;
		pdotp = pxmi*pxnj + pymi*pynj + pzmi*pznj; 
		xmin = (-(pxnj * sxx + pynj * syy + pznj * szz)* pdotp
			+ (pxmi * sxx + pymi * syy + pzmi * szz))
			/ (1.0 - pdotp*pdotp);
		ymin = ((pxmi * sxx + pymi * syy + pzmi * szz)* pdotp
			- (pxnj * sxx + pynj * syy + pznj * szz))
			/ (1.0 - pdotp*pdotp);

		dx = rxnj + ymin*pxnj - rxmi_shift - xmin*pxmi;
		dy = rynj + ymin*pynj - rymi_shift - xmin*pymi;
		dz = rznj + ymin*pznj - rzmi_shift - xmin*pzmi;
		sep = dx*dx + dy*dy + dz*dz;

		ipos = 8;
		yj[0] = rp;
		xi[0] = pxmi*sxx + pymi*syy + pzmi*szz + yj[0] * pdotp; 
		yj[1] = -rp;
		xi[1] = pxmi*sxx + pymi*syy + pzmi*szz + yj[1] * pdotp;
		xi[2] = rp;
		yj[2] = -(pxnj*sxx + pynj*syy + pznj*szz) + xi[2] * pdotp;
		xi[3] = -rp;
		yj[3] = -(pxnj*sxx + pynj*syy + pznj*szz) + xi[3] * pdotp;
		xi[4] = rp;    yj[4] = rp;
		xi[5] = rp;    yj[5] = -rp;
		xi[6] = -rp;   yj[6] = rp;
		xi[7] = -rp;   yj[7] = -rp;
		xi[8] = xmin;  yj[8] = ymin;
		
		// Check if segments are parallel
		if (fabsf(pdotp*pdotp-1.0) <= 1.0e-6) {
			parallel_sort_para(mi, nj, sxx, syy, szz, pxmi, pymi, pzmi,
				pxnj, pynj, pznj, pdotp, rp, &xmin, &ymin);
			sep = (sxx + ymin*pxnj - xmin*pxmi)*(sxx + ymin*pxnj - xmin*pxmi) +
				(syy + ymin*pynj - xmin*pymi)*(syy + ymin*pynj - xmin*pymi) +
				(szz + ymin*pznj - xmin*pzmi)*(szz + ymin*pznj - xmin*pzmi);
		}
		else if (sep < rep_cutoff && (fabsf(xmin) >= rp || fabsf(ymin) >= rp)){
			sep = 1000.0;
			// check which end-side or end-end separation
			// is the smallest
			for (ith = 0; ith < 8; ith++){
				sep_tmp = (sxx + yj[ith] * pxnj - xi[ith] * pxmi)*(sxx + yj[ith] * pxnj - xi[ith] * pxmi) +
					(syy + yj[ith] * pynj - xi[ith] * pymi)*(syy + yj[ith] * pynj - xi[ith] * pymi) +
					(szz + yj[ith] * pznj - xi[ith] * pzmi)*(szz + yj[ith] * pznj - xi[ith] * pzmi);
				if (sep_tmp < sep && fabsf(xi[ith]) <= rp && fabsf(yj[ith]) <= rp){
					sep = sep_tmp;
					ipos = ith;
				}
			}
			xmin = xi[ipos];
			ymin = yj[ipos];
		}		
		gij = sqrtf(sep);	
		if (gij < 2.0){
			atomicAdd(total_overlap, 1); // overs
		}
		if (gij < over_cut){
			gij = over_cut;
		}
		if (sep < contact_cutoff){

			// contact, mi, nj
			oldmi = atomicAdd(ncpf + mi, 1);
			oldnj = atomicAdd(ncpf + nj, 1);
			clist[mi*maxCon + oldmi] = nj;
			clist[nj*maxCon + oldnj] = mi;

			if (mi < nj)
				atomicAdd(total_contact, 1); // total_contacts (broken + unbroken)
		}	
	}

}
__device__ void parallel_sort_para(int mi, int nj, float sx, float sy, float sz, float pxmi, float pymi, float pzmi,
	float pxnj, float pynj, float pznj, float pdotp, float rp,
	float *xmin, float *ymin){

	//printf("accessing parallel_sort_para\n");
	float posneg, pn2, dist, sijp, sijm, sjip, sjim;

	//printf("%4d %4d %15.10f p %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n", mi, nj, pdotp, pxmi, pymi, pzmi, pxnj, pynj, pznj); 

	// The different end point to fiber contact points
	sijp = pxmi*sx + pymi*sy + pzmi*sz + rp*pdotp;
	sijm = pxmi*sx + pymi*sy + pzmi*sz - rp*pdotp;
	sjip = -(pxnj*sx + pynj*sy + pznj*sz) + rp*pdotp;
	sjim = -(pxnj*sx + pynj*sy + pznj*sz) - rp*pdotp;

	//printf("parallel\n"); 
	//printf("%4d %4d sx sy sz %15.10f %15.10f %15.10f sijp sijm sjip sjim %15.10f %15.10f %15.10f %15.10f\n", mi, nj, sx, sy, sz, sijp, sijm, sjip, sjim); 

	// for fiber i
	if (fabsf(sijp) < fabsf(sijm)){
		*xmin = sijp;
		posneg = 1.0;
	}
	else if (fabsf(sijp) > fabsf(sijm)){
		*xmin = sijm;
		posneg = -1.0;
	}
	else{
		*xmin = 0.0;
		posneg = 0.0;
	}
	if (*xmin >= rp){
		*xmin = rp;
	}
	if (*xmin <= -rp){
		*xmin = -rp;
	}
	// for fiber j
	if (fabsf(sjip) < fabsf(sjim)){
		*ymin = sjip;
	}
	else if (fabsf(sjip) > fabsf(sjim)){
		*ymin = sjim;
	}
	else{
		*ymin = 0.0;
		posneg = 0.0;
	}
	if (*ymin >= rp){
		*ymin = rp;
	}
	if (*ymin <= -rp){
		*ymin = -rp;
	}
	//printf("xmin ymin in %12.8f %12.8f\n", *xmin, *ymin);
	//printf("xmin ymin out %4d %4d %16.10f %16.10f\n", mi, nj, *xmin, *ymin); 
	if (fabsf(*xmin) < rp && fabsf(*ymin) < rp){
		if (pdotp > 0.0){
			pn2 = 1.0;
		}
		else{
			pn2 = -1.0;
		}
		dist = (rp + posneg**xmin) / 2.0;
		*xmin = *xmin - posneg*dist;
		*ymin = *ymin + posneg*pn2*dist;
		//printf("xmin ymin in %12.8f %12.8f\n", *xmin, *ymin);
	}
}
