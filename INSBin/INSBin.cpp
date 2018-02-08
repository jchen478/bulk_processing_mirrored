#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;


int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *rxfile, *ryfile, *rzfile;
	FILE *ISV, *ISV_results, *Rfile; 
	
	// open files //
	Parameters = fopen("Parameters.in", "r");

	rxfile = fopen("../../rx.txt", "rb");
	
	if (rxfile == 0){
		printf("fail to open rx file");
		return 1; 
	}

	ryfile = fopen("../../ry.txt", "rb");
	rzfile = fopen("../../rz.txt", "rb");

	ISV = fopen("ISV.txt", "w"); 
	ISV_results = fopen("ISV_results.txt", "w");
	Rfile = fopen("R.txt", "w");

	int nfib, nseg, config_write, nxbin, nybin, nzbin, n_exp, n_conc, idum, INSPhase; 
	float dt, strain, rp, sidex, sidey, sidez;
	float start_INS, start_exp, eq_exp, concFac, expFac, dum; 

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib); //nfib
	fscanf(Parameters, "%*[^\n]%d", &nseg); // nseg
	fscanf(Parameters, "%*[^\n]%f", &rp); // rp
	fscanf(Parameters, "%*[^\n]%f", &dum); // kb
	fscanf(Parameters, "%*[^\n]%f", &dum); // mu_stat
	fscanf(Parameters, "%*[^\n]%f", &dum);	// mu_kin
	fscanf(Parameters, "%*[^\n]%f", &dum);	// contact_cutoff
	fscanf(Parameters, "%*[^\n]%f", &dum);	// rep_cutoff
	fscanf(Parameters, "%*[^\n]%f", &dum); // overlap
	fscanf(Parameters, "%*[^\n]%f", &dt);	// dt
	fscanf(Parameters, "%*[^\n]%f", &strain); // strain
	fscanf(Parameters, "%*[^\n]%f", &sidex); // sidex
	fscanf(Parameters, " %f", &sidey); // sidey
	fscanf(Parameters, " %f", &sidez); // sidez
	fscanf(Parameters, "%*[^\n]%f", &dum); // fraction_rp
	fscanf(Parameters, "%*[^\n]%d", &config_write); // config_write
	fscanf(Parameters, "%*[^\n]%d", &idum); // contact_write
	fscanf(Parameters, "%*[^\n]%f", &dum); // fstar
	fscanf(Parameters, "%*[^\n]%f", &dum); // fact
	fscanf(Parameters, "%*[^\n]%f", &dum); // Astar
	fscanf(Parameters, "%*[^\n]%f", &dum); // decatt
	fscanf(Parameters, "%*[^\n]%f", &dum); // delta_rx
	fscanf(Parameters, "%*[^\n]%f", &dum); // duxdx
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum); // fac
	fscanf(Parameters, "%*[^\n]%f", &dum); // elf
	fscanf(Parameters, "%*[^\n]%f", &dum); // dx
	fscanf(Parameters, "%*[^\n]%f", &dum); // dy
	fscanf(Parameters, "%*[^\n]%f", &dum); // dz
	fscanf(Parameters, "%*[^\n]%d", &idum); // bdimx
	fscanf(Parameters, "%*[^\n]%d", &idum); // bdimy
	fscanf(Parameters, "%*[^\n]%d", &idum); // bdimz
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
	fscanf(Parameters, "%*[^\n]%f", &concFac); // conFac
	fscanf(Parameters, "%*[^\n]%f", &expFac); // expFac
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_conc
	fscanf(Parameters, "%*[^\n]%f", &start_exp); // start_exp
	fscanf(Parameters, "%*[^\n]%f", &dum); // eq_conc
	fscanf(Parameters, "%*[^\n]%f", &eq_exp); // eq_exp
	fscanf(Parameters, "%*[^\n]%d", &n_conc); // n_conc
	fscanf(Parameters, "%*[^\n]%d", &n_exp); // n_exp
	fscanf(Parameters, "%*[^\n]%d", &idum); // dx_entropy
	fscanf(Parameters, "%*[^\n]%d", &idum); // dy_entropy
	fscanf(Parameters, "%*[^\n]%d", &idum);	// dz_entropy
	fscanf(Parameters, "%*[^\n]%d", &nxbin); // nxBin_dist
	fscanf(Parameters, "%*[^\n]%d", &nybin); // nyBin_dist
	fscanf(Parameters, "%*[^\n]%d", &nzbin); // nzBin_dist
	fscanf(Parameters, "%*[^\n]%d", &dum); // concPhase
	fscanf(Parameters, "%*[^\n]%d", &INSPhase); // INSPhase
	fscanf(Parameters, "%*[^\n]%d", &dum); // concPhase
	fscanf(Parameters, "%*[^\n]%f", &start_INS); // start_INS

	if (INSPhase == 1){
		sidex *= powf(concFac, n_conc)*powf(expFac, n_exp); 
		sidey *= powf(concFac, n_conc)*powf(expFac, n_exp);
		sidez *= powf(concFac, n_conc)*powf(expFac, n_exp);
	}

	printf("strain %f %f %f %f %d %d %d\n", strain, sidex, sidey, sidez, nxbin, nybin, nzbin); 

	int n_afterExpand; 
	int step, frames, n_start_calc, m, i, ii, jj, kk, xloc, yloc, zloc, ind, offset, pairInd;
	float start_calc, N_afterExpand, volfrac, pi;
	float dx, dy, dz, *bin, Vseg, Vcell, phi; 
	float *rx, *ry, *rz, rxmi, rymi, rzmi, *intensity;
	float *autox, *autoy, *autoz, *Rx,*Ry, *Rz, prev;
	float sumx, sumy, sumz, avgNseg, binNorm; 
	float *Sx, *Sy, *Sz, *Vx, *Vy, *Vz, dRx, dRy, dRz; 
	float intensitysum, intensityavg, intensitysig; 
	float Sxsum, Sysum, Szsum, Vxsum, Vysum, Vzsum; 
	float Sxavg, Syavg, Szavg, Vxavg, Vyavg, Vzavg;
	float Sxsig, Sysig, Szsig, Vxsig, Vysig, Vzsig;
	float intensitycalc; 
	float Sxcalc, Sycalc, Szcalc, Vxcalc, Vycalc, Vzcalc; 
	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));
	bin = (float*)calloc(nxbin*nybin*nzbin, sizeof(float));
	Rx = (float*)calloc(nxbin / 2, sizeof(float));
	Ry = (float*)calloc(nybin / 2, sizeof(float));
	Rz = (float*)calloc(nzbin / 2, sizeof(float));
	autox = (float*)calloc(nxbin / 2, sizeof(float)); 
	autoy = (float*)calloc(nybin / 2, sizeof(float));
	autoz = (float*)calloc(nzbin / 2, sizeof(float));
	dx = sidex / float(nxbin);
	dy = sidey / float(nybin);
	dz = sidez / float(nzbin);
	intensity = (float*)calloc(200000, sizeof(float));
	Sx = (float*)calloc(200000, sizeof(float)); 
	Sy = (float*)calloc(200000, sizeof(float));
	Sz = (float*)calloc(200000, sizeof(float));
	Vx = (float*)calloc(200000, sizeof(float));
	Vy = (float*)calloc(200000, sizeof(float));
	Vz = (float*)calloc(200000, sizeof(float));
	pi = 3.141592654; 
	n_afterExpand = 0; 
	N_afterExpand = 0.0; 
	intensitysum = 0.0; intensitysig = 0.0; 
	Sxsum = 0.0; Sysum = 0.0; Szsum = 0.0; 
	Vxsum = 0.0; Vysum = 0.0; Vzsum = 0.0; 
	Sxsig = 0.0; Sysig = 0.0; Szsig = 0.0;
	Vxsig = 0.0; Vysig = 0.0; Vzsig = 0.0;
	start_calc = start_exp + float(n_exp)*eq_exp; 	
	if(start_INS != 0){
		start_calc = start_INS;
	}
	n_start_calc = start_calc / dt / float(config_write); 
	frames = strain / dt / float(config_write) + 1; 
	volfrac = 2.0*nfib*nseg*rp*pi / (sidex*sidey*sidez); 
	Vcell = dx*dy*dz; 
	Vseg = 2.0 * pi*rp; 
	avgNseg = float(nseg*nfib) / float(nxbin*nybin*nzbin); 

	dRx = 2.0 / float(nxbin); 
	dRy = 2.0 / float(nybin);
	dRz = 2.0 / float(nzbin);

	for (ii = 0; ii < nxbin / 2; ii++){
		Rx[ii] = 2.0*float(ii) / float(nxbin); 
	}
	for (ii = 0; ii < nybin / 2; ii++){
		Ry[ii] = 2.0*float(ii) / float(nybin);
	}
	for (ii = 0; ii < nzbin / 2; ii++){
		Rz[ii] = 2.0*float(ii) / float(nzbin);
	}

	for (step = 0; step < frames; step++){
		// read data
		fread(&dum, sizeof(float), 1, rxfile);
		fread(&dum, sizeof(float), 1, ryfile);
		fread(&dum, sizeof(float), 1, rzfile);
		fread(rx, sizeof(float), nfib*nseg, rxfile);
		fread(ry, sizeof(float), nfib*nseg, ryfile);
		fread(rz, sizeof(float), nfib*nseg, rzfile);
		// zero variables
		for (ii = 0; ii < nxbin / 2; ii++){
			autox[ii] = 0.0;
		}
		for (ii = 0; ii < nybin / 2; ii++){
			autoy[ii] = 0.0;
		}
		for (ii = 0; ii < nzbin / 2; ii++){
			autoz[ii] = 0.0;
		}
		for (ii = 0; ii < nxbin; ii++){
			for (jj = 0; jj < nybin; jj++){
				for (kk = 0; kk < nzbin; kk++){
					ind = ii + jj*nxbin + kk*nxbin*nybin;
					bin[ind] = 0.0;
				}
			}
		}
		binNorm = 0.0;
		intensitycalc = 0.0; 
		Sxcalc = 0.0; Sycalc = 0.0; Szcalc = 0.0;
		Vxcalc = 0.0; Vycalc = 0.0; Vzcalc = 0.0;
		// binning
		for (m = 0; m < nfib; m++){
			for (i = 0; i < nseg; i++){
				rxmi = rx[m*nseg + i];
				rymi = ry[m*nseg + i];
				rzmi = rz[m*nseg + i];					
				rxmi -= sidex*roundf(rxmi / sidex);
				rymi -= sidey*roundf(rymi / sidey);
				rzmi -= sidez*roundf(rzmi / sidez);
				rxmi += sidex / 2.0;
				rymi += sidey / 2.0;
				rzmi += sidez / 2.0;
				xloc = int(rxmi / dx);
				yloc = int(rymi / dy);
				zloc = int(rzmi / dz);
				if (xloc == nxbin) xloc--; 
				if (yloc == nybin) yloc--;
				if (zloc == nzbin) zloc--;
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
				bin[xloc + yloc*nxbin + zloc*nxbin*nybin] += 1.0; 
			}
		}
		// norm calculation
		for (ii = 0; ii < nxbin; ii++){
			for (jj = 0; jj < nybin; jj++){
				for (kk = 0; kk < nzbin; kk++){
					binNorm += (bin[ii + jj*nxbin + kk*nxbin*nybin] - avgNseg)*(bin[ii + jj*nxbin + kk*nxbin*nybin] - avgNseg);
				}
			}
		}
		binNorm /= float(nxbin*nybin*nzbin);
		// autocorrelation (if nxbin = nybin = nzbin)
		for (offset = 0; offset < nxbin / 2; offset++){
			sumx = 0.0; sumy = 0.0; sumz = 0.0;
			for (ii = 0; ii < nxbin; ii++){
				pairInd = (ii + offset) - nxbin * floorf((ii + offset) / nxbin);
				for (jj = 0; jj < nybin; jj++){
					for (kk = 0; kk < nzbin; kk++){
						sumx += (bin[ii + jj*nxbin + kk*nxbin*nybin] - avgNseg)*(bin[pairInd + jj*nxbin + kk*nxbin*nybin] - avgNseg);
					}
				}
			}	
			for (jj = 0; jj < nybin; jj++){
				pairInd = (jj + offset) - nybin * floorf((jj + offset) / nybin);
				for (ii = 0; ii < nxbin; ii++){
					for (kk = 0; kk < nzbin; kk++){							
						sumy += (bin[ii + jj*nxbin + kk*nxbin*nybin] - avgNseg)*(bin[ii + pairInd*nxbin + kk*nxbin*nybin] - avgNseg);
					}
				}
			}
			for (kk = 0; kk < nzbin; kk++){
				pairInd = (kk + offset) - nzbin * floorf((kk + offset) / nzbin);
				for (jj = 0; jj < nybin; jj++){
					for (ii = 0; ii < nxbin; ii++){
						sumz += (bin[ii + jj*nxbin + kk*nxbin*nybin] - avgNseg)*(bin[ii + jj*nxbin + pairInd*nxbin*nybin] - avgNseg);
					}
				}
			}
			autox[offset] = sumx / float(nybin*nzbin);
			autoy[offset] = sumy / float(nxbin*nzbin);
			autoz[offset] = sumz / float(nxbin*nybin);
		}
		// normalize and print autocorrelation
		fprintf(Rfile, "%f\n", float(step)*dt*config_write);
		for (ii = 0; ii < nxbin/2; ii++){
			autox[ii] /= binNorm*float(nxbin); 
			autoy[ii] /= binNorm*float(nybin);
			autoz[ii] /= binNorm*float(nzbin);
			fprintf(Rfile, "%10.6f %10.6f %10.6f\n", autox[ii], autoy[ii], autoz[ii]);
		}
		// calculate length and volume scale
		for (ii = 1; ii < nxbin / 2; ii++){
			prev = autox[ii - 1];
			if (autox[ii] >= 0.0){
				Sxcalc += prev + autox[ii];
				Vxcalc += prev*Rx[ii - 1] * Rx[ii - 1] + autox[ii] * Rx[ii] * Rx[ii];
			}
			else
				break;
		}
		for (ii = 1; ii < nybin / 2; ii++){
			prev = autoy[ii - 1];
			if (autoy[ii] >= 0.0){
				Sycalc += prev + autoy[ii];
				Vycalc += prev*Ry[ii - 1] * Ry[ii - 1] + autoy[ii] * Ry[ii] * Ry[ii];
			}
			else
				break;
		}
		for (ii = 1; ii < nzbin / 2; ii++){
			prev = autoz[ii - 1];
			if (autoz[ii] >= 0.0){
				Szcalc += prev + autoz[ii];
				Vzcalc += prev*Rz[ii - 1] * Rz[ii - 1] + autoz[ii] * Rz[ii] * Rz[ii];
			}
			else
				break;
		}
		Sxcalc *= dRx / 2.0;
		Vxcalc *= 2.0*pi*dRx / 2.0;
		Sycalc *= dRy / 2.0;
		Vycalc *= 2.0*pi*dRy / 2.0;
		Szcalc *= dRz / 2.0;
		Vzcalc *= 2.0*pi*dRz / 2.0;
		// calculate intensity
		for (ii = 0; ii < nxbin; ii++){
			for (jj = 0; jj < nybin; jj++){
				for (kk = 0; kk < nzbin; kk++){
					ind = ii + jj*nxbin + kk*nxbin*nybin;
					phi = bin[ind] * Vseg / Vcell;
					intensitycalc += (phi - volfrac)*(phi - volfrac);
				}
			}
		}
		intensitycalc /= float(nxbin*nybin*nzbin)*volfrac*(1 - volfrac);
		// print data
		fprintf(ISV, "%8.4f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",
			float(step)*dt*config_write, intensitycalc,
			Sxcalc, Sycalc, Szcalc, Vxcalc, Vycalc, Vzcalc);
		// save calculations after specified step
		if (step >= n_start_calc){			
			Sxsum += Sxcalc; 
			Sysum += Sycalc;
			Szsum += Szcalc;			
			Vxsum += Vxcalc;
			Vysum += Vycalc;
			Vzsum += Vzcalc;
			intensitysum += intensitycalc;
			Sx[n_afterExpand] = Sxcalc;
			Vx[n_afterExpand] = Vxcalc;
			Sy[n_afterExpand] = Sycalc;
			Vy[n_afterExpand] = Vycalc;
			Sz[n_afterExpand] = Szcalc;
			Vz[n_afterExpand] = Vzcalc;
			intensity[n_afterExpand] = intensitycalc;
			n_afterExpand += 1;
			N_afterExpand += 1.0; 
		}
	}
	// calculate average
	intensityavg = intensitysum / N_afterExpand;
	Sxavg = Sxsum / N_afterExpand; 
	Syavg = Sysum / N_afterExpand;
	Szavg = Szsum / N_afterExpand;
	Vxavg = Vxsum / N_afterExpand;
	Vyavg = Vysum / N_afterExpand;
	Vzavg = Vzsum / N_afterExpand;


	for (i = 0; i < n_afterExpand; i++){
		intensitysig = (intensity[i] - intensityavg)*(intensity[i] - intensityavg);
		Sxsig = (Sx[i] - Sxavg)*(Sx[i] - Sxavg); 
		Sysig = (Sy[i] - Syavg)*(Sy[i] - Syavg);
		Szsig = (Sz[i] - Szavg)*(Sz[i] - Szavg);
		Vxsig = (Vx[i] - Vxavg)*(Vx[i] - Vxavg);
		Vysig = (Vy[i] - Vyavg)*(Vy[i] - Vyavg);
		Vzsig = (Vz[i] - Vzavg)*(Vz[i] - Vzavg);
	}

	intensitysig /= N_afterExpand;
	Sxsig /= N_afterExpand; 
	Sysig /= N_afterExpand;
	Szsig /= N_afterExpand;
	Vxsig /= N_afterExpand;
	Vysig /= N_afterExpand;
	Vzsig /= N_afterExpand;

	intensitysig = sqrtf(intensitysig);
	Sxsig = sqrtf(Sxsig);
	Sysig = sqrtf(Sysig);
	Szsig = sqrtf(Szsig);
	Vxsig = sqrtf(Vxsig);
	Vysig = sqrtf(Vysig);
	Vzsig = sqrtf(Vzsig);

	
	fprintf(ISV_results, "----------------------------\n");
	fprintf(ISV_results, "Time averaged values:\n");
	fprintf(ISV_results, "----------------------------\n");
	fprintf(ISV_results, "Intensity\n");
	fprintf(ISV_results, "%10.6f\n", intensityavg);
	fprintf(ISV_results, "Linear scale S\n");
	fprintf(ISV_results, "%10.6f\n%10.6f\n%10.6f\n", Sxavg, Syavg, Szavg);
	fprintf(ISV_results, "Volume scale V\n");
	fprintf(ISV_results, "%10.6f\n%10.6f\n%10.6f\n", Vxavg, Vyavg, Vzavg);

	fprintf(ISV_results, "----------------------------\n");
	fprintf(ISV_results, "Standard deviation:\n");
	fprintf(ISV_results, "----------------------------\n");
	fprintf(ISV_results, "Intensity\n");
	fprintf(ISV_results, "%10.6f\n", intensitysig);
	fprintf(ISV_results, "Linear scale S\n");
	fprintf(ISV_results, "%10.6f\n%10.6f\n%10.6f\n", Sxsig, Sysig, Szsig);
	fprintf(ISV_results, "Volume scale V\n");
	fprintf(ISV_results, "%10.6f\n%10.6f\n%10.6f\n", Vxsig, Vysig, Vzsig);

	fprintf(ISV_results, "----------------------------\n");
	fprintf(ISV_results, "Standard error of the mean:\n");
	fprintf(ISV_results, "----------------------------\n");
	fprintf(ISV_results, "Intensity\n");
	fprintf(ISV_results, "%10.6f\n", intensitysig/sqrtf(N_afterExpand));
	fprintf(ISV_results, "Linear scale S\n");
	fprintf(ISV_results, "%10.6f\n%10.6f\n%10.6f\n", Sxsig/sqrtf(N_afterExpand) ,Sysig/sqrtf(N_afterExpand) , Szsig/sqrtf(N_afterExpand) );
	fprintf(ISV_results, "Volume scale V\n");
	fprintf(ISV_results, "%10.6f\n%10.6f\n%10.6f\n", Vxsig/sqrtf(N_afterExpand) , Vysig/sqrtf(N_afterExpand) , Vzsig/sqrtf(N_afterExpand) );
	
	free(rx); free(ry); free(rz); free(bin); 
	free(autox); free(autoy); free(autoz); 
	free(Rx); free(Ry); free(Rz); 
	free(intensity); free(Sx); free(Sy); free(Sz); 
	free(Vx); free(Vy); free(Vz); 

	fclose(Parameters); fclose(rxfile); fclose(ryfile); fclose(rzfile); 
	fclose(ISV); fclose(ISV_results); fclose(Rfile); 

	return 0; 

}
