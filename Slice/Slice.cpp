#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;


int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *center_mass;
	FILE *Binning_result_xy, *Binning_result_xz, *Binning_result_yz;
	FILE *Binning_parameters, *Binning_3d;
	FILE *rxfile, *ryfile, *rzfile;

	Parameters = fopen("Parameters.in", "r");
	center_mass = fopen("center_mass.txt", "rb");
	rxfile = fopen("rx.txt", "rb");
	ryfile = fopen("ry.txt", "rb");
	rzfile = fopen("rz.txt", "rb");

	Binning_3d = fopen("Binning_3d.txt", "w");
	Binning_result_xy = fopen("Binning_result_xy.txt", "w");
	Binning_result_xz = fopen("Binning_result_xz.txt", "w");
	Binning_result_yz = fopen("Binning_result_yz.txt", "w");
	Binning_parameters = fopen("Binning_parameters.txt", "w");

	if (Parameters == NULL || center_mass == NULL) {
		perror("Error");
	}

	////////////////////////////////////////
	//               variables            //
	////////////////////////////////////////

	int nfib, nseg, config_write;
	int nxbin, nybin, nzbin;
	float dt, strain, sidex, sidey, sidez, dx, dy, dz;

	int idum;
	float dum;

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &dum); // rp
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
	fscanf(Parameters, "%*[^\n]%f", &dum); // fstar
	fscanf(Parameters, "%*[^\n]%f", &dum); // fact
	fscanf(Parameters, "%*[^\n]%f", &dum); // Astar
	fscanf(Parameters, "%*[^\n]%f", &dum); // decatt
	fscanf(Parameters, "%*[^\n]%f", &dum); // delta_rx
	fscanf(Parameters, "%*[^\n]%f", &dum); // duxdx
	fscanf(Parameters, "%*[^\n]%f", &dum); // duydx
	fscanf(Parameters, "%*[^\n]%f", &dum); // duzdx
	fscanf(Parameters, "%*[^\n]%f", &dum); // duxdy
	fscanf(Parameters, "%*[^\n]%f", &dum); // duydy
	fscanf(Parameters, "%*[^\n]%f", &dum); // duzdy
	fscanf(Parameters, "%*[^\n]%f", &dum); // duxdz
	fscanf(Parameters, "%*[^\n]%f", &dum); // duydz
	fscanf(Parameters, "%*[^\n]%f", &dum); // duzdz
	fscanf(Parameters, "%*[^\n]%d", &dum); // fac
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
	fscanf(Parameters, "%*[^\n]%f", &dum); // concFac
	fscanf(Parameters, "%*[^\n]%f", &dum); // expFac
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_conc
	fscanf(Parameters, "%*[^\n]%f", &dum); // start_exp
	fscanf(Parameters, "%*[^\n]%f", &dum); // eq_conc
	fscanf(Parameters, "%*[^\n]%f", &dum); // eq_exp
	fscanf(Parameters, "%*[^\n]%d", &idum); // n_conc 
	fscanf(Parameters, "%*[^\n]%d", &idum); // n_exp 
	fscanf(Parameters, "%*[^\n]%f", &dum); // dx_entropy 
	fscanf(Parameters, "%*[^\n]%f", &dum); // dy_entropy
	fscanf(Parameters, "%*[^\n]%f", &dum); // dz_entropy
	fscanf(Parameters, "%*[^\n]%d", &nxbin); // nxbin 
	fscanf(Parameters, "%*[^\n]%d", &nybin); // nybin 
	fscanf(Parameters, "%*[^\n]%d", &nzbin); // nzbin 


	// vairables

	float *rcmx, *rcmy, *rcmz;
	float *rx, *ry, *rz;
	int frames, step, m, i, j, k;
	float rxmi, rymi, rzmi;
	float rcmxm, rcmym, rcmzm;
	int xloc, yloc, zloc;
	int *bin,*binxy, *binxz, *binyz;

	rcmx = (float*)malloc(nfib*sizeof(float));
	rcmy = (float*)malloc(nfib*sizeof(float));
	rcmz = (float*)malloc(nfib*sizeof(float));
	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));

	frames = int(strain / dt) / config_write + 1;
	
	bin = (int*)calloc(nxbin*nybin*nzbin, sizeof(int)); 
	binxy = (int*)calloc(nxbin*nybin, sizeof(int));
	binxz = (int*)calloc(nxbin*nzbin, sizeof(int));
	binyz = (int*)calloc(nybin*nzbin, sizeof(int));

	dx = sidex / float(nxbin);
	dy = sidey / float(nybin);
	dz = sidez / float(nzbin);

	fprintf(Binning_parameters, "%4d %4d %4d\n", nxbin, nybin, nzbin);
	fprintf(Binning_parameters, "%15.6f %15.6f %15.6f\n", dx, dy, dz);

	for (step = 0; step < frames; step++){

		fread(&dum, sizeof(float), 1, center_mass);
		fread(&dum, sizeof(float), 1, rxfile);
		fread(&dum, sizeof(float), 1, ryfile);
		fread(&dum, sizeof(float), 1, rzfile);
		fread(rcmx, sizeof(float), nfib, center_mass);
		fread(rcmy, sizeof(float), nfib, center_mass);
		fread(rcmz, sizeof(float), nfib, center_mass);
		fread(rx, sizeof(float), nfib*nseg, rxfile);
		fread(ry, sizeof(float), nfib*nseg, ryfile);
		fread(rz, sizeof(float), nfib*nseg, rzfile);

		if (step == frames - 1){
			printf("dum %f\n", dum); 
			for (m = 0; m < nfib; m++){
				/*rcmxm = rcmx[m] + sidex / 2.0;
				rcmym = rcmy[m] + sidey / 2.0;
				rcmzm = rcmz[m] + sidez / 2.0;
				xloc = int(rcmxm / dx);
				yloc = int(rcmym / dy);
				zloc = int(rcmzm / dz);
				if (xloc >= nxbin || yloc >= nybin || zloc >= nzbin)
					printf("%8d %13.6f %13.6f %13.6f\n", m + 1, rcmxm, rcmym, rcmzm);
				if (xloc < 0 || yloc < 0 || zloc < 0)
					printf("%8d %13.6f %13.6f %13.6f\n", m + 1, rcmxm, rcmym, rcmzm);
				bin[xloc + yloc*nxbin + zloc*nxbin*nybin]++; 
				binxy[xloc + yloc*nxbin]++;
				binxz[xloc + zloc*nxbin]++;
				binyz[yloc + zloc*nybin]++;*/
				for (i = 0; i < nseg; i++){
					rxmi = rx[m*nseg+i];
					rymi = ry[m*nseg+i];
					rzmi = rz[m*nseg+i];
					rxmi -= sidex*roundf(rxmi/sidex);
					rymi -= sidey*roundf(rymi/sidey);
					rzmi -= sidez*roundf(rzmi/sidez);
					rxmi += sidex / 2.0;
					rymi += sidey / 2.0;
					rzmi += sidez / 2.0;
					xloc = int(rxmi / dx);
					yloc = int(rymi / dy);
					zloc = int(rzmi / dz);
					if(xloc == nxbin) xloc--;
					if(yloc == nybin) yloc--;
					if(zloc == nzbin) zloc--;
	
					if (xloc >= nxbin || yloc >= nybin || zloc >= nzbin)
						printf("%8d %8d %13.6f %13.6f %13.6f %6d %6d %6d\n", m + 1, i + 1, rxmi, rymi, rzmi, xloc, yloc, zloc);
					if (xloc < 0 || yloc < 0 || zloc < 0)
						printf("%8d %8d %13.6f %13.6f %13.6f %6d %6d %6d\n", m + 1, i + 1, rxmi, rymi, rzmi, xloc, yloc, zloc);
					bin[xloc + yloc*nxbin + zloc*nxbin*nybin]++; 
					binxy[xloc + yloc*nxbin]++;
					binxz[xloc + zloc*nxbin]++;
					binyz[yloc + zloc*nybin]++;
				}

			}

			for (j = 0; j < nybin; j++){
				for (i = 0; i < nxbin; i++){
					fprintf(Binning_result_xy, "%6d", binxy[i + j*nxbin]);

				}
				fprintf(Binning_result_xy, "\n");
			}

			for (j = 0; j < nzbin; j++){
				for (i = 0; i < nxbin; i++){
					fprintf(Binning_result_xz, "%6d", binxz[i + j*nxbin]);

				}
				fprintf(Binning_result_xz, "\n");
			}

			for (j = 0; j < nzbin; j++){
				for (i = 0; i < nybin; i++){
					fprintf(Binning_result_yz, "%6d", binyz[i + j*nybin]);

				}
				fprintf(Binning_result_yz, "\n");
			}
			for (k = 0; k < nzbin; k++){
				for (j = 0; j < nybin; j++){
					for (i = 0; i < nxbin; i++){
						fprintf(Binning_3d, "%2d", bin[i + j*nxbin + k*nxbin*nybin]);
					}
					fprintf(Binning_3d, "\n");
				}
			}


		}
	}

	free(rcmx); free(rcmy); free(rcmz); 
	free(rx); free(ry); free(rz); 
	free(binxz); free(binxy); free(binyz);

	fclose(Parameters); fclose(center_mass);
	fclose(rxfile); fclose(ryfile); fclose(rzfile);
	fclose(Binning_parameters);
	fclose(Binning_result_xy); fclose(Binning_result_xz);
	fclose(Binning_result_yz); fclose(Binning_3d);
	return 0;

}
