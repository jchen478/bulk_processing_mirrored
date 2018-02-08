#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
using namespace std; 

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *center_mass, *pair_dist;

	// open files //
	Parameters = fopen("Parameters.in", "r");
	center_mass = fopen("center_mass.txt", "rb");
	pair_dist = fopen("g.txt", "w");

	if (Parameters == NULL || center_mass == NULL || pair_dist == NULL) {
		perror("Error");
	}

	////////////////////////////////////////
	//               variables            //
	////////////////////////////////////////
	int nfib, nseg, config_write;
	int nBin, nstep, step, nstart, bin;
	float dt, strain, dr, rp;
	float *g, *rcmx, *rcmy, *rcmz;
	float sidex, sidey, sidez; 
	float corx, cory, corz; 
	float dx, dy, dz, rnew, rold, r, L; 
	float start_pair, rho, nid;
	int idum, m, n;
	float dum;

	float pi = 3.141592654;

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rp);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &config_write);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &start_pair);

printf("start_pair %10.6f\n", start_pair); 

	dr = 1;
	L = 2.0*rp*nseg; 
	nBin = int(L / dr); 
	nstep = strain / dt / config_write + 1;
	nstart = start_pair / dt / config_write; 
	rho = float(nfib) / (sidex*sidey*sidez); 

	printf("rho = %10f\n", rho); 

	rcmx = (float*)malloc(nfib*sizeof(float)); 
	rcmy = (float*)malloc(nfib*sizeof(float));	
	rcmz = (float*)malloc(nfib*sizeof(float));
	g = (float*)calloc(nBin, sizeof(float));

	for (step = 0; step < nstep; step++){
		fread(&dum, sizeof(float), 1, center_mass); 
		fread(rcmx, sizeof(float), nfib, center_mass); 
		fread(rcmy, sizeof(float), nfib, center_mass); 
		fread(rcmz, sizeof(float), nfib, center_mass); 
		
		if (step >= nstart){
			// calculate pair distribution function
			for (m = 0; m < nfib - 1; m++){
				for (n = m + 1; n < nfib; n++){
					dx = rcmx[n] - rcmx[m];
					dy = rcmy[n] - rcmy[m];
					dz = rcmz[n] - rcmz[m];
					corx = roundf(dx / sidex);
					cory = roundf(dy / sidey);
					corz = roundf(dz / sidez);
					dx -= corx*sidex; 
					dy -= cory*sidey; 
					dz -= corz*sidez; 
					r = sqrtf(dx*dx + dy*dy + dz*dz); 
					if (r < L){
						bin = int(r / dr);
						if (bin >= nBin){
							printf("bin %d nBin %d\ndx dy dz %6f %6f %6f r %6f L %6f dr %6f\n", bin, nBin, dx, dy, dz, r, L, dr);
							getchar();
						}
						g[bin] += 2.0;
					}
				}
			}
		}
	}

	
	rold = 0.0; 

	for (m = 0; m < nBin; m++){
		rnew = dr *float(m + 1); 
		nid = 4.0 / 3.0*pi*(rnew*rnew*rnew - rold*rold*rold)*rho; 
		rold = rnew; 
		g[m] /= nfib*nid*float(nstep - nstart); 
		if( (float(m)+0.5)*dr > sidex / 2.0) break;
		fprintf(pair_dist, "%6.4f %10f\n", (float(m) + 0.5)*dr / L, g[m]);
	}



	free(g);
	free(rcmx);
	free(rcmy);
	free(rcmz);
	fclose(Parameters);
	fclose(center_mass);
	fclose(pair_dist);


	return 0;
}
