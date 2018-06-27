#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "zeroBin.h"
#include "binning.h"
#include "calcI.h"

using namespace std;

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *rxfile, *ryfile, *rzfile;
	FILE *Intfile, *INSInput; 
	
	// open files //
	Parameters = fopen("Parameters.in", "r");
	INSInput = fopen("INSinput.gen", "r");
	rxfile = fopen("rx.txt", "rb");
	ryfile = fopen("ry.txt", "rb");
	rzfile = fopen("rz.txt", "rb");
	Intfile = fopen("Intensity.txt", "w"); 
	
	int nfib, nseg, config_write, idum;
	float dt, strain, rp, sidex, sidey, sidez, dum;

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

	// Number of configurations
	int frames;
	frames = strain / dt / float(config_write) + 1; 

	// Read INSInput to determine simulation cases //
	int simcase;
	float dxRef, dyRef, dzRef; 	
	fscanf(INSInput, "%d", &simcase);
	fscanf(INSInput, "%*[^\n]%f", &dxRef); 
	fscanf(INSInput, "%*[^\n]%f", &dyRef); 
	fscanf(INSInput, "%*[^\n]%f", &dzRef); 

	// Based on simulation cases, determine simulation
	// box size at every frame
	float *Lx, *Ly, *Lz; 
	Lx = (float*)malloc(frames*sizeof(float)); 
	Ly = (float*)malloc(frames*sizeof(float));
	Lz = (float*)malloc(frames*sizeof(float));
	 
	// case 0 - basis - sidex does not change
	if (simcase == 0){
		for (int f = 0; f < frames; f++){
			Lx[f] = sidex; 
			Ly[f] = sidey; 
			Lz[f] = sidez; 
		}
	}
	// case 1 - redispersion - read sidex from Lbox.txt
	else {
		// Read box info
		FILE *BoxFile;
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

	float pi, volfrac, Vseg;
	float dx, dy, dz;
	int  nxbin, nybin, nzbin; 
	int nxbinMax, nybinMax, nzbinMax;
	pi = 3.141592654; 
	Vseg = 2.0*pi*rp;
	nxbinMax = int(sidex/dxRef);
	nybinMax = int(sidey/dyRef);
	nzbinMax = int(sidez/dzRef);

	float *rx, *ry, *rz, *bin, *intensity; 
	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));
	intensity = (float*)malloc(frames*sizeof(float)); 
	bin = (float*)malloc(nxbinMax*nybinMax*nzbinMax*sizeof(float));

	for (int step = 0; step < frames; step++){
		 
		// read data
		fread(&dum, sizeof(float), 1, rxfile);
		fread(&dum, sizeof(float), 1, ryfile);
		fread(&dum, sizeof(float), 1, rzfile);
		fread(rx, sizeof(float), nfib*nseg, rxfile);
		fread(ry, sizeof(float), nfib*nseg, ryfile);
		fread(rz, sizeof(float), nfib*nseg, rzfile);
		 
		// calculate bin dimension and constants
		nxbin = int(Lx[step]/dxRef);
		nybin = int(Ly[step]/dyRef);
		nzbin = int(Lz[step]/dzRef);
		dx = Lx[step]/float(nxbin);
		dy = Ly[step]/float(nybin);
		dz = Lz[step]/float(nzbin);
		volfrac = 2.0*nfib*nseg*rp*pi / (Lx[step]*Ly[step]*Lz[step]); 

		// zero variables 
		zeroBin(bin,nxbin,nybin,nzbin);
		 
		// binning
		binning(bin, nxbin, nybin, nzbin, dx, dy, dz, 
			rx, ry, rz, nfib, nseg,
			Lx[step], Ly[step], Lz[step]); 
		 
		// calculate intensity
		intensity[step] = calcI(bin, nxbin, nybin, nzbin, volfrac, 
					dx*dy*dz, Vseg);
		intensity[step] /= float(nxbin*nybin*nzbin)*volfrac*(1 - volfrac);
		fprintf(Intfile,"%10.4f %10.6f\n", float(step*config_write)*dt, intensity[step]);
	}
	
	free(rx); free(ry); free(rz); free(bin); 
	free(intensity);  
	free(Lx); free(Ly); free(Lz);

	fclose(Parameters); fclose(rxfile); fclose(ryfile); fclose(rzfile); 
	fclose(Intfile);  fclose(INSInput);   

	return 0; 

}
