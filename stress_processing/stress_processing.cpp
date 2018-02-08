
//
//  main.cpp
//  flexfric
//
//  Created by Jing-Yao Chen on 6/6/17.
//  Copyright (c) 2016 Jing-Yao Chen. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h> /* malloc, calloc, free, exit */
#include <cmath>

using namespace std;

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////
	FILE *Parameters, *Stress_tensor, *Stress_results, *Number_of_Contacts;

	// open files //
	Parameters = fopen("Parameters.in", "r");
	Stress_tensor = fopen("Stress_tensor.txt", "r");
	Number_of_Contacts = fopen("Number_of_Contacts.txt", "r");
	Stress_results = fopen("Stress_results.txt", "w");

	int nfib, nseg, config_write, maxCon, maxGr, maxBin;
	int fac, bdimx, bdimy, bdimz, fiberPerBlock, contact_write;
	float rp, kb, mu_stat, mu_kin;
	float contact_cutoff, rep_cutoff, over_cut;
	float dt, strain, sidex, sidey, sidez;
	float fstar, fact, Astar, decatt, elf, delta_rx;
	float duxdx, duydx, duzdx, duxdy, duydy;
	float duzdy, duxdz, duydz, duzdz;
	float fraction_rp, dx, dy, dz;
	int nfibGrid, nfibBlock;
	int blasGrid, blasBlock;
	int stress_write;
	float dum, start_stress, start_contacts;
	// nfib - number of fibers
	// nseg - number of segments
	// config_write - number of time steps between configuration writes
	// maxCon - max fiber contacting one fiber
	// maxGr - max number of fibers in a group
	// maxBin - max number of fibers in a bin
	// fac - factor of fluid velocity
	// bdimx, bdimy, bdimz - number of bins to look for contacts
	// fiberPerBlock - number of fiber per block
	// rp - segment aspect ratio
	// kb - bending constant
	// mu_stat/mu_kin - static/kinetic coefficient of friction
	// contact_cutoff - cutoff for contacts
	// rep_cutoff - cutoff distance for calculating repulsive forces
	// over_cut - limit to the overlapping of two fibers
	// dt,strain - time step, total time of run
	// sidex, sidey, sidez - simulation box dimensions (x and y directions)
	// fstar,fact - force prefactor, force expfonential factor (decay length)
	// Astar - Prefactor of the attractive force
	// decatt - decay length of the attractive force
	// elf - electric field factor
	// delta_rx - shift variable for sliding periodic images
	// duxdx... - velocity gradient tensor
	// fraction_rp - effective aspect ratio factor
	// dx, dy, dz - size of cell
	// nfibGrid, nfibBlock - grid and block size with total nfib threads
	// blasGrid, blasBlock - grid and block dimension for kernels involving
	//					     cublas functions

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &rp);
	fscanf(Parameters, "%*[^\n]%f", &kb);
	fscanf(Parameters, "%*[^\n]%f", &mu_stat);
	fscanf(Parameters, "%*[^\n]%f", &mu_kin);
	fscanf(Parameters, "%*[^\n]%f", &contact_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &rep_cutoff);
	fscanf(Parameters, "%*[^\n]%f", &over_cut);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &sidex);
	fscanf(Parameters, " %f", &sidey);
	fscanf(Parameters, " %f", &sidez);
	fscanf(Parameters, "%*[^\n]%f", &fraction_rp);
	fscanf(Parameters, "%*[^\n]%d", &config_write);
	fscanf(Parameters, "%*[^\n]%d", &contact_write);
	fscanf(Parameters, "%*[^\n]%f", &fstar);
	fscanf(Parameters, "%*[^\n]%f", &fact);
	fscanf(Parameters, "%*[^\n]%f", &Astar);
	fscanf(Parameters, "%*[^\n]%f", &decatt);
	fscanf(Parameters, "%*[^\n]%f", &delta_rx);
	fscanf(Parameters, "%*[^\n]%f", &duxdx);
	fscanf(Parameters, "%*[^\n]%f", &duydx);
	fscanf(Parameters, "%*[^\n]%f", &duzdx);
	fscanf(Parameters, "%*[^\n]%f", &duxdy);
	fscanf(Parameters, "%*[^\n]%f", &duydy);
	fscanf(Parameters, "%*[^\n]%f", &duzdy);
	fscanf(Parameters, "%*[^\n]%f", &duxdz);
	fscanf(Parameters, "%*[^\n]%f", &duydz);
	fscanf(Parameters, "%*[^\n]%f", &duzdz);
	fscanf(Parameters, "%*[^\n]%d", &fac);
	fscanf(Parameters, "%*[^\n]%f", &elf);
	fscanf(Parameters, "%*[^\n]%f", &dx);
	fscanf(Parameters, "%*[^\n]%f", &dy);
	fscanf(Parameters, "%*[^\n]%f", &dz);
	fscanf(Parameters, "%*[^\n]%d", &bdimx);
	fscanf(Parameters, "%*[^\n]%d", &bdimy);
	fscanf(Parameters, "%*[^\n]%d", &bdimz);
	fscanf(Parameters, "%*[^\n]%d", &fiberPerBlock);
	fscanf(Parameters, "%*[^\n]%d", &maxCon);
	fscanf(Parameters, "%*[^\n]%d", &maxGr);
	fscanf(Parameters, "%*[^\n]%d", &maxBin);
	fscanf(Parameters, "%*[^\n]%d", &nfibGrid);
	fscanf(Parameters, "%*[^\n]%d", &nfibBlock);
	fscanf(Parameters, "%*[^\n]%d", &blasGrid);
	fscanf(Parameters, "%*[^\n]%d", &blasBlock);
	fscanf(Parameters, "%*[^\n]%f", &start_stress);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &stress_write);
	fscanf(Parameters, "%*[^\n]%f", &start_contacts);

	fclose(Parameters); 
	
	float pi = 3.14159265;

	int total_steps, step, start_step;
	int N;
	float nL3, volfrac; 
	float time, Stress11, Stress12, Stress13;
	float Stress22, Stress23, Stress33; 

	float Stress11_sum, Stress12_sum, Stress13_sum;
	float Stress22_sum, Stress23_sum, Stress33_sum;

	float Stress11_sd, Stress12_sd, Stress13_sd;
	float Stress22_sd, Stress23_sd, Stress33_sd;

	nL3 = float(nfib)*powf((float(2.0)*float(nseg)*rp), 3.0) / (sidex*sidey*sidez);
	volfrac = float(nfib)*float(nseg)*pi*(float(2.0)*rp) / (sidex*sidey*sidez);

	total_steps = strain / dt / float(stress_write); 
	start_step = start_stress / dt / float(stress_write) - 2; 


	Stress11_sum = 0.0; Stress12_sum = 0.0; Stress13_sum = 0.0;
	Stress22_sum = 0.0; Stress23_sum = 0.0; Stress33_sum = 0.0;

	Stress11_sd = 0.0; Stress12_sd = 0.0; Stress13_sd = 0.0;
	Stress22_sd = 0.0; Stress23_sd = 0.0; Stress33_sd = 0.0;

	N = 0; 

	for (step = 0; step < total_steps; step++){

		fscanf(Stress_tensor, "%f %f %f %f %f %f %f", &time, &Stress11, &Stress12, &Stress13, &Stress22, &Stress23, &Stress33);

		if (step > start_step){
			Stress11_sum += Stress11;
			Stress12_sum += Stress12;
			Stress13_sum += Stress13;
			Stress22_sum += Stress22;
			Stress23_sum += Stress23;
			Stress33_sum += Stress33;
			N++; 
		}
	}
	
	fclose(Stress_tensor);
	
	Stress11_sum /= float(N);
	Stress12_sum /= float(N);
	Stress13_sum /= float(N);
	Stress22_sum /= float(N);
	Stress23_sum /= float(N);
	Stress33_sum /= float(N);
	
	Stress_tensor = fopen("Stress_tensor.txt", "r");

	N = 0; 
	for (step = 0; step < total_steps; step++){

		fscanf(Stress_tensor, "%f %f %f %f %f %f %f", &time, &Stress11, &Stress12, &Stress13, &Stress22, &Stress23, &Stress33);

		if (step > start_step){
			Stress11_sd += (Stress11 - Stress11_sum)*(Stress11 - Stress11_sum);
			Stress12_sd += (Stress12 - Stress12_sum)*(Stress12 - Stress12_sum);
			Stress13_sd += (Stress13 - Stress13_sum)*(Stress13 - Stress13_sum);
			Stress22_sd += (Stress22 - Stress22_sum)*(Stress22 - Stress22_sum);
			Stress23_sd += (Stress23 - Stress23_sum)*(Stress23 - Stress23_sum);
			Stress33_sd += (Stress33 - Stress33_sum)*(Stress33 - Stress33_sum);
			N++;
		}
	}

	fclose(Stress_tensor);

	Stress11_sd /= float(N-1);
	Stress12_sd /= float(N-1);
	Stress13_sd /= float(N-1);
	Stress22_sd /= float(N-1);
	Stress23_sd /= float(N-1);
	Stress33_sd /= float(N-1);

	Stress11_sd = sqrtf(Stress11_sd); 
	Stress12_sd = sqrtf(Stress12_sd);
	Stress13_sd = sqrtf(Stress13_sd);
	Stress22_sd = sqrtf(Stress22_sd);
	Stress23_sd = sqrtf(Stress23_sd);
	Stress33_sd = sqrtf(Stress33_sd);

	fprintf(Stress_results, "Stress calculation results:\n"); 
	fprintf(Stress_results, "------------------------------------\n");
	fprintf(Stress_results, "Number of fiber / segments: %4d %4d\n", nfib, nseg);
	fprintf(Stress_results, "Aspect ratio of fiber:      %6.2f\n", float(nseg)*rp);
	fprintf(Stress_results, "Box size length:            %6.2f %6.2f %6.2f\n", sidex, sidey, sidez); 
	fprintf(Stress_results, "Coefficients of friction:   %6.2f %6.2f\n", mu_stat, mu_kin); 
	fprintf(Stress_results, "Bending constant:           %6.2f\n", kb);
	fprintf(Stress_results, "Concentration, nL3:         %6.2f\n", nL3);
	fprintf(Stress_results, "Volume fraction:            %10.6f\n", volfrac);

	fprintf(Stress_results, "Overall particle stress tensor:\n");
	fprintf(Stress_results, "%6.3f\n", Stress11_sum);
	fprintf(Stress_results, "%6.3f\n", Stress12_sum);
	fprintf(Stress_results, "%6.3f\n", Stress13_sum);
	fprintf(Stress_results, "%6.3f\n", Stress22_sum);
	fprintf(Stress_results, "%6.3f\n", Stress23_sum);
	fprintf(Stress_results, "%6.3f\n", Stress33_sum);

	fprintf(Stress_results, "First and second normal difference:\n");
	fprintf(Stress_results, "%6.3f\n", Stress11_sum - Stress33_sum);
	fprintf(Stress_results, "%6.3f\n", Stress33_sum - Stress22_sum);

	fprintf(Stress_results, "Standard deviation:\n");
	fprintf(Stress_results, "%6.3f\n", Stress11_sd);
	fprintf(Stress_results, "%6.3f\n", Stress12_sd);
	fprintf(Stress_results, "%6.3f\n", Stress13_sd);
	fprintf(Stress_results, "%6.3f\n", Stress22_sd);
	fprintf(Stress_results, "%6.3f\n", Stress23_sd);
	fprintf(Stress_results, "%6.3f\n", Stress33_sd);

	fprintf(Stress_results, "Standard error of the mean:\n");
	fprintf(Stress_results, "%6.3f\n", Stress11_sd / sqrt(float(N)));
	fprintf(Stress_results, "%6.3f\n", Stress12_sd / sqrt(float(N)));
	fprintf(Stress_results, "%6.3f\n", Stress13_sd / sqrt(float(N)));
	fprintf(Stress_results, "%6.3f\n", Stress22_sd / sqrt(float(N)));
	fprintf(Stress_results, "%6.3f\n", Stress23_sd / sqrt(float(N)));
	fprintf(Stress_results, "%6.3f\n", Stress33_sd / sqrt(float(N)));

	int num_groups, contacts, total_contacts, overs; 
	float contacts_norm, contact_sd;

	// Contact statistics
	total_steps = strain / dt / contact_write; 
	start_step = start_contacts / dt / contact_write-2; 

	total_contacts = 0; 
	N = 0;
	contact_sd = 0.0;  

	for (step = 0; step < total_steps; step++){

		fscanf(Number_of_Contacts, "%f %d %d %f %d", &time, &num_groups, &contacts, &dum, &overs); 

		if (step > start_step){
			total_contacts += contacts; 
			N++;
		}
	}
	
	contacts_norm = float(total_contacts) / float(N*nfib); 
	
	for (step = 0; step < total_steps; step++){

		fscanf(Number_of_Contacts, "%f %d %d %f %d", &time, &num_groups, &contacts, &dum, &overs); 

		if (step > start_step){
			contact_sd += (float(contacts)/float(nfib)-contacts_norm)*(float(contacts)/float(nfib)-contacts_norm); 
		}
	}
	fclose(Number_of_Contacts); 

	contact_sd = sqrt(contact_sd/float(N));

	fprintf(Stress_results, "Number of contacts per fiber\n");
	fprintf(Stress_results, "%6.3f\n", contacts_norm);
	fprintf(Stress_results, "%6.3f\n", contact_sd);
	fprintf(Stress_results, "%6.3f\n", contact_sd / sqrt(float(N)));

	fclose(Stress_results);

	return 0;

}
