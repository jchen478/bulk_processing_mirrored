#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using namespace std;

int main(){

	FILE *Parameters, *Eelastic, *Elastic_Results;

	// open files //
	Parameters = fopen("Parameters.in", "r");
	Eelastic = fopen("Eelastic.txt", "r");
	Elastic_Results = fopen("Elastic_Results.txt", "w");

	int idum, elastic_write, nstep, N;
	int total_steps, step, start_step;
	float dum, start_elastic, strain,dt;
	float *Eb, *Et, *Etot, time;
	float Ebsum, Etsum, Etotsum;
	float Ebsd, Etsd, Etotsd;
	float dum1, dum2, dum3;

	fscanf(Parameters, "%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dt);
	fscanf(Parameters, "%*[^\n]%f", &strain);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, " %f", &dum);
	fscanf(Parameters, " %f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
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
	fscanf(Parameters, "%*[^\n]%d", &idum);
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
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%d", &idum);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%d", &elastic_write);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &start_elastic);
	fclose(Parameters);

	total_steps = strain / dt / float(elastic_write); 
	start_step = start_elastic / dt / float(elastic_write) - 2; 

	Eb = (float*)malloc(total_steps*sizeof(float));
	Et = (float*)malloc(total_steps*sizeof(float));
	Etot = (float*)malloc(total_steps*sizeof(float));

	// calculation for average 
	N = 0; Ebsum = 0.0; Etsum = 0.0; Etotsum = 0.0;
	for (step = 0; step < total_steps; step++){

		fscanf(Eelastic, "%f %f %f %f", &time, &dum1, &dum2, &dum3);

		if (step > start_step){
			Eb[N] = dum1;
			Et[N] = dum2;
			Etot[N] = dum3; 
			Ebsum += dum1;
			Etsum += dum2;
			Etotsum += dum3;
			N++;	
		}
	}
	fclose(Eelastic);

	Ebsum /= float(N);
	Etsum /= float(N);
	Etotsum /= float(N);

	fprintf(Elastic_Results,"Result for elastic energy\n");	
	fprintf(Elastic_Results,"-------------------------------------------------\n");	
	fprintf(Elastic_Results," Average\n");	
	fprintf(Elastic_Results,"-------------------------------------------------\n");	
	fprintf(Elastic_Results,"%10.6f\n%10.6f\n%10.6f\n", Ebsum, Etsum, Etotsum);	
	
	// calculation for standard deviation
	Ebsd = 0.0; Etsd = 0.0; Etotsd = 0.0;
	for (nstep = 0; nstep < N; nstep++){
		Ebsd += (Eb[nstep] - Ebsum)*(Eb[nstep] - Ebsum);
		Etsd += (Et[nstep] - Etsum)*(Et[nstep] - Etsum);
		Etotsd += (Etot[nstep] - Etotsum)*(Etot[nstep] - Etotsum);
	}

	Ebsd = sqrt(Ebsd / float(N-1));
	Etsd = sqrt(Etsd / float(N-1));
	Etotsd = sqrt(Etotsd / float(N-1));

	fprintf(Elastic_Results,"-------------------------------------------------\n");	
	fprintf(Elastic_Results," Standard deviation\n");	
	fprintf(Elastic_Results,"-------------------------------------------------\n");	
	fprintf(Elastic_Results,"%10.6f\n%10.6f\n%10.6f\n", Ebsd, Etsd, Etotsd);	
	
	// standard error of the mean
	
	fprintf(Elastic_Results,"-------------------------------------------------\n");	
	fprintf(Elastic_Results," Standard error of the mean\n");	
	fprintf(Elastic_Results,"-------------------------------------------------\n");	
	fprintf(Elastic_Results,"%10.6f\n%10.6f\n%10.6f\n", Ebsd/sqrt(float(N)), Etsd/sqrt(float(N)), Etotsd/sqrt(float(N)));	

	fclose(Elastic_Results);
	free(Eb); free(Et); free(Etot);

	return 0;

}
