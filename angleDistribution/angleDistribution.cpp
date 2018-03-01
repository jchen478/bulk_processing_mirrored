#define _CRT_SECURE_NO_WARNINGS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

void elasticEnergy(float *px, float *py, float *pz,
	float *rx, float *ry, float *rz,
	float *q0, float *q1, float *q2, float *q3,
	float *R11, float *R12, float *R13,
	float *R21, float *R22, float *R23,
	float *R11eq, float *R12eq, float *R13eq,
	float *R21eq, float *R22eq, float *R23eq,
	float *R31eq, float *R32eq, float *R33eq,
	float kb, int nfibC_ID, int nseg,
	float *Ebend, float *Etwist, float *Eelastic,
	float *dphi, float *dtheta,
	float *phiAvg, float *thetaAvg,
	float *phiMax, float *thetaMax);

int main(){

	////////////////////////////////////////
	//        input and output files      //
	////////////////////////////////////////

	FILE *Parameters, *Equilibrium_Angles;
	FILE *rxfile, *ryfile, *rzfile;
	FILE *pxfile, *pyfile, *pzfile;
	FILE *q0file, *q1file, *q2file, *q3file;
	FILE *Eelastic_file, *Angle;
	Parameters = fopen("Parameters.in", "r");
	Eelastic_file = fopen("Eelastic.txt", "w");
	Angle = fopen("Angle.txt", "w");
	Equilibrium_Angles = fopen("Equilibrium_Angles.in", "r");
	rxfile = fopen("rx.txt", "rb");
	ryfile = fopen("ry.txt", "rb");
	rzfile = fopen("rz.txt", "rb");
	pxfile = fopen("px.txt", "rb");
	pyfile = fopen("py.txt", "rb");
	pzfile = fopen("pz.txt", "rb");
	q0file = fopen("q0.txt", "rb");
	q1file = fopen("q1.txt", "rb");
	q2file = fopen("q2.txt", "rb");
	q3file = fopen("q3.txt", "rb");

	int nfib, nseg, config_write;
	float dt, strain, kb, Eelastic, Ebend, Etwist, dum;
	int m, i, n, j, mi, nj, idum, idum1, idum2;
	// dummy variables

	// Read in Parameters.in //
	fscanf(Parameters, "%d", &nfib);
	fscanf(Parameters, "%*[^\n]%d", &nseg);
	fscanf(Parameters, "%*[^\n]%f", &dum);
	fscanf(Parameters, "%*[^\n]%f", &kb);
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
	fscanf(Parameters, "%*[^\n]%d", &config_write);
	fclose(Parameters);

	int step, nConfig;
	float *thetaeq, *phieq, theta, phi; 
	float thetaAvg, phiAvg, thetaStd, phiStd, thetaMax, phiMax; 
	float *R11eq, *R12eq, *R13eq, *R21eq, *R22eq, *R23eq;
	float *R11, *R12, *R13, *R21, *R22, *R23;
	float *R31eq, *R32eq, *R33eq;
	float *px, *py, *pz;
	float *rx, *ry, *rz;
	float *q0, *q1, *q2, *q3;
	float *dphi, *dtheta; 

	rx = (float*)malloc(nfib*nseg*sizeof(float));
	ry = (float*)malloc(nfib*nseg*sizeof(float));
	rz = (float*)malloc(nfib*nseg*sizeof(float));
	px = (float*)malloc(nfib*nseg*sizeof(float));
	py = (float*)malloc(nfib*nseg*sizeof(float));
	pz = (float*)malloc(nfib*nseg*sizeof(float));
	q0 = (float*)malloc(nfib*nseg*sizeof(float));
	q1 = (float*)malloc(nfib*nseg*sizeof(float));
	q2 = (float*)malloc(nfib*nseg*sizeof(float));
	q3 = (float*)malloc(nfib*nseg*sizeof(float));
	thetaeq = (float*)malloc(nfib*nseg*sizeof(float));
	phieq = (float*)malloc(nfib*nseg*sizeof(float));
	dphi = (float*)malloc(nfib*nseg*sizeof(float));
	dtheta = (float*)malloc(nfib*nseg*sizeof(float));


	R11 = (float*)malloc(nfib*nseg*sizeof(float));
	R12 = (float*)malloc(nfib*nseg*sizeof(float));
	R13 = (float*)malloc(nfib*nseg*sizeof(float));
	R21 = (float*)malloc(nfib*nseg*sizeof(float));
	R22 = (float*)malloc(nfib*nseg*sizeof(float));
	R23 = (float*)malloc(nfib*nseg*sizeof(float));
	R11eq = (float*)malloc(nfib*nseg*sizeof(float));
	R12eq = (float*)malloc(nfib*nseg*sizeof(float));
	R13eq = (float*)malloc(nfib*nseg*sizeof(float));
	R21eq = (float*)malloc(nfib*nseg*sizeof(float));
	R22eq = (float*)malloc(nfib*nseg*sizeof(float));
	R23eq = (float*)malloc(nfib*nseg*sizeof(float));
	R31eq = (float*)malloc(nfib*nseg*sizeof(float));
	R32eq = (float*)malloc(nfib*nseg*sizeof(float));
	R33eq = (float*)malloc(nfib*nseg*sizeof(float));

	nConfig = int(strain / dt / float(config_write)) + 1;

	// Read Equilibrium Angles //
	for (m = 0; m < nfib; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			if (i > 0){
				fscanf(Equilibrium_Angles, "%d %d %f %f", &idum1, &idum2, thetaeq + mi, phieq + mi);
				theta = thetaeq[mi];
				phi = phieq[mi];
				R11eq[mi] = cosf(theta)*cosf(phi);
				R12eq[mi] = cosf(theta)*sinf(phi);
				R13eq[mi] = -sinf(theta);
				R21eq[mi] = -sinf(phi);
				R22eq[mi] = cosf(phi);
				R23eq[mi] = 0.0;
				R31eq[mi] = sinf(theta)*cosf(phi);
				R32eq[mi] = sinf(theta)*sinf(phi);
				R33eq[mi] = cosf(theta);
			}
		}
	}
	fclose(Equilibrium_Angles);

	// read till last frame
	for (step = 0; step < nConfig; step++){

		// r
		fread(&dum, sizeof(float), 1, rxfile);
		fread(rx, sizeof(float), nfib*nseg, rxfile);
		fread(&dum, sizeof(float), 1, ryfile);
		fread(ry, sizeof(float), nfib*nseg, ryfile);
		fread(&dum, sizeof(float), 1, rzfile);
		fread(rz, sizeof(float), nfib*nseg, rzfile);

		// p
		fread(&dum, sizeof(float), 1, pxfile);
		fread(px, sizeof(float), nfib*nseg, pxfile);
		fread(&dum, sizeof(float), 1, pyfile);
		fread(py, sizeof(float), nfib*nseg, pyfile);
		fread(&dum, sizeof(float), 1, pzfile);
		fread(pz, sizeof(float), nfib*nseg, pzfile);

		// q
		fread(&dum, sizeof(float), 1, q0file);
		fread(q0, sizeof(float), nfib*nseg, q0file);
		fread(&dum, sizeof(float), 1, q1file);
		fread(q1, sizeof(float), nfib*nseg, q1file);
		fread(&dum, sizeof(float), 1, q2file);
		fread(q2, sizeof(float), nfib*nseg, q2file);
		fread(&dum, sizeof(float), 1, q3file);
		fread(q3, sizeof(float), nfib*nseg, q3file);
		
		// calculate elastic energy from angles
		thetaAvg = 0.0; 
		phiAvg = 0.0;
		thetaStd = 0.0; 
		phiStd = 0.0;
		thetaMax = 0.0; 
		phiMax = 0.0;
		elasticEnergy(px, py, pz,
			rx, ry, rz,
			q0, q1, q2, q3,
			R11, R12, R13, R21, R22, R23,
			R11eq, R12eq, R13eq, R21eq, R22eq,
			R23eq, R31eq, R32eq, R33eq,
			kb, nfib, nseg, &Ebend, &Etwist, &Eelastic, 
			dphi, dtheta, &phiAvg, &thetaAvg, &phiMax, &thetaMax);
		for (m = 0 ; m < nfib; m++){
			for (i = 1; i < nseg; i++){
				mi = m*nseg+i;
				thetaStd += (dtheta[mi] - thetaAvg)*(dtheta[mi] - thetaAvg); 
				phiStd += (dphi[mi] - phiAvg)*(dphi[mi] - phiAvg); 
			}
		}
		thetaStd /= float(nfib*(nseg-1)-1);	
		phiStd /= float(nfib*(nseg-1)-1); 
		thetaStd = sqrt(thetaStd);	
		phiStd = sqrt(phiStd);	
		fprintf(Angle,"%8.3f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", 
			float(step*config_write)*dt,thetaAvg, phiAvg, thetaStd, phiStd, thetaMax, phiMax);
		fprintf(Eelastic_file,"%8.3f %10.6f %10.6f %10.6f\n", float(step*config_write)*dt,Ebend, Etwist, Eelastic);
	}

	free(rx); free(ry); free(rz);
	free(px); free(py); free(pz);
	free(q0); free(q1); free(q2); free(q3);
	free(thetaeq); free(phieq);
	free(dphi); free(dtheta); 

	free(R11);
	free(R12);
	free(R13);
	free(R21);
	free(R22);
	free(R23);
	free(R11eq);
	free(R12eq);
	free(R13eq);
	free(R21eq);
	free(R22eq);
	free(R23eq);
	free(R31eq);
	free(R32eq);
	free(R33eq);

	fclose(Eelastic_file); fclose(Angle);
	fclose(rxfile); fclose(ryfile); fclose(rzfile);
	fclose(pxfile); fclose(pyfile); fclose(pzfile);
	fclose(q0file); fclose(q1file); fclose(q2file); fclose(q3file);

	return 0;
}

void elasticEnergy(float *px, float *py, float *pz,
	float *rx, float *ry, float *rz,
	float *q0, float *q1, float *q2, float *q3,
	float *R11, float *R12, float *R13,
	float *R21, float *R22, float *R23,
	float *R11eq, float *R12eq, float *R13eq,
	float *R21eq, float *R22eq, float *R23eq,
	float *R31eq, float *R32eq, float *R33eq,
	float kb, int nfibC_ID, int nseg,
	float *Ebend, float *Etwist, float *Eelastic, 
	float *dphi, float *dtheta,
	float *phiAvg, float *thetaAvg,
	float *phiMax, float *thetaMax){

	float cx, cy, cz, ang_theta, ang_phi;
	float zeqx, zeqy, zeqz, zdotz, dum;
	float yitx, yity, yitz, yieqx, yieqy, yieqz;
	float pos_theta, pos_phi;
	int m, i, mi;
	*Ebend = 0.0;
	*Etwist = 0.0;

	for (m = 0; m < nfibC_ID; m++){
		for (i = 0; i < nseg; i++){
			mi = m*nseg + i;
			R11[mi] = 2.0*(q0[mi] * q0[mi] + q1[mi] * q1[mi]) - 1.0;
			R12[mi] = 2.0*(q1[mi] * q2[mi] + q0[mi] * q3[mi]);
			R13[mi] = 2.0*(q1[mi] * q3[mi] - q0[mi] * q2[mi]);
			R21[mi] = 2.0*(q1[mi] * q2[mi] - q0[mi] * q3[mi]);
			R22[mi] = 2.0*(q0[mi] * q0[mi] + q2[mi] * q2[mi]) - 1.0;
			R23[mi] = 2.0*(q3[mi] * q2[mi] + q0[mi] * q1[mi]);
		}
	}

	for (m = 0; m < nfibC_ID; m++){
		for (i = 1; i < nseg; i++){
			ang_theta = 0.0;
			ang_phi = 0.0;
			mi = m*nseg + i;
			// calculations for theta
			zeqx = R11[mi - 1] * R31eq[mi] + R21[mi - 1] * R32eq[mi] + px[mi - 1] * R33eq[mi];
			zeqy = R12[mi - 1] * R31eq[mi] + R22[mi - 1] * R32eq[mi] + py[mi - 1] * R33eq[mi];
			zeqz = R13[mi - 1] * R31eq[mi] + R23[mi - 1] * R32eq[mi] + pz[mi - 1] * R33eq[mi];
			zdotz = px[mi] * zeqx + py[mi] * zeqy + pz[mi] * zeqz;

			if (zdotz <= 1.0 - 1.0E-6){
				if (zdotz < -1.0){
					zdotz = -1.0;
				}
				ang_theta = acosf(zdotz);
			}

			// calculations for phi
			cx = rx[mi] - rx[mi - 1];
			cy = ry[mi] - ry[mi - 1];
			cz = rz[mi] - rz[mi - 1];
			dum = sqrtf(cx*cx + cy*cy + cz*cz);
			cx = cx / dum;
			cy = cy / dum;
			cz = cz / dum;
			zeqx = R11[mi - 1] * R21eq[mi] + R21[mi - 1] * R22eq[mi] + px[mi - 1] * R23eq[mi];
			zeqy = R12[mi - 1] * R21eq[mi] + R22[mi - 1] * R22eq[mi] + py[mi - 1] * R23eq[mi];
			zeqz = R13[mi - 1] * R21eq[mi] + R23[mi - 1] * R22eq[mi] + pz[mi - 1] * R23eq[mi];
			yitx = R21[mi] - cx*(cx*R21[mi] + cy*R22[mi] + cz*R23[mi]);
			yity = R22[mi] - cy*(cx*R21[mi] + cy*R22[mi] + cz*R23[mi]);
			yitz = R23[mi] - cz*(cx*R21[mi] + cy*R22[mi] + cz*R23[mi]);
			dum = sqrtf(yitx*yitx + yity*yity + yitz*yitz);

			yitx = yitx / dum;
			yity = yity / dum;
			yitz = yitz / dum;
			yieqx = zeqx - cx*(cx*zeqx + cy*zeqy + cz*zeqz);
			yieqy = zeqy - cy*(cx*zeqx + cy*zeqy + cz*zeqz);
			yieqz = zeqz - cz*(cx*zeqx + cy*zeqy + cz*zeqz);
			dum = sqrtf(yieqx*yieqx + yieqy*yieqy + yieqz*yieqz);

			yieqx = yieqx / dum;
			yieqy = yieqy / dum;
			yieqz = yieqz / dum;
			zdotz = yitx*yieqx + yity*yieqy + yitz*yieqz;

			if (zdotz <= (1.0 - 1.0E-6)){
				if (zdotz < -1.0){
					zdotz = -1.0;
				}
				ang_phi = acosf(zdotz);
			}
			// calculate elastic energy			
			*Ebend = *Ebend + ang_theta*ang_theta;
			*Etwist = *Etwist + ang_phi*ang_phi;

			// store angle info
			pos_theta = fabsf(ang_theta); 
			pos_phi = fabsf(ang_phi); 
			dtheta[mi] = pos_theta; 
			dphi[mi] = pos_phi;
			*thetaAvg = *thetaAvg + pos_theta;
			*phiAvg = *phiAvg + pos_phi;
			if (pos_theta > *thetaMax) *thetaMax = pos_theta;
			if (pos_phi > *phiMax) *phiMax = pos_phi;
		}
	}

	*Ebend = kb* (*Ebend) / 2.0 / float(nfibC_ID);
	*Etwist = 0.67*kb * (*Etwist) / 2.0 / float (nfibC_ID);
	*Eelastic = *Ebend + *Etwist;
	*thetaAvg = *thetaAvg / float(4.0*nfibC_ID);
	*phiAvg = *phiAvg / float(4.0*nfibC_ID);
}

