/*********************************************************************************************/
/*                               To compile this file                                        */
/*        g++ -o Ising.exe Ising2D.C Functions.C Scripts.C Files.C; time ./Ising.exe         */
/*                                                                                           */
/*********************************************************************************************/
#include <iostream>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Fun.h"

/*********************************************************************************************/
int main(){

srand(time(NULL));//Random Number seeds	randomize();
FILE *DT1, *DT2, *DT3, *DT4, *DT5;
char str[100];

int n, m, N, M, R, nT;
double **S, **SBound, **energy, **magnet;
double beta, kB, Tmin, Tmax, T, dT, J, hbar, MeanE, MeanM, E_T, M_T;

/******PARAMETERS******/

n=40;                                   //Spins by side of square lattice
R=10000;                                  //Total iterations
m=n+2;                                   //For take into account the boundary conditions 
nT=100;                                  //Points for temperature interval   

Tmin=0.075;                                //Initial temperature (K)
Tmax=6;                                  //Final temperature (K)
kB=1;//1.3806488E-23;     						   //Boltzmann constant (J/K)
J=1;                                     //Binding constant between spins (J)

dT=(Tmax-Tmin)/(double)nT;               //Setp increases for temperature (K)
N=pow(n,2);                              //Total number of spins
M=pow(n,2);                              //Total number of spins for boundary conditions

/******FILE CREATION******/
sprintf(str, "Ising2D_Plot0.dat");//To plot the initial system configuration
DT1=fopen(str,"w+");
sprintf(str, "Ising2D_Plot1.dat");//To plot the initial system configuration
DT2=fopen(str,"w+");
sprintf(str, "EnergyPlot.dat");//To plot the energy and magnetization vs time (runs)
DT3=fopen(str,"w");
sprintf(str, "MeanEnergyPlot.dat");//To plot the energy and magnetization vs time (runs)
DT4=fopen(str,"w");

S=CREA(n,n);//Creates a matrix nxm (entries = 0.0)
SBound=CREA(m,m);//Creates a matrix nxm (entries = 0.0)
S=DRandMtx(n,n);//Creates a random matrix nxm (int entries), this is the spin model
energy=CREA(R+1,nT+1);//Creates a matrix nxm (entries = 0.0) for plot energy vs iterations
magnet=CREA(R+1,nT+1);//Creates a matrix nxm (entries = 0.0) for plot magnetization vs iterations

SBound=BoundCond(S, N);//Periodic boundary conditions

for(int i=1; i<=n; i++){//To plots the spin initial configuration
	for(int j=1; j<=n; j++) fprintf(DT1, "%.lf ", SBound[i][j]);
		fprintf(DT1, "\n");}
		
fprintf(DT4, "#  T   |  <E>  |  <M>  |\n");

/******************************************For all T value********************************************************/

for(int l=0; l<=nT; l++){/*********************Starts loop in T***************************************************/
	
	energy[0][l]=TotEnergy(SBound, N, J, &magnet[0][l]);//Calculates total initial energy and magnetization
	
	T=dT*l+Tmin;
	beta=1.0/(kB*T);

	E_T=0;//For average energy after all iterations
	M_T=0;//For average magnetization after all iterations
	
	for(int i=0; i<R; i++){/*********************Starts loop in R***************************************************/
		SBound=Metropolis(SBound, N, J, beta, energy[i][l], magnet[i][l], &energy[i+1][l], &magnet[i+1][l]);
		E_T+=energy[i+1][l];
		M_T+=magnet[i+1][l];}/*********************Ends loop in R*****************************************************/
		
	MeanE=E_T/((double)R*N);//Calculates averages
	MeanM=M_T/((double)R*N);//Calculates averages
	
	fprintf(DT4, "%.3lf %.3lf %.3lf\n", T, MeanE, MeanM);}/*********************Ends loop in T***********************/

for(int i=1; i<=n; i++){//Prints the spines configuration without boundary conditions
	for(int j=1; j<=n; j++) fprintf(DT2, "%.lf ", SBound[i][j]);
		fprintf(DT2, "\n");}

fprintf(DT3, "#R|");
for(int j=0; j<=nT; j++){
	T=dT*j+Tmin;
	fprintf(DT3, "    T=%.3lf    |", T);}
fprintf(DT3, "\n");
for(int i=0; i<=R; i++){//Prints the energy and magnetization values vs total iterations
	fprintf(DT3, "%d ", i);
	for(int j=0; j<=nT; j++) fprintf(DT3, "%.3lf %.3lf ", energy[i][j]/(double)N, magnet[i][j]/(double)N);
		fprintf(DT3, "\n");}

/******CLOSE ALL FILES******/
fclose(DT1);
fclose(DT2);
fclose(DT3);
fclose(DT4);

/******MAKE THE PLOTS******/
Graph1(DT5, R, N, J, nT, Tmin, Tmax);
sprintf(str, "gnuplot ScriptN%d.gp", N);
system(str);

return 0;
}
