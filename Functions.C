#include <iostream>    
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include "Fun.h"

const double PI=acos(-1.);
const double e=exp(1.);

/**********************************************************************************************************************/
double fRand(double fMin, double fMax){//Random number between fMin and fMax
	double f=(double)rand()/RAND_MAX;
 	return fMin+f*(fMax-fMin);}
 	
int fNRand(int fMin, int fMax){//Random number between fMin and fMax
	int f=rand()/RAND_MAX;
 	return fMin+f*(fMax-fMin);}

int **CREA_Int(int n, int m){//Creates matrix nxm (int entries)
		int **X;
		int i, j;
		X=new int *[n];
		for(i=0; i<n; i++) X[i]=new int [m];
		for(i=0; i<n; i++){
			for(j=0; j<m; j++) X[i][j]=0;}
		return X;}

double **CREA(int n, int m){//Creates matrix nxm (double entries)
		double **X;
		int i, j;
		X=new double *[n];
		for(i=0; i<n; i++) X[i]=new double [m];
		for(i=0; i<n; i++){
			for(j=0; j<m; j++) X[i][j]=0;}
		return X;}
		
int **RandMtx(int n, int m){//Creates a random matrix nxm (int entries)
		int **X;
		int i, j, aux;
		X=new int *[n];
		for(i=0; i<n; i++) X[i]=new int [m];
		for(i=0; i<n; i++){
			for(j=0; j<m; j++){
				aux=rand()%2;
				X[i][j]=pow((-1),aux);}}
		return X;}

double **DRandMtx(int n, int m){//Creates a random matrix nxm (double entries)
		double **X, aux;
		int i, j;
		X=new double *[n];
		for(i=0; i<n; i++) X[i]=new double [m];
		for(i=0; i<n; i++){
			for(j=0; j<m; j++){
				aux=fRand(0,1);
				if(aux<=0.5) X[i][j]=-1;
				else X[i][j]=1;}}
		return X;}
		
double **BoundCond(double **S, int N){//For periodic boundary conditions S(n+1)=S(1)
	double **X;
	int n, m;
	n=sqrt(N);
	m=n+2;
	
	X=CREA(m,m);
		
	for(int i=0; i<n; i++){//Boundary conditions
		for(int j=0; j<n; j++){
			X[i+1][j+1]=S[i][j];//Real position for spins 
			X[0][j+1]=S[n-1][j];//if(i==0) bottom side
			X[n+1][j+1]=S[0][j];//if(i==n+1) top side
			X[j+1][0]=S[j][n-1];//if(j==0) left side
			X[j+1][n+1]=S[j][0];}}//if(j==n+1) right side

	return X;}

double TotEnergy(double **X, int N, double J, double *M){//Calculates total initial energy and magnetization
	double aux1=0, aux2=0, E;
	int n=sqrt(N);
	for(int i=1; i<=n; i++){
		for(int j=1; j<=n; j++){
			aux1+=X[i][j]*(X[i][j+1]+X[i+1][j]);
			aux2+=X[i][j];}}
	E=-J*aux1;
	*M=aux2;
	return E;
	}

double **Metropolis(double **X, int N, double J, double beta, double E0, double M0, double *E, double *M){//Metropolis algortihm
	double OldS, NewS, sum, DE, DM, EF, MF, p, Boltz;
	int n=sqrt(N);
	int m=n+2;
	
	for(int k=1; k<=N; k++){

	int j=(int) fRand(1, n);//Select column j randomly for the k-spin
	int i=(int) fRand(1, n);//Select row i randomly for the k-spin
	
	//EF=E0;//System energy not changes
	//MF=M0;//System magnetization not changes
	
	OldS=X[i][j];//Old spin value
	NewS=-OldS;//Split the spin to change its value
	sum=X[i+1][j]+X[i][j+1]+X[i-1][j]+X[i][j-1];//Neighbors energy
	DM=(NewS-OldS);//Changes in magnetization of spin k
	DE=J*(OldS-NewS)*sum;//Changes in energy of spin k
	//printf("\nk=%d i=%d j=%d OldS=%.lf NewS=%.lf sum=%.lf DE=%.lf DM=%.lf\n", k, i, j, OldS, NewS, sum, DE, DM);
	
	p=fRand(0,1);
	Boltz=exp(-beta*DE);
	//printf("p=%.3lf, DE=%.lf, exp=%.3lf", p, DE, Boltz);
	
	if(DE<=0.0 || p<=Boltz){
		X[i][j]=NewS;//Accepts the changes and actualize de system configuration
		E0+=DE;//New system energy
		M0+=DM;//New system magnetization
		if(i==1) X[n+1][j]=NewS;//if(i==1) bottom side
		if(i==n) X[0][j]=NewS;//if(i==n) top side
		if(j==1) X[i][n+1]=NewS;//if(j==1) left side
		if(j==n) X[i][0]=NewS;}//if(j==n) right side
	else{
		X[i][j]=OldS;//Rejects the changes and back to initial value for k-spin
		}
	}
	
	*E=E0;
	*M=M0;

	return X;}

void InitVect(int n, double *v){//Puts in all entries of v the value 0
	for(int i=0; i<n; i++) v[i]=0;}

void InitMtxInt(int n, int m, int **M){//Puts in all entries of M (int) the value 0
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++) M[i][j]=0;}}

void InitMtx(int n, int m, double **M){//Puts in all entries of M (double) the value 0
	for(int i=0; i<n; i++){
		for(int j=0; j<m; j++) M[i][j]=0;}}
