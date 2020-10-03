#include <iostream>    
using namespace std; 

#include <stdio.h>    
#include <stdlib.h>
#include <math.h>

void Graph1(FILE *file, int R, int N, double J, int nT, double Tmin, double Tmax){//Creates script for makes a picture of network
char str[50];
double dT=(Tmax-Tmin)/(double)nT;

sprintf(str,"ScriptN%d.gp", N);
file=fopen(str,"w+");

for(int i=0; i<=nT; i++){
double T=dT*i+Tmin;
fprintf(file,"#To make graph\n");
fprintf(file,"reset\n");
fprintf(file,"set terminal x11 0\n");
fprintf(file,"set term postscript eps size 8cm,8cm enhanced color font \",18\"\n");
fprintf(file,"set output 'Ising2D(T_%.3lf)I.eps'\n", T);
fprintf(file,"set title 'Initial Configuration (R%d,N%d,T%.3lf)'\n", R, N, T);
fprintf(file,"unset key\n");
fprintf(file,"set autoscale xfix\n");
fprintf(file,"set autoscale yfix\n");
fprintf(file,"set autoscale cbfix\n");
fprintf(file,"set palette defined (-1 \"blue\", 1 \"green\")\n");
fprintf(file,"pl 'Ising2D_Plot0.dat' matrix with image notitle \n");
fprintf(file,"#\n");
fprintf(file,"\n");

fprintf(file,"#To make graph\n");
fprintf(file,"reset\n");
fprintf(file,"set terminal x11 0\n");
fprintf(file,"set term postscript eps size 8cm,8cm enhanced color font \",18\"\n");
fprintf(file,"unset key\n");
fprintf(file,"set output 'Ising2D(T_%.3lf)F.eps'\n", T);
fprintf(file,"#\n");
fprintf(file,"set title 'Final Configuration (R%d,N%d,T%.3lf)'\n", R, N, T);
fprintf(file,"set autoscale xfix\n");
fprintf(file,"set autoscale yfix\n");
fprintf(file,"set autoscale cbfix\n");
fprintf(file,"set palette defined (-1 \"blue\", 1 \"green\")\n");
fprintf(file,"pl 'Ising2D_Plot1.dat' matrix with image notitle \n");
fprintf(file,"#\n");
fprintf(file,"\n");

fprintf(file,"#To make E(t)\n");
fprintf(file,"reset\n");
fprintf(file,"set terminal x11 0\n");
fprintf(file,"set term postscript eps size 10cm,7cm enhanced color font \",18\"\n");
fprintf(file,"set output 'Energy(T_%.3lf).eps'\n", T);
fprintf(file,"set title 'Energy vs Time (N=%d spins)'\n", N);
//fprintf(file,"set autoscale\n");
fprintf(file,"set xrange [0:%d]\n", R);
fprintf(file,"set yrange [%.lf:]\n", -2*J);
//fprintf(file,"set format y '%%.3e'\n");
fprintf(file,"set key right bottom\n");
fprintf(file,"set pointsize 0.5\n");
fprintf(file,"set xlabel 'Time (Runs)'\n");
fprintf(file,"set ylabel '{/Symbol e}=E/(JN)'\n");
fprintf(file,"pl \"EnergyPlot.dat\" u 1:%d w points t \"T=%.3lf K\" pt 7\n", 2*(i+1), T);
fprintf(file,"\n");

fprintf(file,"#To make M(t)\n");
fprintf(file,"reset\n");
fprintf(file,"set terminal x11 0\n");
fprintf(file,"set term postscript eps size 10cm,7cm enhanced color font \",18\"\n");
fprintf(file,"set output 'Magnetization(T_%.3lf).eps'\n", T);
fprintf(file,"set title 'Magnetization vs Time (N=%d spins)'\n", N);
//fprintf(file,"set autoscale\n");
fprintf(file,"set xrange [0:%d]\n", R);
fprintf(file,"set yrange [-1:1]\n");
//fprintf(file,"set format y '%%.3e'\n");
fprintf(file,"set key right bottom\n");
fprintf(file,"set pointsize 0.5\n");
fprintf(file,"set xlabel 'Time (Runs)'\n");
fprintf(file,"set ylabel 'm=M/N'\n");
fprintf(file,"pl \"EnergyPlot.dat\" u 1:%d w points t \"T=%.3lf K\" pt 7\n", 2*(i+1)+1, T);
fprintf(file,"\n");}

fprintf(file,"#To make E(T)\n");
fprintf(file,"reset\n");
fprintf(file,"set terminal x11 0\n");
fprintf(file,"set term postscript eps size 10cm,7cm enhanced color font \",18\"\n");
fprintf(file,"set output 'MeanE_vs_T.eps'\n");
fprintf(file,"set title 'Mean Energy vs Temperature (N=%d)\n",N);
fprintf(file,"set autoscale\n");
//fprintf(file,"set yrange [-1:1]\n");
//fprintf(file,"set format y '%%.3e'\n");
fprintf(file,"unset key\n");
fprintf(file,"set pointsize 0.5\n");
fprintf(file,"set xlabel 'Temperature (K)'\n");
fprintf(file,"set ylabel '<{/Symbol e}>'\n");
fprintf(file,"pl \"MeanEnergyPlot.dat\" u 1:2 w lp pt 7\n");
fprintf(file,"\n");

fprintf(file,"#To make M(T)\n");
fprintf(file,"reset\n");
fprintf(file,"set terminal x11 0\n");
fprintf(file,"set term postscript eps size 10cm,7cm enhanced color font \",18\"\n");
fprintf(file,"set output 'MeanM_vs_T.eps'\n");
fprintf(file,"set title 'Mean Magnetization vs Temperature (N=%d)\n",N);
fprintf(file,"set autoscale\n");
//fprintf(file,"set yrange [-1:1]\n");
//fprintf(file,"set format y '%%.3e'\n");
fprintf(file,"unset key\n");
fprintf(file,"set pointsize 0.5\n");
fprintf(file,"set xlabel 'Temperature (K)'\n");
fprintf(file,"set ylabel '<m>'\n");
fprintf(file,"pl \"MeanEnergyPlot.dat\" u 1:3 w lp pt 7\n");
fprintf(file,"\n");
fclose(file);}
















