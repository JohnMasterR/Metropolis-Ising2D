void Graph1(FILE *file, int R, int N, double J, int nT, double Tmin, double Tmax);//Creates script for makes a picture of network
double fRand(double fMin, double fMax);//Random double number between fMin and fMax
int fNRand(int fMin, int fMax);//Random int number between fMin and fMax
int **CREA_Int(int n, int m);//Creates matrix nxm (int entries)
double **CREA(int n, int m);//Creates matrix nxm
int **RandMtx(int n, int m);//Creates a random matrix nxm (int entries)
double **DRandMtx(int n, int m);//Creates a random matrix nxm (double entries)
double **BoundCond(double **S, int N);//For periodic boundary conditions S(n+1)=S(1)
double TotEnergy(double **X, int N, double J, double *M);//Calculates total initial energy and magnetization
double **Metropolis(double **S, int N, double J, double beta, double E0, double M0, double *E, double *M);//Metropolis algortihm
void InitVect(int n, double *v);//Puts in all entries of v the value 0
void InitMtxInt(int n, int m, int **M);//Puts in all entries of M (int) the value 0
void InitMtx(int n, int m, double **M);//Puts in all entries of M (double) the value 0
