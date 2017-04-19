#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>

#define EPS_DDOPRI 1e-11
#define ens 1e-10 //epsilon newton stop
#define eps_gamma 2e-15

using namespace std;

//GLOBAL VARIABLES
double g = 9.8;
ofstream RHOout;
ofstream mout;
//----------------

#include "diploma_func/utilits.cpp"

//------------------------------------------
void fcn(double x, double *y, double *f) {
	double RHO, C1, C2, sinT, cosT;
	
	//------------------------------
	C1 = y[2]*y[3];
	C2 = y[2]*y[4]-g*y[5];
	RHO = sqrt(C1*C1+C2*C2);
	cosT = C1/RHO;
	sinT = C2/RHO;
	//------------------------------
	
	f[0] = y[2]*cosT;
	f[1] = y[2]*sinT;
	f[2] = -g*sinT;
	f[3] = 0;
	f[4] = 0;
	f[5] = -y[3]*cosT-y[4]*sinT;
	
	//RHOout << RHO << endl;
}



//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "diploma_func/ddopri5_brakhistakhrona.cpp"
#include "diploma_func/gauss_m.cpp"

void l(double *beta, double *res, double *y) {
	//double H;
	
	y[0] = 0; // x(0)
	y[1] = 0; // y(0)
	y[2] = 0; // v(0)
	y[3] = beta[0]; //px(0)
	y[4] = beta[1]; //py(0)
	//y[5] = beta[2]; //pv(0)
	y[5] = 1;

	//H = fabs(y[5]);
	ddopri5(6,fcn,0,y,beta[2],1.e-11,beta[2],0.5e0);
	
	res[0] = y[0] - 100;
	res[1] = y[1] + 1;
	res[2] = y[5];
	//res[3] = H - g;
}

void gradf(double *beta, double *dbeta, double *res, double *y) {
	int i, j, N=3;
	double beta_p[N], beta_m[N], res_p[N], res_m[N], res_minus[N], h=1.e-6;
	double **pro;
	
	
	pro=(double**)malloc(N*sizeof(double*));
	for(i=0; i<N; i++) {
		pro[i]=(double*)malloc(N*sizeof(double));
	}

	for(j=0; j<N; j++) {
		for(i=0; i<N; i++) {
			beta_p[i]=beta[i];
			beta_m[i]=beta[i];
		}
		
		beta_p[j]+=h;
		beta_m[j]-=h;
		
		l(beta_p, res_p, y);
		l(beta_m, res_m, y);
		
		for(i=0; i<N; i++) pro[i][j]=(res_p[i]-res_m[i])/(2*h);
	}
	
	//matrix_print(pro, N);
	matrix_print_f(pro, N);
	
	for(i=0; i<N; i++) res_minus[i]=-res[i];
	linear_sys(pro, dbeta, res_minus, N);
	
	for(i=0; i<N; i++) printf("dbeta: %.16f\n",dbeta[i]);
	//printf("TEST!");
	//matrix_print(pro2, N);
}

int NEWTON(double *beta, double *y) {
	int i, j, N=3;
	double res[N], res_w[N], beta_w[N], dbeta[N], gamma;
	bool flag;
	
	l(beta,res,y);
	
	for(j=0; j<15; j++) {
		//cout << "\n-------\nj = " << j << "\n-------\n";
	
		if(fabs(res[0])<ens && fabs(res[1])<ens && fabs(res[2])<ens) {
			cout << endl <<"Ended by " << j << " iteration" << endl; 
			return j;
		}
		gradf(beta,dbeta,res,y);
	
		gamma=1.0;
		flag=false;
		while(gamma>eps_gamma) {
			for(i=0; i<N; i++) 
				{
					beta_w[i]=beta[i]+gamma*dbeta[i];
					//printf("dbeta: %.16f\n",dbeta[i]);
				}
			l(beta_w, res_w, y);
			
			if(norm(res_w,N)<norm(res,N)) {flag=true; break;}
			gamma/=2.0;
		}
		if(flag==false) {cout << endl << "Broken on " << j << "'s iteration" << endl; return -1;}
		for(i=0; i<N; i++) {
			beta[i]=beta_w[i];
			res[i]=res_w[i];
		}
	}
cout << endl << "not enough iteration" << endl;
return -2;
}

int main() {
	double beta[3] = {1, 1, 2.35384};
	double y[6];
	int check;
	
	RHOout.open("RHO.txt");
	mout.open("matrix_pro.txt");
	
	cout << "0 for ddopri 5, 1 for NEWTON: ";
	cin >> check;

	if (check==0) {
			double res[4];
			l(beta,res,y);
	}
	else if (check==1) {
		NEWTON(beta,y);
		cout << "\nT = " << beta[2];
	}

return 0;
}
