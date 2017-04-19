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
double g=9.8;
ofstream fout;
//----------------

#include "diploma_func/utilits.cpp"

//------------------------------------------
void fcn(double x, double *y, double *f) {
	double RHO;
	double a_, b_;
	
	a_ = y[2]*y[2]*y[3]*y[3];
	b_ = (y[2]*y[4]-g*y[5])*(y[2]*y[4]-g*y[5]);
	RHO = sqrt(a_+b_);
	
	f[0] = y[2]*y[2]*y[3]/RHO;
	f[1] = y[2]*(y[4]*y[2]-g*y[5])/RHO;
	f[2] = -g*y[2]*y[3]/RHO;
	f[3] = 0;
	f[4] = 0;
	f[5] = -(y[2]*y[3]*y[3]+y[2]*y[4]*y[4]-g*y[4]*y[5])/RHO;
	
	fout << "\na = " << a_ << "\tb = " << b_ << endl;
	fout << "-------------------------------------"
}



//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "diploma_func/ddopri5_brakhistakhrona.cpp"
#include "diploma_func/gauss_m.cpp"

void l(double *beta, double *res, double *y) {
	double H;

	y[0] = 0; //x(0)
	y[1] = 0; //y(0)
	y[2] = 0; //v(0)
	y[3] = beta[0]; //px(0)
	y[4] = beta[1]; //py(0)
	y[5] = beta[2]; //pv(0)
	
	H = g*fabs(y[5]);//H(0)
	
	ddopri5(6,fcn,0,y,beta[3],1.e-11,1.0e0,0.5e0);
  
    res[0]=y[0]-5; //0.0727930740001; //-0.07279307323(2)1;
	res[1]=y[1]+7; //-6.0309331715847536;//T=5//-8.8869461060126618; //T=2//y_T
	res[2]=y[2];
	res[3]=H-1;
	
}

void gradf(double *beta, double *dbeta, double *res, double *y) {
	int i, j, N=4;
	double beta_p[N], beta_m[N], res_p[N], res_m[N], h=1.e-6;
	double **pro, **pro2;

	pro=(double**)malloc(N*sizeof(double*));
	pro2=(double**)malloc(N*sizeof(double*));
	for(i=0; i<N; i++) {
		pro[i]=(double*)malloc(N*sizeof(double));
		pro2[i]=(double*)malloc(N*sizeof(double));
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
	
	matrix_print(pro, N);
	
	for(i=0; i<N; i++) res_m[i]=-res[i];
	linear_sys(pro, dbeta, res_m, N);

	//printf("TEST!");
	//matrix_print(pro2, N);
}

int NEWTON(double *beta, double *y) {
	int i, j, N=4;
	double res[N], res_w[N], beta_w[N], dbeta[N], gamma;
	bool flag;
	
	l(beta,res,y);
	
	for(j=0; j<15; j++) {
		//cout << "\n-------\nj = " << j << "\n-------\n";
	
		if(fabs(res[0])<ens && fabs(res[1])<ens && fabs(res[2])<ens && fabs(res[3])<ens) {
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
					//printf("beta_w: %.16f\n",beta_w[i]);
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
	double beta[4] = {1, 1, 1, 1};
	double y[6];
	int check;
	
	fout.open("RHO.txt");
	
	cout << "0 for ddopri 5, 1 for NEWTON: ";
	cin >> check;

	if (check==0) {
			double res[4];
			l(beta,res,y);
	}
	else if (check==1) {
		NEWTON(beta,y);
		cout << "\nT = " << beta[4];
	}

return 0;
}
