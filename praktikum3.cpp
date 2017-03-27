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

void matrix_print(double **pro, int N) {
	int i, j;
	cout << "\n--------------------\n";
	for(i=0; i<N; i++) {
		for(j=0; j<N; j++) {
			cout << pro[i][j] << " ";
		}
		cout << endl;
	}
	cout << "--------------------\n";
}

double norm(double *vector, int N) {
	int i;
    double cur, max_W, sum2;

    max_W = 0.0;
    for(i=0; i<N; i++ ){
      cur = fabs(vector[i]);
      if(cur>max_W) max_W=cur;
    }
    if(max_W==0.0) return 0.0;

    sum2=0.0;
    for(i=0; i<N; i++ ){
      cur=vector[i]/max_W;
      sum2+=cur*cur;
    }
return max_W*sqrt(sum2);
}

void fcn(double x, double *y, double *f) {
	/*
	 * y[0] - x1
	 * y[1] - x2
	 * y[2] - p1
	 * y[3] - p2
	*/

	f[0] = y[1];			//x2
	f[1] = y[3];			//p2
	f[2] = -y[0];			//-x1
	f[3] = -y[2]-y[1];		//-p1-x2
}

//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "ddopri5_praktikum3.cpp"
#include "prak3_func/gauss_m.cpp"

void l(double *beta, double *res) {
	double y[4]={0, beta[0], beta[1], 0};
	
	ddopri5(4, fcn, 0.0e0, y, 1.0e0, EPS_DDOPRI, 1.0e0, 0.5e0);
	res[0]=y[0]-0;
	res[1]=y[1]-1;
}

void gradf(double *beta, double *dbeta, double *res) {
	int i, j, N=2;
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
		
		l(beta_p, res_p);
		l(beta_m, res_m);
		
		for(i=0; i<N; i++) pro[i][j]=(res_p[i]-res_m[i])/(2*h);
	}
	
	matrix_print(pro, N);
	
	for(i=0; i<N; i++) res_m[i]=-res[i];
	linear_sys(pro, dbeta, res_m, N);
	//kramer(pro,dbeta,res);
	
	
	printf("TEST!");
	
	matrix_print(pro2, N);
}

int NEWTON(double *beta) {
	int i, k, N=2;
	double res[N], res_w[N], beta_w[N], dbeta[N], gamma;
	bool flag;
	
	l(beta, res);
	
	for(k=0; k<15; k++) {
		cout << "\n-------\nk = " << k << "\n-------\n";
	
		if(fabs(res[0])<ens && fabs(res[1])<ens) {cout << endl <<"Ended by " << k << " iteration" << endl; return k;}
		gradf(beta,dbeta,res);
		//kramer(pro, dbeta, res);
		
		
		printf("\n\n**************lynear_sys***************\n\n");
		for(i=0; i<N; i++) printf("\ndbeta[%d] = %.16f", i, dbeta[i]);
		
		gamma=1.0;
		flag=false;
		while(gamma>eps_gamma) {
			for(i=0; i<N; i++) 
				{
					beta_w[i]=beta[i]+gamma*dbeta[i];
					printf("beta_w: %.16f\n",beta_w[i]);
				}
			l(beta_w, res_w);
			
			if(norm(res_w,N)<norm(res,N)) {flag=true; break;}
			gamma/=2.0;
		}
		if(flag==false) {cout << endl << "Broken on " << k << "'s iteration" << endl; return -1;}
		for(i=0; i<N; i++) {
			beta[i]=beta_w[i];
			res[i]=res_w[i];
		}
	}
cout << endl << "not enough iteration" << endl;
return -2;
}

int main() {
	double beta[2];
	//double y[4]={0, 0, 1, 0};
	printf("input beta: ");
	scanf("%lg %lg", &beta[0], &beta[1]);
	
	//double res[2];
	//double y[4] = {0, 0, beta[0], beta[1]};
	//ddopri5(4, fcn, 0.0e0, y, 1.0e0, 2e-15, 1.0e0, 0.5e0);
	//f(beta, res);
	NEWTON(beta);
return 0;
}
