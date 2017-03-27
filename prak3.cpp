#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>

#define EPS_DDOPRI 1e-11
#define ens 1e-10
#define eps_gamma 2e-15

using namespace std;

#include "prak3_func/utilits.cpp"

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
	double RHO, psi_xv, psi_yvpsi_v;
    
    psi_xv = y[3]*y[3]*y[2]*y[2];
    psi_yvpsi_v = (y[4]*y[2]-y[5])*(y[4]*y[2]-y[5]);
    RHO = sqrt(psi_xv+psi_yvpsi_v);
    
    //правая часть
    //f(v) = v^n, 	n = 1, 	=>	 f'(v) = (v)' = 1;
    f[0]=y[2]*y[2]*y[3]/RHO;									//x^{\cdot}
	f[1]=(y[2]*y[2]*y[4]-y[2]*y[5])/RHO;						//y^{\cdot}
	f[2]=-y[2]-(y[4]*y[2]-y[5])/RHO; 							//v^{\cdot}
	f[3]=0;														//psi_x^{\cdot}
	f[4]=0;														//psi_y^{\cdot}
	f[5]=y[5]-(y[2]*y[3]*y[3]+y[2]*y[4]*y[4]-y[4]*y[5])/RHO;	//psi_v^{\cdot}
}

//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "prak3_func/ddopri5_prak3.cpp"

void f_reserv(double *y, double *beta) {
	//double psi_xv, psi_yvpsi_v, RHO;

	y[0]=0; //x_0
	y[1]=0; //y_0
	y[2]=0; //v_0
	y[3]=beta[0];
	y[4]=beta[1];              
	y[5]=beta[2];
	ddopri5(6,fcn,0,y,beta[3],1.e-11,2.0e0,1e0);
}

void f(double *beta, double *res, double *y, double hres2) { //added double *y, double hres2
	double psi_xv, psi_yvpsi_v, RHO;

	y[0]=0; //x_0
	y[1]=10; //y_0
	y[2]=0; //v_0
	y[3]=beta[0];
	y[4]=beta[1];              
	y[5]=beta[2];
	ddopri5(6,fcn,0,y,beta[3],1.e-11,1.0e0,0.5e0);
	psi_xv = y[3]*y[3]*y[2]*y[2];
    psi_yvpsi_v = (y[4]*y[2]-y[5])*(y[4]*y[2]-y[5]);
    RHO = sqrt(psi_xv+psi_yvpsi_v); //RHO(T)
    
	res[0]=y[0]-1;//-0.1386315299798196; //x_T
	res[1]=y[1];//-8.8869461060126618; //y_T
	res[2]=y[2]-0.8583067592798990;//-0.8583067592798990; //v_T
	res[3]=RHO-y[5]*y[2]-1;	//f(v) = v
return;
}

void gradf(double *beta, double **pro, double *y, double hres2) { //added double *y, double hres2
	int i, j, N=3;
	double beta_p[N], beta_m[N], res_p[N], res_m[N], h=1.e-6;

	for(j=0; j<N; j++) {
		for(i=0; i<N; i++) {
			beta_p[i]=beta[i];
			beta_m[i]=beta[i];
		}
		
		beta_p[j]+=h;
		beta_m[j]-=h;
		
		f(beta_p, res_p, y, hres2);
		f(beta_m, res_m, y, hres2);
		
		for(i=0; i<N; i++) pro[i][j]=(res_p[i]-res_m[i])/(2*h);
	}
}

#include "prak3_func/gauss_m.cpp"


int NEWTON(double *beta, double *y, double hres2) { //added double *y, double hres2
	int i, k, flag, N=3;
	double res[N], res_w[N], beta_w[N], dbeta[N], gamma, **pro, res_m[N];
	
	pro=(double**)malloc(N*sizeof(double*));
	for(i=0; i<N; i++) pro[i]=(double*)malloc(N*sizeof(double));
	
	f(beta, res, y, hres2);
	
	for(k=0; k<15; k++) {
		if(fabs(res[0])<ens && fabs(res[1])<ens && fabs(res[2])<ens && fabs(res[3])<ens) {cout << endl <<"Ended by " << k << " iteration" << endl; return k;}
		gradf(beta, pro, y, hres2);
		for(i=0; i<N; i++) res_m[i]=-res[i];
		linear_sys(pro, dbeta, res_m, 4);
		
		gamma=1.0;
		flag=-1;
		
		while(gamma>sqrt(eps_gamma)) {
			for(i=0; i<N; i++) beta_w[i]=beta[i]+gamma*dbeta[i];
			f(beta_w, res_w, y, hres2);
			if(norm(res_w,N)<norm(res,N)) {flag=0; break;}
			gamma/=2.0;
		}
		if(flag==-1) {cout << endl << "Broken on " << k << "'s iteration" << endl; return -1;}
		for(i=0; i<N; i++) {
			beta[i]=beta_w[i];
			res[i]=res_w[i];
		}
	}
cout << endl << "not enough iteration" << endl;
return -2;
}

void mprpp(double *y, double *beta) { //метод продолжения решения по параметру
	double epsilon=0.1, hres2=0.0;
	int kol=0;
	
	NEWTON(beta, y, hres2);
	beta[3] = 2;
	hres2+=1.e-5;
	NEWTON(beta, y, hres2);
	
	while(fabs(y[5])>epsilon) { //неведома фигня, короче!
		if(y[5]<0) {
			hres2=-hres2;
			hres2-=1.e-5;
		}
		else hres2+=1.e-5;
		NEWTON(beta, y, hres2);
		kol++;
		printf("\n\n\n\t\tkol: %d", kol);
	}
	printf("\n\n\n\t\tkol: %d", kol);
	
	cout << "T=2" << endl;
		cout <<fixed<<setprecision(16)<< "\t\tx(2) = " << y[0] << endl;
		cout <<fixed<<setprecision(16)<< "\t\ty(2) = " << y[1] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tv(2) = " << y[2] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tpsi_x(2) = " << y[3] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tpsi_y(2) = " << y[4] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tpsi_v(2) = " << y[5] << endl;
		cout << "-------------------------------" << endl;
}

int main() {
	double hres2=1;
	double y[6], beta[4] = {1, 1, 1, 2};
	//f_reserv(y, beta);
	NEWTON(beta, y, hres2);
	//mprpp(y, beta);
	//printf("\n\n\n");
	//for(int i=0; i<4; i++) printf("beta[%d] = %.16f\n", i, beta[i]);
return 0;
}
