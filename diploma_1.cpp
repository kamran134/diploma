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

double hres;

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

//resistance function
double fr(double k, double v) {
return k*v*v;
}

double f_r(double k, double v) {
return 2*k*v;
}

void fcn(double x, double *y, double *f) {
	double RHO, psi_xv, psi_yvpsi_v;
	double k;
	
	psi_xv = y[3]*y[3]*y[2]*y[2];
    psi_yvpsi_v = (y[4]*y[2]-y[5])*(y[4]*y[2]-y[5]);
    RHO = sqrt(psi_xv+psi_yvpsi_v);
    
    k=0;
    
    //правая часть
    //f(v) = kv^2, 	=>	 f'(v) = 2kv;
    f[0]=y[2]*y[2]*y[3]/RHO;									//x^{\cdot}
	f[1]=(y[2]*y[2]*y[4]-y[2]*y[5])/RHO;						//y^{\cdot}
	f[2]=-fr(k,y[2])-(y[4]*y[2]-y[5])/RHO; 							//v^{\cdot}
	f[3]=0;														//psi_x^{\cdot}
	f[4]=0;														//psi_y^{\cdot}
	f[5]=f_r(k, y[2])*y[5]-(y[2]*y[3]*y[3]+y[2]*y[4]*y[4]-y[4]*y[5])/RHO;	//psi_v^{\cdot}
}

//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "ddopri5_praktikum3.cpp"
#include "prak3_func/gauss_m.cpp"

void l(double *beta, double *res, double *y) {
	double psi_xv, psi_yvpsi_v, RHO;
	//double y[6];
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
    
    res[0]=y[0]-0; //0.0727930740001; //-0.07279307323(2)1;
	res[1]=y[1]-1; //-6.0309331715847536;//T=5//-8.8869461060126618; //T=2//y_T
	res[2]=y[2]-0; //-0.9927347427368552;//T=5//-0.8583067592798990; //T=2//v_T
	res[3]=RHO-y[5]*y[2]-1;	//f(v) = v
/*
	res[0]=y[0]-0.0727930732377216; //0.0727930740001; //-0.07279307323(2)1;
	res[1]=y[1]-6.0309331715847536; //-6.0309331715847536;//T=5//-8.8869461060126618; //T=2//y_T
	res[2]=y[2]-0.9927347427368552; //-0.9927347427368552;//T=5//-0.8583067592798990; //T=2//v_T
	res[3]=RHO-y[5]*y[2]-1;	//f(v) = v
*/
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
	//kramer(pro,dbeta,res);
	
	/*
	pro2[0][0] = sqrt(sqrt(5)/2-0.5)*sinh(sqrt(sqrt(5)/2-0.5))/sqrt(5)+sqrt(sqrt(5)/2+0.5)*sin(sqrt(sqrt(5)/2+0.5))/sqrt(5);
	pro2[0][1] = -sqrt(sqrt(5)/2+0.5)*sinh(sqrt(sqrt(5)/2-0.5))/sqrt(5)+sqrt(sqrt(5)/2-0.5)*sin(sqrt(sqrt(5)/2+0.5))/sqrt(5);
	pro2[1][0] =(sqrt(5)/2-0.5)*cosh(sqrt(sqrt(5)/2-0.5))/sqrt(5)+(sqrt(5)/2+0.5)*cos(sqrt(sqrt(5)/2+0.5))/sqrt(5) ;
	pro2[1][1] =-cosh(sqrt(sqrt(5)/2-0.5))/sqrt(5)+cos(sqrt(sqrt(5)/2+0.5))/sqrt(5) ;
	*/
	
	//printf("TEST!");
	
	//matrix_print(pro2, N);
}

int NEWTON(double *beta, double *y) {
	int i, k, N=4;
	double res[N], res_w[N], beta_w[N], dbeta[N], gamma;
	bool flag;
	
	l(beta, res, y);
	
	for(k=0; k<15; k++) {
		cout << "\n-------\nk = " << k << "\n-------\n";
	
		if(fabs(res[0])<ens && fabs(res[1])<ens && fabs(res[2])<ens && fabs(res[3])<ens) {cout << endl <<"Ended by " << k << " iteration" << endl; return k;}
		gradf(beta,dbeta,res,y);
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
			l(beta_w, res_w, y);
			
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

void pvn(void) { //попадание в 0
	
}

int main() {
	double beta[4] = {1, 1, 1, 5}, y[6];
	int kol=0;
	//double res[4];
	NEWTON(beta,y);
	//l(beta,res,y);
	printf("\n\n\nPSI_v(T) = %.16f\nkol = %d", y[5], kol);
	printf("\n\n--------\n%.16f\n%.16f\n%.16f\n%.16f\n%.16f\n", y[0],y[1],y[2],y[3],y[4]);

return 0;
}
