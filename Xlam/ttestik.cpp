#include <iostream>
#include <cmath>

using namespace std;

void fcn(double x, double *y,double *f) {
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

double dsign(double a, double b) {
	if (b<0) { a=fabs(a)*(-1.0); return a; }
	if (b>0) return fabs(a);
	return 0.0;
}

double max_root(double x, double lambda) {
	return max(-0.5+0.5*lambda*x, 0.5-0.5*lambda*x);
}

//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "ddopri5.cpp"


int main() {
	double beta[4] = {1, 1, 1, 5};
	double y[6], Y[3], eps[3];
	double igrek1, igrek2, igrek3;
	double delta1, delta2;
	int i;
	double GLOBAL_ERROR[3];
	
	
	eps[0]=1e-7;
	eps[1]=1e-9;
	eps[2]=1e-11;
	
	
	//В точке t = 0.0260000000000000
	cout << "t = 0.26" << endl;
	for(i=0; i<3; i++) {
		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=0; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[i] = ddopri5(6,fcn,0,y,0.026,eps[i],1.0e0,0.5e0);
		Y[i] = y[0];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[i], GLOBAL_ERROR[i]);
	}
	cout << "**********************" << endl;
	cout << "eps = 1e-7 \t : \t y1 = " << Y[0] << endl;
	cout << "eps = 1e-9 \t : \t y2 = " << Y[1] << endl;
	cout << "eps = 1e-11 \t : \t y3 = " << Y[2] << endl;
	cout << "**********************" << endl;
	igrek1 = Y[0] - Y[1];
	igrek2 = Y[1] - Y[2];
	igrek3 = igrek1/igrek2;
	delta1 = GLOBAL_ERROR[0]/GLOBAL_ERROR[1];
	delta2 = GLOBAL_ERROR[1]/GLOBAL_ERROR[2];
	printf("y1 - y2 = %.16f\ny2 - y3 = %.16f\n(y1-y2)/(y2-y3) = %.16f\n",igrek1,igrek2,igrek3);
	printf("d1 / d2 = %.16f\nd2/d3 = %.16f\n",delta1,delta2);
	cout << "============================================================" << endl << endl;
	
	//В точке t = 0.5
	cout << "t = 0.5" << endl;
	for(i=0; i<3; i++) {
		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=0; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[i] = ddopri5(6,fcn,0,y,2.5,eps[i],1.0e0,0.5e0);
		Y[i] = y[0];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[i], GLOBAL_ERROR[i]);
	}
	cout << "**********************" << endl;
	cout << "eps = 1e-7 \t : \t y1 = " << Y[0] << endl;
	cout << "eps = 1e-9 \t : \t y2 = " << Y[1] << endl;
	cout << "eps = 1e-11 \t : \t y3 = " << Y[2] << endl;
	cout << "**********************" << endl;
	igrek1 = Y[0] - Y[1];
	igrek2 = Y[1] - Y[2];
	igrek3 = igrek1/igrek2;
	delta1 = GLOBAL_ERROR[0]/GLOBAL_ERROR[1];
	delta2 = GLOBAL_ERROR[1]/GLOBAL_ERROR[2];
	printf("y1 - y2 = %.16f\ny2 - y3 = %.16f\n(y1-y2)/(y2-y3) = %.16f\n",igrek1,igrek2,igrek3);
	printf("d1 / d2 = %.16f\nd2/d3 = %.16f\n",delta1,delta2);
	cout << "============================================================" << endl << endl;
	
	//В точке t = 1
	cout << "t = 1" << endl;
	for(i=0; i<3; i++) {
		y[0]=0; //x_0
		y[1]=10; //y_0
		y[2]=0; //v_0
		y[3]=beta[0];
		y[4]=beta[1];
		y[5]=beta[2];
		GLOBAL_ERROR[i] = ddopri5(6,fcn,0,y,5.0,eps[i],1.0e0,0.5e0);
		Y[i] = y[0];
		printf("\nepsilon: %lg -------- GLOBAL ERROR: %.16f\n", eps[i], GLOBAL_ERROR[i]);
	}
	cout << "**********************" << endl;
	cout << "eps = 1e-7 \t : \t y1 = " << Y[0] << endl;
	cout << "eps = 1e-9 \t : \t y2 = " << Y[1] << endl;
	cout << "eps = 1e-11 \t : \t y3 = " << Y[2] << endl;
	cout << "**********************" << endl;
	igrek1 = Y[0] - Y[1];
	igrek2 = Y[1] - Y[2];
	igrek3 = igrek1/igrek2;
	delta1 = GLOBAL_ERROR[0]/GLOBAL_ERROR[1];
	delta2 = GLOBAL_ERROR[1]/GLOBAL_ERROR[2];
	printf("y1 - y2 = %.16f\ny2 - y3 = %.16f\n(y1-y2)/(y2-y3) = %.16f\n",igrek1,igrek2,igrek3);
	printf("d1 / d2 = %.16f\nd2/d3 = %.16f\n",delta1,delta2);
	cout << "============================================================" << endl << endl;
	
	
	return 0;
}

