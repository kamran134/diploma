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

double hres_x, hres_y;

#include "diploma_func/utilits.cpp"

void MY_PRINT(double *y, double t) {
	cout << endl << t << " : " << endl;
	cout << y[0] << " , " << y[1] << " , " << y[2] << " , " << y[3] << endl;
return;
}

//------------------------------------------
void fcn(double x, double *y, double *f) {
	double RHO;
	double g=9.8;

    RHO = sqrt(y[6]*y[6]+y[7]*y[7]);
	
	f[0] = y[2];
	f[1] = y[3];
	
	if(RHO<=1e-12) {
		RHO = sqrt(y[4]*y[4]+y[5]*y[5]);
		f[2] = (g/2.)*(y[4]/RHO);
		f[3] = (g/2.) + (g/2.)*(y[5]/RHO);
	}
	else {
		f[2] = (g/2.)*(y[6]/RHO);
		f[3] = (g/2.) + (g/2.)*(y[7]/RHO);
	}
	f[4] = 0;
	f[5] = 0;
	f[6] = -y[4];
	f[7] = -y[5];
	
	//f[2] = (g*y[6])/(2*RHO);
	//f[3] = g/2. + (g*y[7])/(2*RHO);
}

double H(double *y) {
	double g=9.8, RHO;
	RHO = sqrt(y[6]*y[6]+y[7]*y[7]);
	return y[4]*y[2]+y[5]*y[3]+(g*RHO)/2. + (y[7]*g)/2.;
}

//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "diploma_func/ddopri5_brakhistakhrona.cpp"
#include "diploma_func/gauss_m.cpp"

void l(double *beta, double *res, double *y) {
	double RHO2;
	double px_0, py_0, pvx_0, pvy_0;
	double g=9.8;

	
	//---------------------
	RHO2 = 2./(g*beta[1]-g*beta[1]*cos(2.*beta[0]));
	px_0 = -RHO2*sin(2.*beta[0]);
	py_0 = RHO2*cos(2.*beta[0]);
	pvx_0 = px_0*beta[1];
	pvy_0 = py_0*beta[1];
	//---------------------

	y[0]=0; //x_0
	y[1]=0; //y_0
	y[2]=0; //vx_0
	y[3]=0; //vy_0
	y[4]=px_0; //px_0
	y[5]=py_0; //py_0
	y[6]=pvx_0; //pvx_0
	y[7]=pvy_0; //pvy_0
	
	//cout << "\npvx_0/px_0 = " << pvx_0/px_0 <<"\npvy_0/py_0 = " << pvy_0/py_0 << endl;
	
	ddopri5(8,fcn,0,y,beta[1],1.e-11,1.0e0,0.5e0);
	res[0]=y[0] - hres_x;
	res[1]=y[1] - hres_y;
	
}

#include "diploma_func/brakh_newtoon.cpp"

int mprpp(double *beta, double *y) {
	double hres_x_new;
	hres_x_new = hres_x;
	hres_x = 5;
	//hres_y = 7;
	double h = 0.1;
	NEWTON(beta,y);
	while(hres_x<hres_x_new) {
		hres_x+=h;
		NEWTON(beta,y);
		//cout << foo << endl;
		cout << beta[0] << " , " << beta[1] << endl;
	}

return 1;
}

int main() {
	double beta[2] = {0.9, 2}; //beta[0] = \teta, beta[1] = T;
	//double beta[2] = {-1.49486, 15.7541};
	double y[8];
	double res[2];
	
	cout << "Input x and y: ";
	cin >> hres_x >> hres_y;
	/*
	hres_x = 93;
	hres_y = 7;
	*/
	//l(beta,res,y);
	//NEWTON(beta,y);
	mprpp(beta,y);
	cout << "\nT = " << beta[1] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tx(" << beta[1] << ") = " << y[0] << endl;
		cout <<fixed<<setprecision(16)<< "\t\ty(" << beta[1] << ") = " << y[1] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tv_x(" << beta[1] << ") = " << y[2] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tv_y(" << beta[1] << ") = " << y[3] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tp_x(" << beta[1] << ") = " << y[4] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tp_y(" << beta[1] << ") = " << y[5] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tp_v_x(" << beta[1] << ") = " << y[6] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tp_v_y(" << beta[1] << ") = " << y[7] << endl;
	
	cout << "\n res[0] = " << res[0] << "\n res[1] = " << res[1];

return 0;
}
