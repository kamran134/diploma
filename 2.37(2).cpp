#include <iostream>
#include <cmath>

using namespace std;

void fcn(double x, double *y,double *f) {
    f[0]=y[1]; //y_2
	f[1]=-x*y[0]*y[2]; //x = t
	f[2]=0;
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


//Поиск минимальной лямбды
double lambdaSearch(double eps, double *y, double (*ddopri5)(int, void (*fcn)(double, double*, double*), double, double*, double, double, double, double)) {
	double hl;
	hl = 1.0;
		
	while(fabs(y[0])>1.e-10) {
		if(y[0]>0) y[2]+=hl;
		else if(y[0]<0) {hl = hl/2; y[2]-=hl;}
		
		y[0]=0;
		y[1]=1;
		
		ddopri5(3,fcn,0,y,1.0e0,eps,1.0e0,0.5e0);
		cout << endl << "<------------------------------>" << endl;
	}
	return y[2];
}


int main() {
	double lambda;
	double y[3], Y[3], eps[3];
	double igrek1, igrek2, igrek3;
	double delta1, delta2;
	int i;
	double GLOBAL_ERROR[3];
	
	//lambda = 18.956265591084957122802734375;
	lambda = 0;
	eps[0]=1e-7;
	eps[1]=1e-9;
	eps[2]=1e-11;
	
	//Начальные данные для eps3  --------------------------------
	y[0]=0;
	y[1]=1;
	y[2]=lambda;
	
	ddopri5(3,fcn,0,y,1.0e0,eps[2],1.0e0,0.5e0);
	lambda = lambdaSearch(eps[2], y, ddopri5);
	printf("\nMinimalnaya lambda: %.40f\n", lambda);
	//============================================
	
	//В точке t = 0.0260000000000000
	cout << "t = 0.26" << endl;
	for(i=0; i<3; i++) {
		y[0]=0;
		y[1]=1;
		y[2]=lambda;
		GLOBAL_ERROR[i] = ddopri5(3,fcn,0,y,0.026,eps[i],1.0e0,0.5e0);
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
		y[0]=0;
		y[1]=1;
		y[2]=lambda;
		GLOBAL_ERROR[i] = ddopri5(3,fcn,0,y,0.5000128147791552,eps[i],1.0e0,0.5e0);
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
		y[0]=0;
		y[1]=1;
		y[2]=lambda;
		GLOBAL_ERROR[i] = ddopri5(3,fcn,0,y,1.0e0,eps[i],1.0e0,0.5e0);
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

