#include <iostream>
#include <cmath>

using namespace std;

void fcn(double x, double *y,double *f) {
    f[0] = y[1];
    f[1] = y[3];
    f[2] = -y[0];
    f[3] = -y[2]-y[1];
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


int main()
{
	//double lambda;
	double eps[3], y[4];
	//double igrek1, igrek2, igrek3;
	//double delta1, delta2;
	//int i;
	//double GLOBAL_ERROR[3];
	
	//lambda = 18.956265591084957122802734375;
	//lambda = 0;
	eps[0]=1e-7;
	eps[1]=1e-9;
	eps[2]=1e-11;
	
	//Начальные данные для eps3  --------------------------------
	y[0]=0;
	y[1]=1;
	y[2]=0;
	y[3]=0;
	
	ddopri5(4,fcn,0,y,1.0e0,eps[2],1.0e0,0.5e0);
	
	
	return 0;
}

