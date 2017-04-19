#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "diploma_func/gauss_m.cpp"

int main(){
	double **A;
	double x[4], b[4];
	
	A = (double **) malloc(4 * sizeof(double *)); //копия матрицы
    for (int i = 0; i < 4; i++) A[i] = (double *) malloc(4 * sizeof(double));
	
	A[0][0] =1;		A[0][1] =1;		A[0][2] =2;		A[0][3] =2;
	
	
	A[1][0] =2;		A[1][1] =7;		A[1][2] =2;		A[1][3] =2;
	
	
	A[2][0] =2;		A[2][1] =2;		A[2][2] =23;		A[2][3] =3;
	
	
	A[3][0] =4;		A[3][1] =2;		A[3][2] =1;		A[3][3] =8;
	
	
	b[0] = 1;
	b[1] = 3;
	b[2] = 2;
	b[3] = 0;
	
	linear_sys(A, x, b, 4);
	
	for(int i=0; i<4; i++) {
		cout << x[i] << "\n";
	}
return 0;
}
