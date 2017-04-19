#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <iomanip>

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

#include "prak3_func/gauss_m.cpp"

int main() {
	double **a;
	double b[3] = {1, 2, 3}, x[3];
	
	a=(double**)malloc(3*sizeof(double*));
	for(int i=0; i<3; i++) a[i]=(double*)malloc(3*sizeof(double));
	
	a[0][0]=5; a[0][1]=4; a[0][2]=3;
	a[1][0]=2; a[1][1]=5; a[1][2]=4;
	a[2][0]=0; a[2][1]=2; a[2][2]=5;
	
	linear_sys(a, x, b, 3);
	
return 0;
}
