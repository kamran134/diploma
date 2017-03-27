#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>

using namespace std;

#include "gauss_m.cpp"
int main(){

double M[3][3];
double b[3], x[3];
int i,j;


for(i=0; i<3; i++) {
	for(j=0; j<3; j++) {
		scanf("%lf", &M[i][j]);
	}
	scanf("%lg", &b[i]);
}

linear_sys(M, x, b, 3);

for(i=0; i<3; i++) printf("%.16f\n", x[i]);

return 0;
}
