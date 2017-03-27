#include <iostream>
#include <cmath>
#include <iomanip>
//#define eps 2.e-15
using namespace std;

void matrix_print(double **pro, int N) {
	int i, j;
	cout << endl;
	for(j=0; j<N; j++) {
		for(i=0; i<N; i++) {
			cout << pro[i][j] << " ";
		}
		cout << endl;
	}
}

double norm(double *vector, int N) {
	int i;
    float cur, max_W, sum2;

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

void masswap(double **a, double *b, int n, int i, int j) {
    for (int k = 0; k < n; k++) swap(a[i][k], a[j][k]);
    swap(b[i], b[j]);
}

void gauss(double **a, double *x, double *b, int n) {
    int i, j, k, mxi;
    double mx;
    for (i=0; i<n; i++) {
        mx=a[i][i];
        mxi=i;
        for (j=i+1; j<n; j++)
            if (a[j][i]>mx) {
                mx=a[j][i];
                mxi=j;
            }
        masswap(a,b,n,i,mxi);
        for (k=i+1; k<n; k++) {
            for (j =i+1; j<n; j++) {
                if (a[k][i]>0.0 || a[k][i]<0.0) a[k][j]=a[k][j]/a[k][i]-a[i][j]/a[i][i];
            }
            if (a[k][i]>0.0 || a[k][i]<0.0) b[k]=b[k]/a[k][i]-b[i]/a[i][i];
            a[k][i]=0.0;
        }
        b[i]=b[i]/a[i][i];
        for (j=i+1; j<n; j++) a[i][j]/=a[i][i];
        a[i][i] = 1.0;
    }
    for (i=n-1; i>=0; i--) {
        for (j=n-1; j>i; j--) {
            b[i]-=a[i][j]*x[j];
        }
        x[i]=b[i];
    }
}

void matrix_multiplication(double **a, double *x, double *b, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        b[i] = 0;
        for (j = 0; j < n; j++) {
            b[i] += (a[i][j] * x[j]);
        }
    }
}

void linear_sys(double **a, double *x, double *b, int n) {
    double **a_cpy, *b_cpy, *x_nq, *b_new, *r;
    int i, j;
    
    //-----------------------------------------------------------------------
    a_cpy = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) a_cpy[i] = (double *) malloc(n * sizeof(double));
    b_cpy = (double *) malloc(n * sizeof(double));
    b_new = (double *) malloc(n * sizeof(double));
    x_nq = (double *) malloc(n * sizeof(double));
    r = (double *) malloc(n * sizeof(double));
    //------------------------------------------------------------------------
    
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a_cpy[i][j] = a[i][j];
    for (i = 0; i < n; i++) b_cpy[i] = b[i];
    gauss(a_cpy, x_nq, b_cpy, n);
    matrix_multiplication(a, x_nq, b_new, n);
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a_cpy[i][j] = a[i][j];
    for (i = 0; i < n; i++) b_cpy[i] = b[i] - b_new[i];
    gauss(a_cpy, r, b_cpy, n);
    for (i = 0; i < n; i++) x[i] = x_nq[i] - r[i];
}


void f(double *beta, double *res) {
	//double y[4];
	
	/* 6 5 4 3
	 * 4 6 5 4
	 * 2 4 6 5
	 * 0 2 4 6
	*/
	
	res[0] = 6*beta[0] + 5*beta[1] + 4*beta[2] + 3*beta[3] - 18;
	res[1] = 4*beta[0] + 6*beta[1] + 5*beta[2] + 4*beta[3] - 19;
	res[2] = 2*beta[0] + 4*beta[1] + 6*beta[2] + 5*beta[3] - 17;
	res[3] = 0*beta[0] + 2*beta[1] + 4*beta[2] + 6*beta[3] - 12;
return;
}

void gradf(double *beta, double *dbeta, double *res) {
	double beta_p[4], beta_n[4], res_p[4], res_n[4], **pro;
	int i, j;
	
	pro=(double**)malloc(4*sizeof(double*));
    for (i=0; i<4; i++) pro[i]=(double*)malloc(4*sizeof(double));
	
	for(j=0; j<4; j++) {
		for(i=0; i<4; i++) {
			beta_p[i]=beta[i];
			beta_n[i]=beta[i];
		} //сомнительное вложение
		beta_p[j]+=1e-5;
		beta_n[j]-=1e-5;
		//cout <<fixed<<setprecision(16)<< endl << "-----------------------------" << endl <<"beta_p = " << beta_p[j] << endl;
		//cout <<fixed<<setprecision(16)<< endl << "beta_n = " << beta_n[j] << endl << "-----------------------------" << endl;
		f(beta_p, res_p);
		f(beta_n, res_n);
		for(i=0; i<4; i++) {
			pro[i][j]=(res_p[i]-res_n[i])/2.e-5; //сомнительное отсутствие j
			cout <<fixed<<setprecision(16)<< pro[i][j] << " ";
		}
		cout << endl;
	}
	matrix_print(pro, 4);
	linear_sys(pro, res, dbeta, 4);
	cout << "after linear_sys";
	matrix_print(pro, 4);
return;
}

int newton(double *beta) {
	/*
	 * return >=0; - всё посчитано, вернули число итераций
	 * return -1; - одна из итераций не работает
	 * return -2; - итераций не хватило
	 * */
	
	double res[4], dbeta[4], gamma=1, beta_w[4], res_w[4];
	double eps = 2.e-15;
	int i, k, flag;
	
	f(beta, res);
	
	cout << "\n**********\n";
	printf("res[0] = %.16f\nres[1] = %.16f\nres[2] = %.16f\nres[3] = %.16f\n----------------\n", res[0], res[1], res[2], res[3]);
	//return 30;
	
	for(k=0; k<15; k++) {
		if(fabs(res[0])<eps && fabs(res[1])<eps && fabs(res[2])<eps) return k;
		gradf(beta, dbeta, res);
		gamma*=2;
		if(gamma>1) gamma=1.0;
		
		flag=-1;
		while(gamma>eps) {
			for(i=0; i<4; i++) beta_w[i]=beta[i]+gamma*dbeta[i];
			f(beta_w, res_w);
			if(norm(res_w, 4)<norm(res, 4)) {
				flag=0;
				break;
			}
			gamma/=2;
		}
		if(flag==-1) {
			printf("\nITERATION %d NOT WORK!\n", k);
			return -1;
		}
		for(i=0; i<4; i++) {
			beta[i]=beta_w[i];
			res[i]=res_w[i];
		}
		return -2;
	}
	return 20;
}

int main(void) {
	int newton_chek;
	double beta[4] = {-15, 17, 39, 18};
	newton_chek = newton(beta);
	cout << endl << newton_chek;
return 0;
}
