#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#define EPS 1e-11
#define eps2 1e-11

void fcn(double x, double *y,double *f)
{   
	f[0]=y[1];
	f[1]=y[2];
	f[2]=y[3];
	f[3]=y[0];
}
 
double dsign(double a, double b) 
{
	if (b<0) { a=fabs(a)*(-1.0); return a; }
	if (b>=0) return fabs(a);
	return 0.0;
}


void rgk(int n,void (*fcn)(double, double*,double*),double x,double *y,double xend, double eps, double hmax,double h)
{
    double k1[10],k2[10],k3[10],k4[10],k5[10],y1[10];
    bool reject;
    double xph,err,denom,fac,hnew,posneg;
    int nmax=30000,i;
    double uround=2.2205e-16;
    posneg=dsign(1.e0,xend-x);
    hmax=fabs(hmax);
    h=min(max(1.e-4,fabs(h)),hmax);
    h=dsign(h,posneg);
    eps=max(eps,7.e0*uround);
    reject=false;
    int naccpt=0;
    int nrejct=0;
    int nfcn=1;
    int nstep=0;    
    fcn(x,y,k1);    
    
    while (1)
     {
		if ( nstep > nmax || x+.1e0*h==x ) break;
		if ( (x-xend)*posneg+uround > 0.e0) break;
		if ( (x+h-xend)*posneg > 0.e0) h=xend-x;
		nstep++;
		for (i=0; i<n; i++) y1[i]=y[i]+h*.2e0*k1[i];
		fcn(x+h*.2e0,y1,k2);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((3.e0/40.e0)*k1[i]+(9.e0/40.e0)*k2[i]);
		fcn(x+h*.3e0,y1,k3);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((44.e0/45.e0)*k1[i]-(56.e0/15.e0)*k2[i]+(32.e0/9.e0)*k3[i]);
		fcn(x+h*.8e0,y1,k4);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((19372.e0/6561.e0)*k1[i]-(25360.e0/2187.e0)*k2[i]+(64448.e0/6561.e0)*k3[i]-(212.e0/729.e0)*k4[i]);
		fcn(x+h*(8.e0/9.e0),y1,k5);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((9017.e0/3168.e0)*k1[i]-(355.e0/33.e0)*k2[i]+(46732.e0/5247.e0)*k3[i]+(49.e0/176.e0)*k4[i]-(5103.e0/18656.e0)*k5[i]);
		xph=x+h;
		fcn(xph,y1,k2);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((35.e0/384.e0)*k1[i]+(500.e0/1113.e0)*k3[i]+(125.e0/192.e0)*k4[i]-(2187.e0/6784.e0)*k5[i]+(11.e0/84.e0)*k2[i]);
		for (i=0; i<n; i++) k2[i]=(71.e0/57600.e0)*k1[i]-(71.e0/16695.e0)*k3[i]+(71.e0/1920.e0)*k4[i]-(17253.e0/339200.e0)*k5[i]+(22.e0/525.e0)*k2[i];
		fcn(xph,y1,k3);
		for (i=0; i<n; i++) k4[i]=(k2[i]-(1.e0/40.e0)*k3[i])*h;
		nfcn+=6;
		err=0;

   for (i=0; i<n; i++) {
			denom=max(1.e-5,max(fabs(y1[i]),max(fabs(y[i]),2.e0*uround/eps)));
			err+=pow(k4[i]/denom,2);
		}
  
      err=sqrt(err/double(n));
		fac=max( .1e0, min( 5.e0, pow( err/eps,0.2e0 )/.9e0) );
		hnew=h/fac;
		
   if(err<=eps) {
			naccpt++;
			for (i=0; i<n; i++) {
				k1[i]=k3[i];
				y[i]=y1[i];
			}
    x=xph;
			if(fabs(hnew)>hmax) hnew=posneg*hmax;
			if(reject) hnew=posneg*min(fabs(hnew),fabs(h)),reject=false;
			else reject=true;
			if(naccpt >= 1) nrejct++;
	}
		h=hnew;
}
   printf("%lf: \n\t%.16f\n\t%.16f\n\t%.16f\n\t%.16f\n", x, y[0], y[1], y[2], y[3]);
}

void l(double *beta, double *res)
{
//double eps=1.e-11;
double y[4];
   y[0]=0;
   y[1]=0;
   y[2]=beta[0];
   y[3]=beta[1];
   rgk(4,fcn,0.0e0,y,1.0e0,EPS,1.0e0,0.5e0);
   res[0]=y[0]-5.0;
   res[1]=y[1]-2.0;
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
    
    a_cpy = (double **) malloc(n * sizeof(double *)); //копия матрицы
    for (i = 0; i < n; i++) a_cpy[i] = (double *) malloc(n * sizeof(double));
    b_cpy = (double *) malloc(n * sizeof(double)); //копия вектора
    b_new = (double *) malloc(n * sizeof(double)); //новый вектор
    x_nq = (double *) malloc(n * sizeof(double));
    r = (double *) malloc(n * sizeof(double));
    
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a_cpy[i][j] = a[i][j]; //копирование матрицы
    for (i = 0; i < n; i++) b_cpy[i] = b[i]; //копирование вектора
        
    gauss(a_cpy, x_nq, b_cpy, n);
        
    matrix_multiplication(a, x_nq, b_new, n);
    
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a_cpy[i][j] = a[i][j];
    for (i = 0; i < n; i++) b_cpy[i] = b[i] - b_new[i];
    
    gauss(a_cpy, r, b_cpy, n);
    for (i = 0; i < n; i++) x[i] = x_nq[i] - r[i];

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

void kramer(double **pro,double *dbeta, double *res)
{
	double delta_1,delta_2,delta;
	delta_1=res[1]*pro[0][1]-res[0]*pro[1][1];
	delta_2=res[0]*pro[1][0]-res[1]*pro[0][0];
	delta=pro[0][0]*pro[1][1]-pro[0][1]*pro[1][0];
	dbeta[0]=delta_1/delta;
	dbeta[1]=delta_2/delta;
	}



void gradf(double *beta,double *dbeta, double *res)
{
 int j,i,N=2;
 double beta_p[N],beta_m[N],res_p[N],res_m[N];
 double h=1.e-6;
 double res_new[N];
 double **pro;
  pro=(double**)malloc(N*sizeof(double*));
  for(i=0;i<N;i++) pro[i]=(double*)malloc(N*sizeof(double));
 
for(j=0;j<N;j++)
 {
 
 for(i=0;i<N;i++) {
		beta_p[i]=beta[i];
		beta_m[i]=beta[i];
	 }
	 beta_p[j]+=h;
	 beta_m[j]-=h;
	 
	 l(beta_p,res_p);
	 l(beta_m,res_m);
	 for(i=0;i<N;i++)pro[i][j]=(res_p[i]-res_m[i])/(2*h);
}
  //kramer(pro,dbeta,res);
  for(i=0; i<N; i++) res_new[i]=-res[i];
  linear_sys(pro,dbeta,res_new,N);
  
  /*
  for(i=0;i<2;i++){for(j=0;j<2;j++)printf("pro[%d][%d]=%lf \n",i,j,pro[i][j]);}
     printf("dbeta[0] = %.16f ,  dbeta[1]=%.16f \n", dbeta[0], dbeta[1]);
	*/
}


double newton(double *beta) {
	int i, k, N=2;
	double dbeta[N],beta_w[N],gamma=1,res_w[N], res[N];
	bool flag;

l(beta,res);

for(k=0;k<15;k++)
{
  if(fabs(res[0])<eps2 && fabs(res[1])<eps2) {printf("kolichestvo itecaciy %d \n",k);return k;}//повезло нам
  
  gradf(beta,dbeta,res);	  

  gamma=1.0;
  flag=false;
  while(gamma>eps2) {
	for(i=0;i<N;i++) {
		beta_w[i]=beta[i]+gamma*dbeta[i]; 
		printf("beta_w[%d]=%lf \n", i,beta_w[i]);
	}
	
	l(beta_w,res_w);
		
	if(norm(res_w,N)<norm(res,N)) {flag=true; break;}
	gamma/=2;

	if(flag==false){printf("slomalsya na %d iteracii\n",k);return -1;} //iteraciya ne rabotaet
  }
  	
  for(i=0;i<N;i++) {beta[i]=beta_w[i],res[i]=res_w[i];}
}
	
return -2; //iteracii zakonchilis

}

int main() {
 double beta[2];
 
	printf("vvesti beta_1= \n");	
	scanf("%lf",&beta[0]);
	printf("vvesti beta_2= \n");
	scanf("%lf",&beta[1]);

	newton(beta);
return 0;	
}
