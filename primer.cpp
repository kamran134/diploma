#include <stdio.h>
#include <iostream>
#include <math.h>
using namespace std;

#define EPS_DDOPRI 1e-11
#define ens 1e-10 //epsilon newton stop
#define eps_gamma 2e-15

void fcn(double x, double y[],double f[])
{   
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
	f[5]=y[5]-(y[2]*y[3]*y[3]+y[2]*y[4]*y[4]+y[4]*y[5])/RHO;	//psi_v^{\cdot}
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

void l(double *betta, double *res)
{
double y[4];
	y[0]=0; //x_0
	y[1]=10; //y_0
	y[2]=0; //v_0
	y[3]=betta[0];
	y[4]=betta[1];              
	y[5]=betta[2];
   rgk(6,fcn,0,y,betta[3],EPS_DDOPRI,1,0.5);
   res[0]=y[0]-5.0;
   res[1]=y[1]-2.0;
	// printf("y[0]=%lf y[1]=%lf y[2]=%lf y[3]=%lf \n",y[0],y[1],y[2],y[3]);
	// printf("res[0]=%lf res[1]=%lf \n",res[0],res[1]);
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



void kramer(double **pro,double *dbetta, double *res)
{
	double delta_1,delta_2,delta;
	delta_1=res[1]*pro[0][1]-res[0]*pro[1][1];
	delta_2=res[0]*pro[1][0]-res[1]*pro[0][0];
	delta=pro[0][0]*pro[1][1]-pro[0][1]*pro[1][0];
	dbetta[0]=delta_1/delta;
	dbetta[1]=delta_2/delta;
	}



void gradf(double *betta,double *dbetta, double *res, double **pro, double *y)
{
 double betta_p[2],betta_m[2],res_p[2],res_m[2];
 int j,i;
 //double **pro;
 double res_new[2];
  //pro=(double**)malloc(2*sizeof(double*));
  //for(i=0;i<2;i++) pro[i]=(double*)malloc(2*sizeof(double));
 for(j=0;j<2;j++)
 {
 
 for(i=0;i<2;i++) 
	 {
		betta_p[i]=betta[i];
		betta_m[i]=betta[i];
	 }
	 betta_p[j]+=1.e-5;
	 betta_m[j]-=1.e-5;
	 l(betta_p,res_p);
	 l(betta_m,res_m);
	 for(i=0;i<2;i++)pro[i][j]=(res_p[i]-res_m[i])/2.e-5;
	 }
  //kramer(pro,dbetta,res);
  for(i=0; i<2; i++) res_new[i]=-res[i];
  linear_sys(pro,dbetta,res_new,2);
  for(i=0;i<2;i++){for(j=0;j<2;j++)printf("pro[%d][%d]=%lf \n",i,j,pro[i][j]);}
     printf("dbetta[0] = %.16f ,  dbetta[1]=%.16f \n", dbetta[0], dbetta[1]);
}


double newton(double *betta, double *y)
{
double dbetta[2],betta_w[2],gamma=1,res_w[2], res[2];
int i,k;
bool flag;
double **pro;

pro=(double**)malloc(2*sizeof(double*));
for(i=0;i<2;i++) pro[i]=(double*)malloc(2*sizeof(double));
l(betta,res);

for(k=0;k<15;k++)
{
  if(fabs(res[0])<ens && fabs(res[1])<ens) {printf("kolichestvo itecaciy %d \n",k);return k;}//повезло нам
  gradf(betta,dbetta,res,pro,y);	  

  flag=false;
  gamma=1.0;
  while(gamma>eps_gamma) {
		for(i=0;i<2;i++)
		{
			betta_w[i]=betta[i]+gamma*dbetta[i]; 
			 printf("betta_w[%d]=%lf \n", i,betta_w[i]);
		}
			l(betta_w,res_w);		
		
  if(sqrt(res_w[0]*res_w[0]+res_w[1]*res_w[1])<sqrt(res[0]*res[0]+res[1]*res[1])){flag=true; break;}
  gamma/=2;

		if(flag==false) {
			printf("slomalsya na %d iteracii\n",k);
			return -1;
		} //iteraciya ne rabotaet
	}
		
	for(i=0;i<2;i++){betta[i]=betta_w[i],res[i]=res_w[i];}
	}
	
return -2; //iteracii zakonchilis

}

int main()
{
 double betta[4]={1,1,1,2};
 //int i,j;
 //double dbetta[2];
 double y[6];
 	 // printf("vvesti betta_1= \n");	
	 // scanf("%lf",&betta[0]);
	 // printf("vvesti betta_2= \n");
	 // scanf("%lf",&betta[1]);
   y[0]=0; //x_0
	y[1]=10; //y_0
	y[2]=0; //v_0
	y[3]=betta[0];
	y[4]=betta[1];              
	y[5]=betta[2];
  rgk(6,fcn,0,y,betta[3],EPS_DDOPRI,1,0.5);
	 //  l(betta,res,y);
   //newton(betta,y);	  
 //  gradf(betta,dbetta,res,pro,y);
   //printf("dbetta_1 =%ld dbetta_2 =%ld \n",dbetta[0],dbetta[1]);
  // printf("res[0]=%lf res[1]=%lf \n",res[0],res[1]);
 // dx1=0*sin(1)-0.5*cos(1)+0.25*exp(1)+0.25*exp(-1)-y[0];
  // printf("y[0]=%.16f\n y[1]=%.16f \n y[2]=%.16f \n y[3]=%.16f \n",y[0],y[1],y[2],y[3]);
//   printf("dx1=%.16f \n",dx1);

 //printf("otvet integrala = %lf \n",y[3]*y[3]-y[1]*y[1]-y[0]*y[0]);
return 0;	
}
