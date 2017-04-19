#include <iostream>
#include <cmath>
#define gamma 10
using namespace std;

void fcn(double x, double *y,double *f) {
    f[1]=-1+(1+gamma*sin(x))*cos(y[1])*cos(y[1]);
}

double dsign(double a, double b) {
	if (b<0) { a=fabs(a)*(-1.0); return a; }
	if (b>0) return fabs(a);
	return 0.0;
}

void ddopri5(int n,void (*fcn)(double, double*,double*),double x,double *y,double xend,double hmax,double h) {
	double eps=1e-10;
    double k1[51],k2[51],k3[51],k4[51],k5[51],y1[51];
    bool reject;
    double xph,err,denom,fac,hnew,posneg;
    int nmax=30000,i;
    double uround=2.2205e-16;
    posneg=dsign(1.e0,xend-x);
    hmax=fabs(hmax);
    h=min(max(1.e-4,fabs(h)),hmax);
    h=dsign(h,posneg);
    reject=false;
    int naccpt=0;
    int nrejct=0;
    int nfcn=1;
    int nstep=0;
    fcn(x,y,k1);
    while (1) {
		if ( nstep > nmax ) break;
		if ( (x-xend)*posneg+uround > 0.e0) break;
		if ( (x+h-xend)*posneg > 0.e0) h=xend-x;
		nstep++;
		for (i=1; i<=n; i++) y1[i]=y[i]+h*.2e0*k1[i];
		fcn(x+h*.2e0,y1,k2);
		for (i=1; i<=n; i++) y1[i]=y[i]+h*((3.e0/40.e0)*k1[i]+(9.e0/40.e0)*k2[i]);
		fcn(x+h*.3e0,y1,k3);
		for (i=1; i<=n; i++) y1[i]=y[i]+h*((44.e0/45.e0)*k1[i]-(56.e0/15.e0)*k2[i]+(32.e0/9.e0)*k3[i]);
		fcn(x+h*.8e0,y1,k4);
		for (i=1; i<=n; i++) y1[i]=y[i]+h*((19372.e0/6561.e0)*k1[i]-(25360.e0/2187.e0)*k2[i]+(64448.e0/6561.e0)*k3[i]-(212.e0/729.e0)*k4[i]);
		fcn(x+h*(8.e0/9.e0),y1,k5);
		for (i=1; i<=n; i++) y1[i]=y[i]+h*((9017.e0/3168.e0)*k1[i]-(355.e0/33.e0)*k2[i]+(46732.e0/5247.e0)*k3[i]+(49.e0/176.e0)*k4[i]-(5103.e0/18656.e0)*k5[i]);
		xph=x+h;
		fcn(xph,y1,k2);
		for (i=1; i<=n; i++) y1[i]=y[i]+h*((35.e0/384.e0)*k1[i]+(500.e0/1113.e0)*k3[i]+(125.e0/192.e0)*k4[i]-(2187.e0/6784.e0)*k5[i]+(11.e0/84.e0)*k2[i]);
		for (i=1; i<=n; i++) k2[i]=(71.e0/57600.e0)*k1[i]-(71.e0/16695.e0)*k3[i]+(71.e0/1920.e0)*k4[i]-(17253.e0/339200.e0)*k5[i]+(22.e0/525.e0)*k2[i];
		fcn(xph,y1,k3);
		for (i=1; i<=n; i++) k4[i]=(k2[i]-(1.e0/40.e0)*k3[i])*h;
		nfcn+=6;
		err=0;
		for (i=1; i<=n; i++) {
			denom=max(1.e-5,max( fabs(y1[i]) , max(fabs(y[i]),2.e0*uround/eps)));
			err+=pow(k4[i]/denom,2);
		}
		err=sqrt(err/double(n));
		fac=max( .1e0, min( 5.e0, pow( err/eps,0.2e0 )/.9e0) );
		hnew=h/fac;
		if(err <= eps) {
			naccpt++;
			for (i=1; i<=n; i++) {
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
		//cout << x << " " << y[1] << " " << y[2] << " " << y[3] << " " << y[4] << endl;
    }
	//x1=x;
}

int main()
{
	double i;
	double y[3];
	for (i=-M_PI/2; i<=M_PI/2; i+=0.001)
	{
	y[1]=i;
	ddopri5(1,fcn,0,y,2*M_PI,2*M_PI,0.1);
	cout << i << " " << y[1] << endl;
	}
	return 0;
}
