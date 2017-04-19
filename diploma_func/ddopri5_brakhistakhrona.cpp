double dsign(double a, double b) {
	if (b<0) { a=fabs(a)*(-1.0); return a; }
	if (b>0) return fabs(a);
	return 0.0;
}

double max_root(double t, double lambda) {
	return max(-0.5+0.5*lambda*t, 0.5-0.5*lambda*t);
}

double ddopri5(int n,void (*fcn)(double, double*,double*),double t,double *y,double tend, double eps, double hmax,double h) {	
//double ddopri5(int n,double t,double *y,double tend, double eps, double hmax,double h) {    
    double k1[11],k2[11],k3[11],k4[11],k5[11],y1[11];
    bool reject;
    double tph,err,denom,fac,hnew,posneg;
    int nmax=30000,i;
    double uround=2.2205e-16;
    double gerror=0;
    
    posneg=dsign(1.e0,tend-t);
    hmax=fabs(hmax);
    h=min(max(1.e-4,fabs(h)),hmax);
    h=dsign(h,posneg);
    reject=false;
    int naccpt=0;
    int nrejct=0;
    int nfcn=1;
    int nstep=0;
    eps=max(eps,7.e0*uround);
    fcn(t,y,k1);
    while (1) {
		if ( nstep > nmax ) break;
		if ( (t-tend)*posneg+uround > 0.e0) break;
		if ( (t+h-tend)*posneg > 0.e0) h=tend-t;
		nstep++;
		for (i=0; i<n; i++) y1[i]=y[i]+h*.2e0*k1[i];
		fcn(t+h*.2e0,y1,k2);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((3.e0/40.e0)*k1[i]+(9.e0/40.e0)*k2[i]);
		fcn(t+h*.3e0,y1,k3);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((44.e0/45.e0)*k1[i]-(56.e0/15.e0)*k2[i]+(32.e0/9.e0)*k3[i]);
		fcn(t+h*.8e0,y1,k4);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((19372.e0/6561.e0)*k1[i]-(25360.e0/2187.e0)*k2[i]+(64448.e0/6561.e0)*k3[i]-(212.e0/729.e0)*k4[i]);
		fcn(t+h*(8.e0/9.e0),y1,k5);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((9017.e0/3168.e0)*k1[i]-(355.e0/33.e0)*k2[i]+(46732.e0/5247.e0)*k3[i]+(49.e0/176.e0)*k4[i]-(5103.e0/18656.e0)*k5[i]);
		tph=t+h;
		fcn(tph,y1,k2);
		for (i=0; i<n; i++) y1[i]=y[i]+h*((35.e0/384.e0)*k1[i]+(500.e0/1113.e0)*k3[i]+(125.e0/192.e0)*k4[i]-(2187.e0/6784.e0)*k5[i]+(11.e0/84.e0)*k2[i]);
		for (i=0; i<n; i++) k2[i]=(71.e0/57600.e0)*k1[i]-(71.e0/16695.e0)*k3[i]+(71.e0/1920.e0)*k4[i]-(17253.e0/339200.e0)*k5[i]+(22.e0/525.e0)*k2[i];
		fcn(tph,y1,k3);
		for (i=0; i<n; i++) k4[i]=(k2[i]-(1.e0/40.e0)*k3[i])*h;
		nfcn+=6;
		err=0;
		for (i=0; i<n; i++) {
			denom=max(1.e-5,max( fabs(y1[i]) , max(fabs(y[i]),2.e0*uround/eps)));
			err+=pow(k4[i]/denom,2);
		}
		err=sqrt(err/double(n));
		fac=max( .1e0, min( 5.e0, pow( err/eps,0.2e0 )/.9e0) );
		hnew=h/fac;
		if(err <= eps) {
			naccpt++;
			for (i=0; i<n; i++) {
				k1[i]=k3[i];
				y[i]=y1[i];
			}
			t=tph;
			if(fabs(hnew)>hmax) hnew=posneg*hmax;
			if(reject) hnew=posneg*min(fabs(hnew),fabs(h)),reject=false;
			else reject=true;
			if(naccpt >= 1) nrejct++;
			double omega = max_root(t, y[3]);
			gerror = err + gerror*pow(M_E, h*omega); //gerror это то, что в теории называется лямбдой маленькой, а err это лямбда большая (треугольная)
		}
		h=hnew;
    }
		
		cout << t << endl;
		cout << y[0] << " , " << y[1] << " , " << y[2] << endl;
		cout << y[3] << " , " << y[4] << " , " << y[5] << endl;
		cout << "-------------------------------" << endl;
    return gerror;
}
