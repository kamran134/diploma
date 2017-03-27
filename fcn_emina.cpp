void fcn(double t, double *y, double *f) {
    f[0] = y[1];
    f[1] = -y[0] * exp(-alpha * y[0]) + y[3];
    f[2] = y[3] * exp(-alpha * y[0]) * (1 - alpha * y[0]);
    f[3] = -y[2];
}

void gradf(double *x, double *dbeta, double *res, int n) {
//void gradf(double **df, double *x, int n, void f(double *, double *, int)) {
    int i, j;
    double eps = 2.e-15;
    double **pro;
    double *f_p, *f_m, *x_n, h = sqrt(eps);
    
    
    f_p = (double *) malloc(n * sizeof(double));
    f_m = (double *) malloc(n * sizeof(double));
    x_n = (double *) malloc(n * sizeof(double));
    
    
	pro=(double**)malloc(n*sizeof(double*));
    for (i=0; i<n; i++) pro[i]=(double*)malloc(n*sizeof(double));
    
    for (i = 0; i < n; i++) x_n[i] = x[i];
    for (i = 0; i < n; i++) {
        x_n[i] = x[i] + h;
        fcn(f_p, x_n, n);
        x_n[i] = x[i] - h;
        fcn(f_m, x_n, n);
        x_n[i] = x[i];
        for (j = 0; j < n; j++) pro[j][i] = (f_p[j] - f_m[j]) / (2 * h);
    }
    linear_sys(pro, dbeta, res, N);
}
