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
