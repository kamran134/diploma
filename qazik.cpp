#include <cstdio>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

#define epsilon 1e-10
#define epsilon_newton 2e-15

using namespace std;

double alpha;

double lambda(double *y, double t) {
    double a = -1 / (1 + alpha * t * t) + 1;
    return (sqrt(pow(a, 2) + 1 / 2 + sqrt(pow(a, 2) + 1 / 4))) / 2;
}

double dsign(double a, double b) {
    return ((b < 0.0) ? fabs(a) * (-1.0) :
             ((b > 0.0) ? fabs(a) :
               0.0));
}

void fcn(double t, double *y, double *f) {
    f[0] = y[1];			//x2
	f[1] = y[2];			//x3
	f[2] = y[3];			//x4
	f[3] = y[0];			//x1
    /*
    f[0] = y[1];
    f[1] = -y[0] / (1 + alpha * t * t) + y[3];
    f[2] = y[3] / (1 + alpha * t * t);
    f[3] = -y[2];
    */
}

double ddopri5(int n, void f(double, double *, double *), double t, double *y,
               double tend, double hmax, double h) {

    bool reject = false;
    double k1[51], k2[51], k3[51], k4[51], k5[51], y1[51];
    double lamb, g_error = 0, xph, err, denom, fac, hnew, posneg, uround = 2.2205e-16;
    int i, nfcn = 1, nstep = 0, nrejct = 0, naccpt = 0, nmax = 30000;
    posneg = dsign(1.e0, tend - t);
    hmax = fabs(hmax);
    h = min(max(1.e-4, fabs(h)), hmax);
    h = dsign(h, posneg);

    f(t, y, k1);

    while (1) {
        if (nstep > nmax) break;
        if ((t - tend) * posneg + uround > 0.e0) break;
        if ((t + h - tend) * posneg > 0.e0) h = tend - t;
        nstep++;
        for (i = 0; i < n; i++) y1[i] = y[i] + h * .2e0 * k1[i];
        f(t + h * .2e0, y1, k2);
        for (i = 0; i < n; i++) y1[i] = y[i] + h * ((3.e0 / 40.e0) * k1[i] + (9.e0 / 40.e0) * k2[i]);
        f(t + h * .3e0, y1, k3);
        for (i = 0; i < n; i++)
            y1[i] = y[i] + h * ((44.e0 / 45.e0) * k1[i] - (56.e0 / 15.e0) * k2[i] + (32.e0 / 9.e0) * k3[i]);
        f(t + h * .8e0, y1, k4);
        for (i = 0; i < n; i++)
            y1[i] = y[i] + h * ((19372.e0 / 6561.e0) * k1[i] - (25360.e0 / 2187.e0) * k2[i] +
                                (64448.e0 / 6561.e0) * k3[i] - (212.e0 / 729.e0) * k4[i]);
        f(t + h * (8.e0 / 9.e0), y1, k5);
        for (i = 0; i < n; i++)
            y1[i] = y[i] + h * ((9017.e0 / 3168.e0) * k1[i] - (355.e0 / 33.e0) * k2[i] + (46732.e0 / 5247.e0) * k3[i] +
                                (49.e0 / 176.e0) * k4[i] - (5103.e0 / 18656.e0) * k5[i]);
        xph = t + h;
        f(xph, y1, k2);
        for (i = 0; i < n; i++)
            y1[i] = y[i] + h * ((35.e0 / 384.e0) * k1[i] + (500.e0 / 1113.e0) * k3[i] + (125.e0 / 192.e0) * k4[i] -
                                (2187.e0 / 6784.e0) * k5[i] + (11.e0 / 84.e0) * k2[i]);
        for (i = 0; i < n; i++)
            k2[i] = (71.e0 / 57600.e0) * k1[i] - (71.e0 / 16695.e0) * k3[i] + (71.e0 / 1920.e0) * k4[i] -
                    (17253.e0 / 339200.e0) * k5[i] + (22.e0 / 525.e0) * k2[i];
        f(xph, y1, k3);
        for (i = 0; i < n; i++) k4[i] = (k2[i] - (1.e0 / 40.e0) * k3[i]) * h;
        nfcn += 6;
        err = 0;
        for (i = 0; i < n; i++) {
            denom = max(1.e-5, max(fabs(y1[i]), max(fabs(y[i]), 2.e0 * uround / epsilon)));
            err += pow(k4[i] / denom, 2);
        }
        err = sqrt(err / double(n));
        fac = max(.1e0, min(5.e0, pow(err / epsilon, 0.2e0) / .9e0));
        hnew = h / fac;
        if (err <= epsilon) {
            naccpt++;
            for (i = 0; i < n; i++) {
                k1[i] = k3[i];
                y[i] = y1[i];
            }
            t = xph;
            if (fabs(hnew) > hmax) hnew = posneg * hmax;
            if (reject) hnew = posneg * min(fabs(hnew), fabs(h)), reject = false;
            else reject = true;
            if (naccpt >= 1) nrejct++;

            lamb = lambda(y, t);
            g_error = err + g_error * pow(M_E, h * lamb);
        }
        h = hnew;
        
        cout << t << endl;
		
		cout <<fixed<<setprecision(16)<< "\t\tx(" << t << ") = " << y[0] << endl;
		cout <<fixed<<setprecision(16)<< "\t\ty(" << t << ") = " << y[1] << endl;
		
		cout <<fixed<<setprecision(16)<< "\t\tv(" << t << ") = " << y[2] << endl;
		cout <<fixed<<setprecision(16)<< "\t\tpsi_x(" << t << ") = " << y[3] << endl;
        
    }

    return g_error;
}


void func(double *res, double *beta) {
    double y[4] = {0, 0, beta[0], beta[1]};  //если {0, beta[0], beta[1], 0} тогда сходится за 3 итерации, но Илья хочет так: x_1=0, x_2=0, p_1=1, p_2=1
    ddopri5(4, fcn, 0.0, y, 1.0e0, 1.0e0, 0.5e0);
    res[0] = y[0]-5;
    res[1] = y[1]-2;
}

double norm(double *vector) {
	int i, N=2;
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
//return max_W;
}

double norm2(double *f) {
    double norm = 0.0, n = 2;
    for (int i = 0; i < n; i++) norm += (f[i] * f[i]);
    return sqrt(norm);
}

void gradf(double *beta, double **proizvodnaya) {
    int i, j, n = 2;
    double beta_p[n], beta_m[n], res_p[n], res_m[n], h = 1.e-6;

    for(j = 0; j < n; j++){
        for(i = 0; i < n; i++){
            beta_p[i] = beta[i];
            beta_m[i] = beta[i];
        }

        beta_p[j] += h;
        beta_m[j] -= h;

        func(res_p, beta_p);
        func(res_m, beta_m);

        for (i = 0; i < n; i++) { proizvodnaya[i][j] = (res_p[i] - res_m[i]) / (2 * h); }
    }
}

void kramer(double **proizvodnaya, double *dbeta, double *res) {
    double delta_1, delta_2, delta;

    delta_1 = res[1] * proizvodnaya[0][1] - res[0] * proizvodnaya[1][1];
    delta_2 = res[0] * proizvodnaya[1][0] - res[1] * proizvodnaya[0][0];
    delta = proizvodnaya[0][0] * proizvodnaya[1][1] - proizvodnaya[0][1] * proizvodnaya[1][0];
    dbeta[0] = delta_1 / delta;
    dbeta[1] = delta_2 / delta;
}

#include "prak3_func/gauss_m.cpp"

int newton_method(double *beta) {
    double res[2], res_w[2], beta_w[2], dbeta[2], gamma, **proizvodnaya, res_m[2];
    int i, k, flag, n = 2;

    proizvodnaya = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) proizvodnaya[i] = (double *) malloc(n * sizeof(double));

    func(res, beta);

    for(k = 0; k < 15; k++){
        if(fabs(res[0]) < epsilon_newton && fabs(res[1]) < epsilon_newton) {cout << endl << "ended by " << k << " iteration" << endl; return k; }

        gradf(beta, proizvodnaya);
        //kramer(proizvodnaya, dbeta, res);
        for(i=0; i<n; i++) res_m[i] = -res[i];
        linear_sys(proizvodnaya, dbeta, res_m, n);
        gamma = 1.0; flag = -1;

        while(gamma > epsilon_newton){
            for(i = 0; i < n; i++){ beta_w[i] = beta[i] + gamma * dbeta[i]; }
                func(res_w, beta_w);
                if(norm(res_w) < norm(res)) {flag = 0; break;}
            gamma /= 2.0;
        }

        if(flag == -1){cout << endl << "broked on " << k << "'s iteration" << endl; return -1; }

        for(i = 0; i < n; i++){ beta[i] = beta_w[i]; res[i] = res_w[i]; }
    }

    return -2;
}


int main(int argc, const char *argv[]) {
	//double y[4];
	double beta[2] = {0,0};
    
    //ddopri5(4, fcn, 0, y, 1, 1.0e0, 0.5e0);
    newton_method(beta);
    
    return 0;
    /*
    double beta[2] = {1, 2};
    int i, j, newton_result, alpha_mas_count;
    double alpha_mas[] = {0.0, 0.01, 0.1, 1, 1.02, 2.21, 2.5, 5, 7.3, 8.5, 10, 10.01, 11};
    alpha_mas_count = sizeof(alpha_mas)/sizeof(*alpha_mas);

    for (j = 0; j < alpha_mas_count; j++) {
        alpha = alpha_mas[j];

        beta[0] = 1; beta[1] = 2;

        newton_result = newton_method(beta);

        printf("alpha = %.8lf\n", alpha);

        for(i = 0; i < 2; i++){
            printf("beta[%d] = %.10lf\n", i, beta[i]);
        }
        printf("newton_result = %d\n", newton_result);
        //double y[4] = {0, beta[0], beta[1], 0};
        //double dz = ddopri5(4, fcn, 0, y, M_PI_2, M_PI_2, 0.01);

        //printf("x2(0)=%.16lf, p1(0)=%.16lf\n", beta[0], beta[1]);
	*/
        /*//printf("%.8lf\\\\%.8lf &\\centering\n", b[0], b[1]);
        printf("(x1(pi/2), x2(pi/2), p1(pi/2), p2(pi/2)) = (%.8lf, %.8lf, %.8lf, %.8lf)\n", y[0], y[1], y[2], y[3]);*/

    //    cout << "=========================================================" << endl << endl;
    //}
    
    
}
