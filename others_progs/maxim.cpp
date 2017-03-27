/**
 * Численное решение вычислительного практикума #3
 * Задача Лагранжа, под №3 из задачника И.С. Григорьева 2005г.
 *
 * Эмин Яхьяев
 * МГУ 2016
 */

#include <iostream>
#include <cmath>
#include <algorithm>

#define EPS 1e-10

using namespace std;

double alpha = 1;

// максимальное собственное значение симметрической матрицы производных
double lambda(double *y) {
    double a = -exp(-alpha * y[0]) * (1 - alpha * y[0]);
    double b = alpha * y[3] * exp(-alpha * y[0]) * (-2 + alpha * y[0]);
    return (sqrt(pow(a + 1, 2) + b * b / 2 + 1 / 2 + fabs(b - 1) * sqrt(pow(a + 1, 2) + pow(b + 1, 2) / 4))) / 2;
}

// уравнения краевой задачи
void fcn(double t, double *y, double *f) {
    f[0] = y[1];
    f[1] = -y[0] * exp(-alpha * y[0]) + y[3];
    f[2] = y[3] * exp(-alpha * y[0]) * (1 - alpha * y[0]);
    f[3] = -y[2];
}

double dsign(double a, double b) {
    if (b < 0.0) {
        a = fabs(a) * (-1.0);
        return a;
    }
    if (b > 0.0) return fabs(a);
    return 0.0;
}

double ddopri5(int n, void f(double, double *, double *), double t, double *y,
               double tend, double hmax, double h, int print = false) {

    double k1[51], k2[51], k3[51], k4[51], k5[51], y1[51];
    bool reject;
    double xph, err, denom, fac, hnew, posneg;
    int nmax = 30000, i;
    double uround = 2.2205e-16;
    posneg = dsign(1.e0, tend - t);
    hmax = fabs(hmax);
    h = min(max(1.e-4, fabs(h)), hmax);
    h = dsign(h, posneg);
    reject = false;
    int naccpt = 0;
    int nrejct = 0;
    int nfcn = 1;
    int nstep = 0;

    double g_error = 0;

    if (print) printf("\n");

    f(t, y, k1);
    while (1) {
        //if (print) printf("%.4lf %.4lf %.4lf %.2lf\n",t,y[0],y[1],u);
        // if (print) printf("%.4lf %.4lf\n",t,y[0]);
        if (print) printf("%.4lf \t\t %.4lf %.4lf %.4lf %.4lf\n", t, y[0], y[1], y[2], y[3]);
        //if (print) printf("%.4lf %.2lf\n",t,u);
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
            denom = max(1.e-5, max(fabs(y1[i]), max(fabs(y[i]), 2.e0 * uround / EPS)));
            err += pow(k4[i] / denom, 2);
        }
        err = sqrt(err / double(n));
        fac = max(.1e0, min(5.e0, pow(err / EPS, 0.2e0) / .9e0));
        hnew = h / fac;
        if (err <= EPS) {
            if (y[1] >= EPS && y1[1] <= -EPS) {
                hnew /= 2.0;
                naccpt++;
                for (i = 0; i < n; i++) {
                    k1[i] = k3[i];
                    y[i] = y1[i];
                }
                t = xph;
                reject = !reject;
                if (naccpt >= 1) nrejct++;
            }
            else if (y[1] <= -EPS && y1[1] >= EPS) {
                hnew /= 2.0;
                naccpt++;
                for (i = 0; i < n; i++) {
                    k1[i] = k3[i];
                    y[i] = y1[i];
                }
                t = xph;
                reject = !reject;
                if (naccpt >= 1) nrejct++;
            }
            else {
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
            }

            double lamb = lambda(y);
            g_error = err + g_error * pow(M_E, h * lamb);
        }
        h = hnew;
        //if (y[1]>0.0) u=2.0;
        //else u=-2.0;
    }

    return g_error;
}

void func(double *f, double *x, int n) {
    double y[4] = {x[0], 0, 0, x[1]};
    ddopri5(4, fcn, 0.0, y, M_PI_2, M_PI_2, 0.01, false);
    f[0] = y[0] - 0.8;
    f[1] = y[3];
}

double norm(double *f, int n) {
    double norm = 0.0;
    for (int i = 0; i < n; i++) norm += (f[i] * f[i]);
    return sqrt(norm);
}

void masswap(double **a, double *b, int n, int i, int j) {
    for (int k = 0; k < n; k++) swap(a[i][k], a[j][k]);
    swap(b[i], b[j]);
}

void gauss(double **a, double *x, double *b, int n) {
    int i, j, k, mxi;
    double mx;
    for (i = 0; i < n; i++) {
        mx = a[i][i];
        mxi = i;
        for (j = i + 1; j < n; j++)
            if (a[j][i] > mx) {
                mx = a[j][i];
                mxi = j;
            }
        masswap(a, b, n, i, mxi);
        for (k = i + 1; k < n; k++) {
            for (j = i + 1; j < n; j++) {
                if (a[k][i] > 0.0 || a[k][i] < 0.0) a[k][j] = a[k][j] / a[k][i] - a[i][j] / a[i][i];
            }
            if (a[k][i] > 0.0 || a[k][i] < 0.0)b[k] = b[k] / a[k][i] - b[i] / a[i][i];
            a[k][i] = 0.0;
        }
        b[i] = b[i] / a[i][i];
        for (j = i + 1; j < n; j++) a[i][j] /= a[i][i];
        a[i][i] = 1.0;
    }
    for (i = n - 1; i >= 0; i--) {
        for (j = n - 1; j > i; j--) {
            b[i] -= a[i][j] * x[j];
        }
        x[i] = b[i];
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
    a_cpy = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) a_cpy[i] = (double *) malloc(n * sizeof(double));
    b_cpy = (double *) malloc(n * sizeof(double));
    b_new = (double *) malloc(n * sizeof(double));
    x_nq = (double *) malloc(n * sizeof(double));
    r = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a_cpy[i][j] = a[i][j];
    for (i = 0; i < n; i++) b_cpy[i] = b[i];
    gauss(a_cpy, x_nq, b_cpy, n);
    matrix_multiplication(a, x_nq, b_new, n);
    for (i = 0; i < n; i++) for (j = 0; j < n; j++) a_cpy[i][j] = a[i][j];
    for (i = 0; i < n; i++) b_cpy[i] = b[i] - b_new[i];
    gauss(a_cpy, r, b_cpy, n);
    for (i = 0; i < n; i++) x[i] = x_nq[i] - r[i];
}

void div_mtx(double **df, double *x, int n, void f(double *, double *, int)) {
    int i, j;
    double *f_p, *f_m, *x_n, h = sqrt(EPS);
    f_p = (double *) malloc(n * sizeof(double));
    f_m = (double *) malloc(n * sizeof(double));
    x_n = (double *) malloc(n * sizeof(double));
    for (i = 0; i < n; i++) x_n[i] = x[i];
    for (i = 0; i < n; i++) {
        x_n[i] = x[i] + h;
        f(f_p, x_n, n);
        x_n[i] = x[i] - h;
        f(f_m, x_n, n);
        x_n[i] = x[i];
        for (j = 0; j < n; j++) df[j][i] = (f_p[j] - f_m[j]) / (2 * h);
    }
}

void newton(double *x, int n, void f(double *, double *, int)) {
    double *fn, **dfn, *x_t, *x_n, *x_n1, nrm, nrmn, sigma;
    int i, j;
    dfn = (double **) malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) dfn[i] = (double *) malloc(n * sizeof(double));
    x_t = (double *) malloc(n * sizeof(double));
    x_n = (double *) malloc(n * sizeof(double));
    x_n1 = (double *) malloc(n * sizeof(double));
    fn = (double *) malloc(n * sizeof(double));
    f(fn, x, n);
    while (norm(fn, n) >= sqrt(EPS)) {
        //cout << x[0] << " " << x[1] << " " << x[2] << " " << norm(fn,n) << endl;
        div_mtx(dfn, x, n, f);
        linear_sys(dfn, x_t, fn, n);
        f(fn, x, n);
        nrm = norm(fn, n);
        sigma = 1.0;
        for (j = 1; j <= 12; j++) {
            for (i = 0; i < n; i++) x_n1[i] = sigma * (x[i] - x_t[i]);
            f(fn, x_n1, n);
            nrmn = norm(fn, n);
            if (nrmn < nrm) {
                for (i = 0; i < n; i++) x_n[i] = x_n1[i];
                break;
            }
            sigma /= 2.0;
        }
        for (i = 0; i < n; i++) x[i] = x_n[i];
        f(fn, x, n);
    }
}

// значение функционала
double b_func(double *b, double t) {
    double y[4] = {b[0], 0, 0, b[1]};
    ddopri5(4, fcn, 0, y, t, t, 0.01, false);

    double p2 = y[3];
    return p2 * p2;
}

// численное интегрирование
double integrate(double(*f)(double *, double), double *b, double h_start, double h_end, double eps) {
    double result = 0, h, hi, hn, x1, x2, x3, ch, del, del_sr = 0, del_gar = 0, I1, I2, I3;
    int i = 0, j = 0;
    h = (h_end + h_start) / 2;

    x1 = h_start;
    while (x1 < h_end - eps) {
        x2 = x1 + h;
        if (x2 > h_end) {
            h = h_end - x1;
            x2 = h_end;
        }
        x3 = (x2 + x1) / 2;
        I1 = (h / 6) * (f(b, x1) + 4 * f(b, x3) + f(b, x2));
        I2 = (h / 12) * (f(b, x1) + 4 * f(b, (x2 + 3 * x1) / 4) + f(b, x3));
        I3 = (h / 12) * (f(b, x3) + 4 * f(b, (x1 + 3 * x2) / 4) + f(b, x2));
        del = (I1 - I2 - I3) / 15;
        ch = fabs(del);
        hi = ch / eps;
        hi = pow(hi, 0.2);

        if (hi > 10.) hi = 10.;
        if (hi < 0.1) hi = 0.1;
        hn = 0.95 * (h / hi);
        if (ch < eps) {
            x1 = x2;
            i++;
            del_sr = del_sr + del;
            del_gar = del_gar + ch;
            result = result + I2 + I3;
        }
        else j++;
        h = hn;
    }
    return result;
}

int main(int argc, const char *argv[]) {
    // начальные параметры пристрелки
//    double b[2] = {-3.2 / M_PI, 3.2 / M_PI}; // x2, p1
    double b[2] = {-3.2 / M_PI, 3.2 / M_PI}; // x1, p2

    for (alpha = 0; alpha <= 25; alpha += 5) {
        newton(b, 2, func);
        double y[4] = {b[0], 0, -0, b[1]};
        double dz = ddopri5(4, fcn, 0, y, M_PI_2, M_PI_2, 0.01, false);

        printf("\\hline\n");
        printf("%lg &\\centering\n", alpha);
        //printf("x2(0)=%.8lf, p1(0)=%.8lf\n", b[0], b[1]);
        printf("%.8lf\\\\%.8lf &\\centering\n", b[0], b[1]);

        //printf("(x1(pi/2), x2(pi/2), p1(pi/2), p2(pi/2)) = (%.8lf, %.8lf, %.8lf, %.8lf)\n", y[0], y[1], y[2], y[3]);
        printf("%.8lf\\\\- &\\centering\n", integrate(b_func, b, 0, M_PI_2, 1.e-11));

        printf("%.3lg\\\\-\n", dz);
        printf("&\\\\\n");
//        printf("=========================================================\n\n");
    }
}

// g++ main.cpp; clear; printf '\e[3J'; ./a.out
