////////////////////////////////////////////////////////////////////////////////
// File: embedded_prince_dormand_v1_4_5.c                                     //
// Routines:                                                                  //
//    Embedded_Prince_Dormand_v1_4_5                                          //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  Description:                                                              //
//     This Runge-Kutta-Prince-Dormand method is an adaptive procedure for    //
//     approximating the solution of the differential equation y'(x) = f(x,y) //
//     with initial condition y(x0) = c.  This implementation evaluates       //
//     f(x,y) six times per step using embedded fourth order and fifth order  //
//     Runge-Kutta estimates to estimate the not only the solution but also   //
//     the error.                                                             //
//     The next step size is then calculated using the preassigned tolerance  //
//     and error estimate.                                                    //
//     For step i+1,                                                          //
//        y[i+1] = y[i] +  h * ( 31 / 540 * k1 + 190 / 297 * k3               //
//                           - 145 / 108 * k4 + 351 / 220 * k5 + 1/20 * k6 )  //
//     where                                                                  //
//     k1 = f( x[i],y[i] ),                                                   //
//     k2 = f( x[i]+h/5, y[i] + h*k1/5 ),                                     //
//     k3 = f( x[i]+3h/10, y[i]+(h/40)*(3 k1 + 9 k2) ),                       //
//     k4 = f( x[i]+3h/5, y[i]+(h/10)*(3 k1 - 9 k2 + 12 k3) ),                //
//     k5 = f( x[i]+2h/3, y[i]+(h/729)*(226 k1 - 675 k2 + 880 k3 + 55 k4) )   //
//     k6 = f( x[i]+h, y[i]+(h/2970)*(-1991 k1 + 7425 k2 - 2660 k3            //
//                                                  - 10010 k4 + 10206 k5) )  //
//     x[i+1] = x[i] + h.                                                     //
//                                                                            //
//     The error is estimated to be                                           //
//        err = h*( 77 k1 - 400 k3 + 1925 k4 - 1701 k5 + 99 k6 ) / 2520       //
//     The step size h is then scaled by the scale factor                     //
//         scale = 0.8 * | epsilon * y[i] / [err * (xmax - x[0])] | ^ 1/4     //
//     The scale factor is further constrained 0.125 < scale < 4.0.           //
//     The new step size is h := scale * h.                                   //
////////////////////////////////////////////////////////////////////////////////

#include <math.h>

#define max(x, y) ( (x) < (y) ? (y) : (x) )
#define min(x, y) ( (x) < (y) ? (x) : (y) )

#define ATTEMPTS 12
#define MIN_SCALE_FACTOR 0.125
#define MAX_SCALE_FACTOR 4.0

static double Runge_Kutta(double (*f)(double, double), double *y, double x, double h);

////////////////////////////////////////////////////////////////////////////////
// int Embedded_Prince_Dormand_v1_4_5( double (*f)(double, double),           //
//       double y[], double x, double h, double xmax, double *h_next,         //
//                                                        double tolerance )  //
//                                                                            //
//  Description:                                                              //
//     This function solves the differential equation y'=f(x,y) with the      //
//     initial condition y(x) = y[0].  The value at xmax is returned in y[1]. //
//     The function returns 0 if successful or -1 if it fails.                //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y0) corresponding to the //
//                initial condition y(x0) = y0.                               //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at xmax.                               //
//     double x   The initial value of x.                                     //
//     double h   Initial step size.                                          //
//     double xmax The endpoint of x.                                         //
//     double *h_next   A pointer to the estimated step size for successive   //
//                      calls to Runge_Kutta_Prince_Dormand_v1_4_5.           //
//     double tolerance The tolerance of y(xmax), i.e. a solution is sought   //
//                so that the relative error < tolerance.                     //
//                                                                            //
//  Return Values:                                                            //
//     0   The solution of y' = f(x,y) from x to xmax is stored y[1] and      //
//         h_next has the value to the next size to try.                      //
//    -1   The solution of y' = f(x,y) from x to xmax failed.                 //
//    -2   Failed because either xmax < x or the step size h <= 0.            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
int dopri5(double (*f)(double, double), double *y, double x, double h, double xmax, double *h_next, double tolerance) {

    double scale;
    double temp_y[2];
    double err;
    double yy;
    int i;
    int last_interval = 0;

    // Verify that the step size is positive and that the upper endpoint //
    // of integration is greater than the initial enpoint.               //

    if (xmax < x || h <= 0.0) return -2;

    // If the upper endpoint of the independent variable agrees with the //
    // initial value of the independent variable.  Set the value of the  //
    // dependent variable and return success.                            //

    *h_next = h;
    y[1] = y[0];
    if (xmax == x) return 0;

    // Insure that the step size h is not larger than the length of the //
    // integration interval.                                            //

    h = min(h, xmax - x);

    // Redefine the error tolerance to an error tolerance per unit    //
    // length of the integration interval.                            //

    tolerance /= (xmax - x);

    // Integrate the diff eq y'=f(x,y) from x=x to x=xmax trying to  //
    // maintain an error less than tolerance * (xmax-x) using an     //
    // initial step size of h and initial value: y = y[0]            //

    temp_y[0] = y[0];
    while (x < xmax) {
        scale = 1.0;
        for (i = 0; i < ATTEMPTS; i++) {
            err = fabs(Runge_Kutta(f, temp_y, x, h));
            if (err == 0.0) {
                scale = MAX_SCALE_FACTOR;
                break;
            }
            yy = (temp_y[0] == 0.0) ? tolerance : fabs(temp_y[0]);
            scale = 0.8 * sqrt(sqrt(tolerance * yy / err));
            scale = min(max(scale, MIN_SCALE_FACTOR), MAX_SCALE_FACTOR);
            if (err < (tolerance * yy)) break;
            h *= scale;
            if (x + h > xmax) h = xmax - x;
            else if (x + h + 0.5 * h > xmax) h = 0.5 * h;
        }
        if (i >= ATTEMPTS) {
            *h_next = h * scale;
            return -1;
        };
        temp_y[0] = temp_y[1];
        x += h;
        h *= scale;
        *h_next = h;
        if (last_interval) break;
        if (x + h > xmax) {
            last_interval = 1;
            h = xmax - x;
        }
        else if (x + h + 0.5 * h > xmax) h = 0.5 * h;
    }
    y[1] = temp_y[1];
    return 0;
}


////////////////////////////////////////////////////////////////////////////////
//  static double Runge_Kutta(double (*f)(double,double), double *y,          //
//                                                       double x0, double h) //
//                                                                            //
//  Description:                                                              //
//     This routine uses Prince-Dormand's embedded 4th and 5th order methods  //
//     to approximate the solution of the differential equation y'=f(x,y)     //
//     with the initial condition y = y[0] at x = x0.  The value at x + h is  //
//     returned in y[1].  The function returns err / h ( the absolute error   //
//     per step size ).                                                       //
//                                                                            //
//  Arguments:                                                                //
//     double *f  Pointer to the function which returns the slope at (x,y) of //
//                integral curve of the differential equation y' = f(x,y)     //
//                which passes through the point (x0,y[0]).                   //
//     double y[] On input y[0] is the initial value of y at x, on output     //
//                y[1] is the solution at x + h.                              //
//     double x   Initial value of x.                                         //
//     double h   Step size                                                   //
//                                                                            //
//  Return Values:                                                            //
//     This routine returns the err / h.  The solution of y(x) at x + h is    //
//     returned in y[1].                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

static double Runge_Kutta(double (*f)(double, double), double *y, double x0, double h) {
    static const double two_thirds = 2.0 / 3.0;
    static const double one_seventwoninths = 1.0 / 729.0;
    static const double one_twoninesevenzero = 1.0 / 2970.0;
    static const double one_twofivetwozero = 1.0 / 2520.0;
    static const double one_fiveninefourzero = 1.0 / 5940.0;

    double k1, k2, k3, k4, k5, k6;
    double h5 = 0.2 * h;

    k1 = (*f)(x0, *y);
    k2 = (*f)(x0 + h5, *y + h5 * k1);
    k3 = (*f)(x0 + 0.3 * h, *y + h * (0.075 * k1 + 0.225 * k2));
    k4 = (*f)(x0 + 0.6 * h, *y + h * (0.3 * k1 - 0.9 * k2 + 1.2 * k3));
    k5 = (*f)(x0 + two_thirds * h, *y + one_seventwoninths * h * (226.0 * k1 - 675.0 * k2 + 880.0 * k3 + 55.0 * k4));
    k6 = (*f)(x0 + h, *y + one_twoninesevenzero * h * (-1991 * k1 + 7425.0 * k2
                                                       - 2660.0 * k3 - 10010.0 * k4 + 10206.0 * k5));
    *(y + 1) = *y + one_fiveninefourzero * h * (341.0 * k1 + 3800.0 * k3
                                                - 7975.0 * k4 + 9477.0 * k5 + 297.0 * k6);
    return one_twofivetwozero * (77.0 * k1 - 400.0 * k3 + 1925.0 * k4
                                 - 1701.0 * k5 + 99.0 * k6);
}

