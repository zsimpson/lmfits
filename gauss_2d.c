#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include "time.h"
#include "string.h"
#include "levmar.h"

#ifndef LM_DBL_PREC
#error Demo program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif

// A Gaussian 2D fitter using Levenberg-Marquardt algorithm
// implemented with levmar-2.6. See http://users.ics.forth.gr/~lourakis/levmar/

void gauss_2d(double *p, double *dst_x, int m, int n, void *data) {
    // Arguments:
    //   p: parameters of the 2D Gaussian: array [amp, sig_x, sig_y, pos_x, pos_y, rho, offset]
    //   dst_x: Destination buffer that will contain the function evaluation given the parameters
    //   m: number of parameters (length of p)
    //   n: number of data points
    //   data: data

    double amp = p[0];
    double sig_x = p[1];
    double sig_y = p[2];
    double pos_x = p[3];
    double pos_y = p[4];
    double rho = p[5];
    double offset = p[6];

    double pi2 = 2.0 * M_PI;
    double sgxs = sig_x * sig_x;
    double sgys = sig_y * sig_y;
    double rs = rho * rho;
    double omrs = 1.0 - rs;
    double tem_a = 1.0 / (sig_x * sig_y * omrs);
    double denom = 2.0 * (rho - 1.0) * (rho + 1.0) * sgxs * sgys;
    double numer_const = -2.0 * rho * sig_x * sig_y;
    double linear_term = amp * tem_a * sqrt(omrs);

    int mea = (int)sqrt(n);
    double *dst = dst_x;
    for (int i=0; i<mea; i++) {
        double x = (double)i;
        double xmpx = x - pos_x;
        for (int j=0; j<mea; j++) {
            double y = (double)j;
            double ympy = y - pos_y;
            *dst++ = (
                offset + linear_term * exp(
                    (
                        numer_const * xmpx * ympy
                        + sgxs * ympy * ympy
                        + sgys * xmpx * xmpx
                    ) / denom
                ) / pi2
            );
        }
    }
}

void jac_gauss_2d(double *p, double *dst_jac, int m, int n, void *data) {
    // Arguments:
    //   p: parameters array [amp, sig_x, sig_y, pos_x, pos_y, rho, offset]
    //   dst_jac: destination buffer to hold the Jacobian
    //   m: number of parameters (length of p)
    //   n: number of data points
    //   data: data

    double amp = p[0];
    double sig_x = p[1];
    double sig_y = p[2];
    double pos_x = p[3];
    double pos_y = p[4];
    double rho = p[5];
    double offset = p[6];

    double pi2 = 2.0 * M_PI;
    double sgxs = sig_x * sig_x;
    double sgys = sig_y * sig_y;
    double rs = rho * rho;
    double omrs = 1.0 - rs;
    double tem_a = 1.0 / (sig_x * sig_y * omrs);
    double denom = 2.0 * (rho - 1) * (rho + 1) * sgxs * sgys;
    double linear_term = -2.0 * rho * sig_x * sig_y;
    int mea = (int)sqrt((double)n);

    double *dst = dst_jac;
    for (int i =0; i<mea; i++) {
        double x = (double)i;
        double xmpx = x - pos_x;
        for (int j =0; j<mea; j++) {
            double y = (double)j;
            double ympy = y - pos_y;
            double tem_b = tem_a * sqrt(omrs) * exp(
                (
                    linear_term * xmpx * ympy
                    + sgxs * ympy * ympy
                    + sgys * xmpx * xmpx
                ) / denom
            ) / pi2;

            *dst++ = tem_b;

            double tem_ab_amp = tem_a * tem_b * amp;
            double xmpxy = xmpx * ympy;

            *dst++ = tem_ab_amp * (
                - omrs * sgxs * sig_y
                - rho * sig_x * xmpxy
                + sig_y * xmpx * xmpx
            ) / sgxs;

            *dst++ = - tem_ab_amp * (
                  omrs * sig_x * sgys
                + rho * sig_y * xmpxy
                - sig_x * ympy * ympy
            ) / sgys;

            *dst++ = tem_a * (
                - rho * sig_x * ympy
                + sig_y * xmpx
            ) * tem_b * amp / sig_x;

            *dst++ = - tem_ab_amp * (
                  rho * sig_y * xmpx
                - sig_x * ympy
            ) / sig_y;

            *dst++ = -tem_a * tem_b * (
                rho * (
                    - omrs * sgys
                    + ympy * ympy
                ) * sgxs
                - sig_y * (2.0 - omrs) * xmpxy * sig_x
                + rho * sgys * xmpx * xmpx
            ) * amp / (sig_x * omrs * sig_y);

            *dst++ = 1.0;
        }
    }
}

double dlevmar_opts[LM_OPTS_SZ];
dlevmar_opts[0] = LM_INIT_MU;
    // scale factor for initial mu
dlevmar_opts[1] = 1E-15;
    // stopping threshold for ||J^T e||_inf
dlevmar_opts[2] = 1E-15;
    // stopping threshold for ||Dp||^2
dlevmar_opts[3] = 1E-20;
    // stopping threshold for ||e||^2
dlevmar_opts[4] = LM_DIFF_DELTA;
    // step used in difference approximation to the Jacobian.
    // If delta<0, the Jacobian is approximated  with central differences
    // which are more accurate (but slower!) compared to the forward differences
    // employed by default. Set to NULL for defaults to be used.

char *dlevmar_stop_reasons[] = {
    "Unknown",
    "stopped by small gradient J^T e",
    "stopped by small Dp",
    "stopped by itmax",
    "singular matrix. Restart from current p with increased mu",
    "no further error reduction is possible. Restart with increased mu",
    "stopped by small ||e||_2",
    "stopped by invalid (i.e. NaN or Inf) 'func' values. This is a user error",
};


#define N_GAUSSIAN_2D_PARAMETERS (7)
#define N_MAX_ITERATIONS (1000)

int fit_gauss_2d(double *pixels, int mea, double params[7], double *info, double *covar) {
    /*
    Fit a 2D Gaussian given a square array of pixels with length of size "mea"

    Arguments:
        pixels:
            Square array of pixels (mea x mea)
        mea:
            number of pixels on each side of the square pixels
        params:
            An array [7] of the initial guess of the parameters.
            Will be overwritten with the fit parameters.
            Order:
                amplitude, sigma_x, sigma_y, pos_x, pos_y, rho, offset
        info:
            An array [10] (or NULL if you don't need this info) which contains the fitting info:
                 info[0] = ||e||_2 at initial p.
                 info[1-4] = [ ||e||_2, ||J^T e||_inf,  ||Dp||_2, \mu/max[J^T J]_ii ], all computed at estimated p.
                 info[5] = number of iterations,
                 info[6] = reason for terminating:
                    (See dlevmar_stop_reasons for string constants)
                    0 - Unknown
                    1 - stopped by small gradient J^T e
                    2 - stopped by small Dp
                    3 - stopped by itmax
                    4 - singular matrix. Restart from current p with increased \mu
                    5 - no further error reduction is possible. Restart with increased mu
                    6 - stopped by small ||e||_2
                    7 - stopped by invalid (i.e. NaN or Inf) "func" values; a user error
                 info[7] = number of function evaluations
                 info[8] = number of Jacobian evaluations
                 info[9] = number of linear systems solved, i.e. number of attempts for reducing error
        covar:
            A array [7x7] (or NULL) will be filled in with the fitter's covariance
            estimate on parameters.
    */

    int n_pixels = mea * mea;
    double info[LM_INFO_SZ];
    int ret = 0;

    ret = dlevmar_der(
        gauss_2d, jac_gauss_2d,
        fit_p,
        pixels, N_GAUSSIAN_2D_PARAMETERS, n_pixels,
        N_MAX_ITERATIONS, opts, info,
        NULL,
            // This is a working buffer that will be allocated if NULL
            // I timed it and it made no difference so it seemed easier
            // to let levmar handle malloc/free.
        covar, NULL
    );

    // Without jacobian, useful for sanity checking
    /*
    ret = dlevmar_dif(
        gauss_2d,
        fit_p,
        pixels, N_GAUSSIAN_2D_PARAMETERS, n_pixels,
        N_MAX_ITERATIONS, opts, info,
        NULL,  // See above about working buffer
        covar, NULL
    );
    */

    return ret;
}
