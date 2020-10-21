#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "levmar.h"
#ifndef LM_DBL_PREC
#error Demo program assumes that levmar has been compiled with double precision, see LM_DBL_PREC!
#endif


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
    double pi2 = M_PI * 2.0;
    double sgxs = sig_x * sig_x;
    double sgys = sig_y * sig_y;
    double rs = rho * rho;
    double omrs = 1.0 - rho * rho;
    double tem_a = 1.0 / (sig_x * sig_y * omrs);
    double const_numerator = -2.0 * rho * sig_x * sig_y;
    double denominator = 2.0 * (rho - 1.0) * (rho + 1.0) * sgxs * sgys;
    int mea = (int)sqrt((double)n);
    int k = 0;
    for(int i=0; i<mea; i++) {
        double x = (double)i;
        double xmpx = x - pos_x;
        for(int j=0; j<mea; j++) {
            double y = (double)j;
            double ympy = y - pos_y;
            dst_x[k++] = (
                offset
                + amp * tem_a * sqrt(omrs) * exp(
                    (
                        const_numerator
                        + sgxs * ympy * ympy
                        + sgys * xmpx * xmpx
                    ) / denominator
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

    double sgxs = sig_x * sig_x;
    double sgys = sig_y * sig_y;
    double rs = rho * rho;
    double omrs = 1.0 - rs;
    double tem_a = 1.0 / (sig_x * sig_y * omrs);
    int mea = (int)sqrt((double)n);
    int k = 0;
    for (int i =0; i<mea; i++) {
        double x = (double)i;
        double xmpx = x - pos_x;
        for (int j =0; j<mea; j++) {
            double y = (double)j;
            double ympy = y - pos_y;
            double tem_b = tem_a * sqrt(omrs) * exp(
                (
                    -2.0 * rho * sig_x * sig_y * xmpx * ympy
                    + sgxs * ympy * ympy
                    + sgys * xmpx * xmpx
                ) / (
                    2.0 * (rho - 1) * (rho + 1) * sgxs * sgys
                )
            ) / (2.0 * M_PI);

            dst_jac[k++] = tem_b;

            double pr = amp * tem_b;

            dst_jac[k++] = tem_a * tem_b * amp * (
                - omrs * sgxs * sig_y
                - rho * sig_x * xmpx * ympy
                + sig_y * xmpx * xmpx
            ) / sgxs;

            dst_jac[k++] = -tem_a * tem_b * amp * (
                  omrs * sig_x * sgys
                + rho * sig_y * xmpx * ympy
                - sig_x * ympy * ympy
            ) / sgys;

            dst_jac[k++] = tem_a * (
                - rho * sig_x * ympy
                + sig_y * xmpx
            ) * tem_b * amp / sig_x;

            dst_jac[k++] = -tem_a * tem_b * amp * (
                  rho * sig_y * xmpx
                - sig_x * ympy
            ) / sig_y;

            dst_jac[k++] = -tem_a * tem_b * (
                rho * (
                    - omrs * sgys
                    + ympy * ympy
                ) * sgxs
                - sig_y * (2.0 - omrs) * ympy * xmpx * sig_x
                + rho * sgys * xmpx * xmpx
            ) * amp / (sig_x * omrs * sig_y);

            dst_jac[k++] = 1.0;
        }
    }
}

double rand_double(double _min, double _max) {
    double scale = rand() / (double)RAND_MAX;
    return _min + scale * (_max - _min);
}

int main(int argc, char **argv) {
    /*
    TODO
    if(argc == 1) {
        printf("lmfitter\n");
        printf("  mode: String[dif|der]\n");
        printf("        dif=finite difference mode\n");
        printf("        der=analytic derivative mode\n");
        return 1;
    }
    */

    // GENERATE an example 2D gaussian with the following hard-coded parameters
    // into an 11x11 pixel grid.
    #define N_GAUSSIAN_2D_PARAMETERS (7)
    double true_p[] = {
        100.0, // amp
        1.8, // sig_x
        1.2, // sig_y
        5.0, // pos_x
        4.8, // pos_y
        0.5, // rho
        50.0, // offset
    };
    int mea = 11;
    int n_pixels = mea * mea;
    double *true_pixels = (double *)alloca(sizeof(double) * n_pixels);
    gauss_2d(true_p, true_pixels, N_GAUSSIAN_2D_PARAMETERS, n_pixels, true_pixels);

    // ADD some noise to the true_pixels to make noisy_pixels
    double *noisy_pixels = (double *)alloca(sizeof(double) * n_pixels);
    for(int i=0; i<n_pixels; i++) {
        noisy_pixels[i] = true_pixels[i] + rand_double(-0.1, 0.1);
    }

    // START the fitter some offset form the true parameters
    // Note these initial guess parameters will be overwritten with the fit values
    double guess_p[] = {
        110.0, // amp
        2.0, // sig_x
        2.0, // sig_y
        4.0, // pos_x
        4.0, // pos_y
        0.3, // rho
        40.0, // offset
    };

    // CALL the LM fitter on the noisy_pixels
    double opts[LM_OPTS_SZ], info[LM_INFO_SZ];
    opts[0] = LM_INIT_MU;
        // scale factor for initial mu
    opts[1] = 1E-15;
        // stopping threshold for ||J^T e||_inf
    opts[2] = 1E-15;
        // stopping threshold for ||Dp||^2
    opts[3] = 1E-20;
        // stopping threshold for ||e||^2
    opts[4] = LM_DIFF_DELTA;
        // step used in difference approximation to the Jacobian.
        // If delta<0, the Jacobian is approximated  with central differences
        // which are more accurate (but slower!) compared to the forward differences
        // employed by default. Set to NULL for defaults to be used.

    #define N_MAX_ITERATIONS (1000)
    int ret = 0;

    // Without jacobain
    ret = dlevmar_dif(
        gauss_2d,
        guess_p,
        noisy_pixels, N_GAUSSIAN_2D_PARAMETERS, n_pixels,
        N_MAX_ITERATIONS, opts, info, NULL, NULL, NULL
    );

    // With jacobian
    /*
    ret = dlevmar_der(
        gauss_2d, jac_gauss_2d,
        guess_p,
        noisy_pixels, N_GAUSSIAN_2D_PARAMETERS, n_pixels,
        N_MAX_ITERATIONS, opts, info, NULL, NULL, NULL
    );
    */

    char *reasons[] = {
        "Unknown",
        "stopped by small gradient J^T e",
        "stopped by small Dp",
        "stopped by itmax",
        "singular matrix. Restart from current p with increased mu",
        "no further error reduction is possible. Restart with increased mu",
        "stopped by small ||e||_2",
        "stopped by invalid (i.e. NaN or Inf) 'func' values. This is a user error",
    };

    if(ret < 0) {
        printf("Levenberg-Marquardt Failed\n");
    }

    printf("Levenberg-Marquardt returned %d in %g iter, reason %g '%s'\n", ret, info[5], info[6], reasons[(int)info[6]]);
    printf("Solution:\n");
    for(int i=0; i<N_GAUSSIAN_2D_PARAMETERS; ++i) {
        printf("%.7g ", guess_p[i]);
    }
    printf("\n\n");
    printf("Minimization info:\n");
    for(int i=0; i<LM_INFO_SZ; ++i) {
        printf("%g ", info[i]);
    }
    printf("\n");

    return 0;
}
