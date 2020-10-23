#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "float.h"
#include "time.h"
#include "string.h"
#include "levmar.h"
#include "gauss_2d.h"
#include "assert.h"


double rand_double(double _min, double _max) {
    double scale = rand() / (double)RAND_MAX;
    return _min + scale * (_max - _min);
}


int main(int argc, char **argv) {
    // GENERATE an example 2D gaussian with the following hard-coded parameters
    // into an 11x11 pixel grid.
    #define N_GAUSSIAN_2D_PARAMETERS (7)
    double true_p[N_GAUSSIAN_2D_PARAMETERS] = {
        1000.0, // amp
        1.8, // sig_x
        1.2, // sig_y
        5.0, // pos_x
        4.8, // pos_y
        0.05, // rho
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

    printf("PIXELS:\n");
    for(int y=0; y<mea; y++) {
        printf("[ ");
        for(int x=0; x<mea; x++) {
            printf("%4.2f, ", true_pixels[y*mea + x]);
        }
        printf(" ],\n");
    }

    double guess_p[N_GAUSSIAN_2D_PARAMETERS] = {
        1050.0, // amp
        2.0, // sig_x
        2.0, // sig_y
        4.0, // pos_x
        4.0, // pos_y
        0.3, // rho
        40.0, // offset
    };

    double fit_p[N_GAUSSIAN_2D_PARAMETERS];
    memcpy(fit_p, guess_p, sizeof(fit_p));

    double info[10];
    int ret = fit_gauss_2d(noisy_pixels, mea, fit_p, info, NULL);

    if(ret < 0) {
        printf("Levenberg-Marquardt Failed\n");
    }

    printf("Levenberg-Marquardt returned %d in %g iter, reason %g '%s'\n", ret, info[5], info[6], get_dlevmar_stop_reason((int)info[6]));
    printf("Solution:\n");
    for(int i=0; i<N_GAUSSIAN_2D_PARAMETERS; ++i) {
        printf("  %.7g diff=(%.7g, %2.2f%%)\n", fit_p[i], fit_p[i]-true_p[i], (100.0 * (fit_p[i]-true_p[i])) / true_p[i]);
    }
    printf("\n\n");
    printf("Minimization info:\n");
    for(int i=0; i<LM_INFO_SZ; ++i) {
        printf("%g ", info[i]);
    }
    printf("\n");

    return ret;
}
