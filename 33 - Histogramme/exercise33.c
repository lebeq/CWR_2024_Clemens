#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>

#include "../MyNumerics/my_numerics.h"

/*--------------- Params --------------------*/

const double bin_number = 50.0;
const double min = -5.0;
const double max = 5.0;

/*------------------------ gsl_rng setup and co --------------------------*/

gsl_rng* get_gsl_gen(){
    gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(generator, time(NULL));
    return generator;
}

void free_gsl(gsl_rng* generator){
    if(generator != NULL){
        gsl_rng_free(generator);
    }
}

double random_uniform(gsl_rng* generator, double min, double max){
    double extent = max - min;
    double r = gsl_rng_uniform(generator)*extent;
    return r + min;
}



int main(){
    FILE* hist_data = fopen("hist_data.csv", "w");
    double i = 0.0;
    double i_max = 5e6; //wir wollen 1e7 zahlen, aber mn_random_gaussian generiert 2, also brauchen wir nur die hälfte an durchlaufen

    // Hiatogramm erzeugen
    Hist* hist = hist_malloc(bin_number,min,max);
    hist_reset(hist);

    //gsl_histogram* hist = gsl_histogram_alloc(bin_number);
    //gsl_histogram_set_ranges_uniform(hist, min, max);


    gsl_rng* r = get_gsl_gen();
    while(i<i_max){
        //get tuple of rand numbers using polarmenthode
        Tuple res = mn_random_gaussian(r);
        //update histogram
        hist_add(hist, res.x);
        //hist_add(hist, res.y);
        i++;
    }

    //printf to file
    hist_fprintf_rel(hist, hist_data);
    printf("Erwartungswert von dem Histogramm ist %.16f\n", hist_mean(hist));
    printf("Standardabweichung von dem Histogramm ist %.16f\n", hist_std(hist));

    //calculate mean and st. dev
    /*double mean = gsl_histogram_mean(hist);
    double std = gsl_histogram_sigma(hist);

    printf("Mean of histogram is %g\n", mean);
    printf("Standard deviaton of histogram is %g\n", std);*/

    hist_free(hist);
    free_gsl(r);
    fclose(hist_data);
    return 0;
}