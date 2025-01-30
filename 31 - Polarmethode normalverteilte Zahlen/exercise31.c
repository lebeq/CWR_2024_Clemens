#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "../MyNumerics/my_numerics.h"

int main(){
    gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(generator, time(NULL));
    Tuple rand_numbers = mn_random_gaussian(generator);
    printf("The random tuple is %.16f, %.16f\n", rand_numbers.x, rand_numbers.y);
    gsl_rng_free(generator);
    return 0;
}

