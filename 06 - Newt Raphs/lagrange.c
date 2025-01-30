#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

const double mu_earth = 3.986e14;
const double mu_sun = 1.327e20;
const double omega = 1.991e-7;
const double radius_orbit_earth = 1.496e11;

double acceleration(const double x){
    double a_earth = -1*mu_earth/((x - radius_orbit_earth)*(x - radius_orbit_earth));
    double a_sun = -1*mu_sun/(x*x);
    double a_centrifugal = omega*omega*x;
    return a_earth + a_sun + a_centrifugal;
}

int main(){
    double root = mn_find_root(acceleration,1.1*radius_orbit_earth,1000,0.0001,1000);
    double dist_earth = root - radius_orbit_earth;
    printf("%.16e m\n", dist_earth);
    return EXIT_SUCCESS;
}