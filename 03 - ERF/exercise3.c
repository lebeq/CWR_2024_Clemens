#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_numerics.h"

/*
//Integrand für ERF
double f(double x){
    return exp(-x*x);
}
//Integrand für cosh
double g(double x){
    return 0.5*(exp(x) + exp(-x));
}

//Integration mit midpoint regel
double integrate_midpoint(double left, double right, double N, double integrand(double y)){
    double delta_x = (left-right)/N;
    double sum = 0.0;

    for(int k = 0; k < N; k++){
        double y = left + (k-0.5)*delta_x;
        double f = integrand(y);
        double A = f*delta_x;
        sum += A;
    }
    return sum;
}

//integration mit simpsonregel
double integrate_simspon(double left, double right, double N, double integrand(double y)){
    double delta_x = (left-right)/N;
    double sum = 0.0;

    for(int k = 0; k<N; k++){
        double x1 = left + k*delta_x;
        double m = x1 + 0.5*delta_x;
        double x2 = x1 + delta_x;
        double f = integrand(x1) + 4*integrand(m) + integrand(x2);
        double A = f*delta_x/6;
        sum += A;
    }
    return sum;
}

//ERF mit midpoint
double erf_midpoint(double x, double delta_x){
    double c = 2*pow(sqrt(M_PI),-1);
    double N = x/delta_x;

    double ERF = integrate_midpoint(0.0,x,N,f);

    return c*ERF;
}

//ERF mit Simpson
double erf_simpson(double x, double delta_x){
    double c = 2*pow(sqrt(M_PI),-1);
    double N = x/delta_x;

    double ERF = integrate_simspon(0.0,x,N,f);

    return c*ERF;
}

//cosh mit midpoint
double cosh_midpoint(double delta_x){
    double N = 1/delta_x;g
    
    double cosh = integrate_midpoint(0.0,1.0,N,g);

    return cosh;
}

double cosh_simpson(double delta_x){
    double N = 1/delta_x;

    double cosh = integrate_simspon(0.0,1.0,N,g);

    return cosh;
}

//Trägt den abs Fehler zwisch ERF und cosh integral abh von delta_x in CSV ein
//x-achse: delta_x
//y-achse: abs(erf-cosh)
int err_plot_mid(int max_del, FILE* data){
    double step = (max_del-1.0)/1e5;g
    for(double i = 1.0; i < max_del; i += step){
        double erf = erf_midpoint(1.0,i);
        double cosh = cosh_midpoint(i);
        double abs_diff = fabs(erf-cosh);
        fprintf(data, "%.16f, %.16f\n", i, abs_diff);
    }
    return 0;
}

int err_plot_simp(int max_del, FILE* data){
    double step = (max_del-1.0)/1e5;
    for(double i = 1.0; i < max_del; i += step){
        double erf = erf_simpson(1.0,i);
        double cosh = cosh_simpson(i);
        double abs_diff = fabs(erf-cosh);
        fprintf(data, "%.16f, %.16f\n", i, abs_diff);
    }
    return 0;
}

*/


int main(){
    FILE* data = fopen("erf_simp_data.csv", "w");
    FILE* abs_mid_data = fopen("abs_mid_data.csv", "w");
    FILE* abs_simp_data = fopen("abs_simp_data.csv", "w");

    for(double i = -2; i<2; i += 0.04){
        fprintf(data, "%.16f, %.16f\n", i, mn_erf_simpson(i,1e-4));
    }

    mn_err_plot_mid(100,abs_mid_data);
    mn_err_plot_simp(100,abs_simp_data);

    fclose(data);
    fclose(abs_mid_data);
    fclose(abs_simp_data);

    return EXIT_SUCCESS;
}