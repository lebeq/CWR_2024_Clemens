#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

/*
We set mu = plank_const = 1
and m = 0.5
*/

const double mass = 0.5;
const double mu = 1.0;
const double C = 1.74804;

double data_trafo(double x){
    return mn_sign(x)*log10(fabs(x)+1);
}

double E_analytical(int n){
    double factor = M_PI/C;
    return pow(factor*(n-0.5),4.0/3.0); 
}

int schroedinger_harmonic_ODE(double x, const double y[], double f[], void *params){
    double E = *(double *)params;
    f[0] = y[1];
    f[1] = -2.0*mass*(E - mu*pow(x,4.0))*y[0];
    return 0;
}

double test_energy(double E){
    double x_0 = -5.0;
    double x_1 = 0.0;
    double delta_x = 1e-3;
    double psi[2] = {1.0, 0.0};
    double dpsi[2] = {0.0, 0.0};
    double x = x_0;
    while(x + delta_x < x_1){
        mn_rk4_step(x,delta_x,psi,dpsi,schroedinger_harmonic_ODE,2,&E);
        x += delta_x;
    }
    //exakt x_1 treffen sicherstellen
    mn_rk4_step(x,x_1-x,psi,dpsi,schroedinger_harmonic_ODE,2,&E);
    return psi[0]*psi[1];
}

int main(){
    /*
    HIER dann mn_find_root für die verschiedenen 
    dann noch mal simulation mi den n um zu plotten
    */
    FILE* psi_1_data = fopen("psi_1_data.csv", "w");
    FILE* psi_2_data = fopen("psi_2_data.csv", "w");
    FILE* psi_3_data = fopen("psi_3_data.csv", "w");

    double root = 0.0;
    double delta_root = 1.0;
    double n_1 = mn_find_root(test_energy,root,1e-4,1e-10,3000);
    double n_2 = mn_find_root(test_energy,n_1+delta_root,1e-4,1e-10,3000);
    double n_3 = mn_find_root(test_energy,n_2+2*delta_root,1e-4,1e-10,3000);
    printf("the three numerical eigenvals are %.16f, %.16f, %.16f\n",n_1,n_2,n_3);
    printf("the three analytical eigenvals are %.16f, %.16f, %.16f\n",E_analytical(1),E_analytical(2),E_analytical(3));

    //Die zugehörigen Psi_i berechnen und in CSV speichern
    double x_0 = -5.0;
    double x_1 = 0.0;
    double delta_x = 1e-3;
    double psi_1[2]= {0.0,delta_x};
    double psi_2[2] = {0.0,delta_x};
    double psi_3[2] = {0.0,delta_x};
    double dpsi[2] = {0.0,0.0};
    while(x_0 + delta_x < x_1){
        mn_rk4_step(x_0,delta_x,psi_1,dpsi,schroedinger_harmonic_ODE,2,&n_1);
        fprintf(psi_1_data, "%.16f, %.16f\n", x_0, psi_1[0]);
        mn_rk4_step(x_0,delta_x,psi_2,dpsi,schroedinger_harmonic_ODE,2,&n_2);
        fprintf(psi_2_data, "%.16f, %.16f\n", x_0,psi_2[0]);
        mn_rk4_step(x_0,delta_x,psi_3,dpsi,schroedinger_harmonic_ODE,2,&n_3);
        fprintf(psi_3_data, "%.16f, %.16f\n", x_0,psi_3[0]);
        x_0 += delta_x;
    }
    //endpunkt exakt treffen
    mn_rk4_step(x_0,x_1-x_0,psi_1,dpsi,schroedinger_harmonic_ODE,2,&n_1);
    fprintf(psi_1_data, "%.16f, %.16f\n", x_0, psi_1[0]);
    mn_rk4_step(x_0,x_1-x_0,psi_2,dpsi,schroedinger_harmonic_ODE,2,&n_2);
    fprintf(psi_2_data, "%.16f, %.16f\n", x_0, psi_2[0]);
    mn_rk4_step(x_0,x_1-x_0,psi_3,dpsi,schroedinger_harmonic_ODE,2,&n_3);
    fprintf(psi_3_data, "%.16f, %.16f\n", x_0, psi_3[0]);

    fclose(psi_1_data);
    fclose(psi_2_data);
    fclose(psi_3_data);
    return 0;
}