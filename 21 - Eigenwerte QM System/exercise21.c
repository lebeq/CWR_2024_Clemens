#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

/*
We set m=k=plank_const = 1
*/

double data_trafo(double x){
    return mn_sign(x)*log10(fabs(x)+1);
}

int schroedinger_harmonic_ODE(double x, const double y[], double f[], void *params){
    double E = *(double *)params;
    f[0] = y[1];
    f[1] = -2.0*(E - (x*x)/2.0)*y[0];
    return 0;
}

double test_energy(double E){
    double x_0 = -10.0;
    double x_1 = 0.0;
    double delta_x = 1e-3;
    double psi[2] = {1.0, 0.0};
    double dpsi[2] = {0.0, 0.0};
    double x = x_0;
    while(x < x_1){
        mn_rk4_step(x,delta_x,psi,dpsi,schroedinger_harmonic_ODE,2,&E);
        x += delta_x;
    }
    //exakt x_1 treffen sicherstellen
    mn_rk4_step(x,x_1-x,psi,dpsi,schroedinger_harmonic_ODE,2,&E);
    return psi[0]*psi[1];
}

int main(){
    //Datei für Output
    FILE* energy_data = fopen("energy_data.csv", "w");
    FILE* psi_0_data = fopen("psi_0_data.csv", "w");
    FILE* psi_1_data = fopen("psi_1_data.csv", "w");

    //Energieintervall
    double E_0 = 0.1;
    double E_1 = 10.0;
    double delta_E = 1e-3;
    double E = 0.0;
    
    //Energien plotten
    double energy_output = 0.0;
    for(int i = 0; i < 9.9*1e3; i++){
        E = E_0 + i*delta_E;
        //energie berechnene und Transformieren
        energy_output = data_trafo(test_energy(E));
        fprintf(energy_data, "%.16f, %.16f\n", E, energy_output);
    }

    //Nulstellensuche
    double n_0 = mn_find_root(test_energy, 0.0, 1e-3, 1e-10, 300);
    printf("die erste Energie Nullstelle ist %.16f\n", n_0);
    double n_1 = mn_find_root(test_energy, 0.7, 1e-3, 1e-10, 300);
    printf("die zweie Energie Nullstelle ist %.16f\n", n_1);

    //Die zugehörigen Psi_i berechnen und in CSV speichern
    double x_0 = -10.0;
    double x_1 = 0.0;
    double delta_x = 1e-3;
    double psi_0[2]= {0.0,delta_x};
    double psi_1[2] = {0.0,delta_x};
    double dpsi[2] = {0.0,0.0};
    while(x_0 < x_1){
        mn_rk4_step(x_0,delta_x,psi_0,dpsi,schroedinger_harmonic_ODE,2,&n_0);
        fprintf(psi_0_data, "%.16f, %.16f\n", x_0, psi_0[0]);
        mn_rk4_step(x_0,delta_x,psi_1,dpsi,schroedinger_harmonic_ODE,2,&n_1);
        fprintf(psi_1_data, "%.16f, %.16f\n", x_0,psi_1[0]);
        x_0 += delta_x;
    }
    //endpunkt exakt treffen
    mn_rk4_step(x_0,x_1-x_0,psi_0,dpsi,schroedinger_harmonic_ODE,2,&n_0);
    fprintf(psi_0_data, "%.16f, %.16f\n", x_0, psi_0[0]);
    mn_rk4_step(x_0,x_1-x_0,psi_1,dpsi,schroedinger_harmonic_ODE,2,&n_1);
    fprintf(psi_1_data, "%.16f, %.16f\n", x_0, psi_1[0]);


    fclose(energy_data);
    fclose(psi_0_data);
    fclose(psi_1_data);
    return 0;
}