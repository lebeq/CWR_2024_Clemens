#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

/*----------------- Konstanten -----------------*/
const double my_gamma = 4.0/3.0;
const double my_kappa = 3.85*1e9;
const double G = 6.67384*1e-11; //Gravconst m^3/(kg*s^2)
const double rad_sun = 7*1e8; //Radius of the sun in m
double R = 1.2*rad_sun; //Radius of alpha centauri A in m
double epsilon = 1e-4; //in m


double mass_eps(double rho_0){
    //mass of ball with radius epsilon
    return 4.0/3.0*M_PI*rho_0*pow(epsilon,3.0);
}

//Zustandsvektor ist rho, m
int density_mass_ODE(double r, const double y[], double dy[], void* params){
    dy[0] = -3.0/(2.0*my_kappa*my_gamma) * pow(y[0],2 - my_gamma) * G * y[1]/(r*r);
    dy[1] = 4*M_PI*y[0]*r*r;
    return 0;
}

double sternen_radius(double rho_0){
    double r = epsilon;
    double delta_r = 1e3;
    // init zustand
    double rho[2] = 
    {
        rho_0,
        mass_eps(rho_0),
    };
    double drho[2] = 
    {
        0.0,
        0.0,
    };

    while(rho[0] > 0.1){
        mn_rk4_step(r,delta_r,rho,drho,density_mass_ODE,2,NULL);
        r += delta_r;
    }

    return r;
}

double differenz(double rho_0){
    double rad = sternen_radius(rho_0);
    return fabs(rad - R);
}


int main(){
    double rho_0 = 2*1e4; //in kg/m^3
    double delta_rho = 1e3;

    double res_rho = mn_find_root(differenz, rho_0, delta_rho, 0.01, 100);
    
    //Simulation mit dem gefundenen rho = res_rho
    double r = epsilon;
    double dr = 1e3;
    int counter = 0;
    FILE* density_data = fopen("density_data.csv","w");

    //Anfangswerte
    double rho[2] = {res_rho,mass_eps(res_rho)};
    double drho[2] = {0.0,0.0};    
    
    while(rho[0]>0.1){
        if(counter%10000 == 0){
            fprintf(density_data, "%.16f, %.16f\n",r,rho[0]);
        }
        mn_rk4_step(r,dr,rho,drho,density_mass_ODE,2,NULL);
        r += dr;
        counter++;
    }

    fclose(density_data);
    return 0;
}