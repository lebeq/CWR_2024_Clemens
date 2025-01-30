#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

const double a = 1.0/6.1;
const double mygamma = 1.0/3.7;
const double N = 83*1e6;


//RHS von DGL für S
//double S_RHS(double delta_x, double)

//ein Step von Euler für S(x) 
double S_step(double delta_x, double I, double R_0, double S_0){
    //S_{i+1} = (-gamma*R_0*I_{i}*S_{i}/N)*delta_x + S_i
    return (-mygamma*R_0*I*S_0*delta_x)/N + S_0;
} 

//ein Step von Euler für E(x)
double E_step(double delta_x, double I, double R_0, double S, double E){
    //E_{i+1} = (gamma*R_0*I_{i}*S_{i}/N + a*E_{i})*delta_x + E_{i}
    return (mygamma*R_0*I*S/N - a*E)*delta_x +E;
}

//ein step von Euler  für I
double I_step(double delta_x, double E, double I){
    return (a*E - mygamma*I)*delta_x + I;
}

//ein step von Euler für R
double R_step(double delta_x, double I, double R){
    return mygamma*I*delta_x + R;
}

int main(){
    //Anfangswerte
    double S_0 = 83*1e6; //Anfangswert von Susceptible
    double E_0 = 3*1e4; //Anfangswert von Exposed
    double I_0 = 9*1e3; //Anfangswert von Infectious
    double R_anf = 0.0; //Anfangswert von Removed
    double delta = 0.01; //delta t
    double R0 = 2;  //Basisreproduktionszahl
    double t_max = 300.0; //max laufzeit der Simulation

    
    FILE* data = fopen("SEIRdata.csv", "w");

    double S = S_0;
    double E = E_0;
    double I = I_0;
    double R = R_anf;

    for(double t = 0.0; t < t_max; t+= delta){
        fprintf(data, "%f, %.16f, %.16f, %.16f, %.16f\n", t,S,E,I,R);
        S = S_step(delta,I,R0,S);
        E = E_step(delta,I,R0,S,E);
        I = I_step(delta,E,I);
        R = R_step(delta,I,R);    
    }


    fclose(data);
    return EXIT_SUCCESS;
}