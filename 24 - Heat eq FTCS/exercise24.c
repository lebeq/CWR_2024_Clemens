#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

//Aufgabenstellung
const double delta_x = 1e-2;
const double lenght = 1.0;
const double D = 0.1;
const int N = lenght/delta_x;
//Dirichlet Randbedingungen
double T_left_bound = 1.0;
double T_right_bound = -1.0;

void heatFTCS(double t, double y[], double dt){
    double *y_new = malloc(N*sizeof(double));
    
    for(size_t i = 0; i<N; i++){
        double T = y[i%N];
        double TW = y[(i-1)%N];
        double TO = y[(i+1)%N];
        if(i == 0){
            y_new[i] = y[i] + dt*(D * (T_left_bound - 2.0*T + TO))/(delta_x * delta_x);
        }
        else if(i == N-1){
            y_new[i] = y[i] + dt*(D * (TW - 2.0*T + T_right_bound))/(delta_x * delta_x);

        }
        else{
            y_new[i] = y[i] + dt*(D * (TW - 2.0*T + TO))/(delta_x * delta_x);
        }
    }
    //y_new in y speichern für nächsten schritt
    for(size_t i = 0; i<N; i++){
        y[i] = y_new[i];
    }
    free(y_new);
}

int main(){
    FILE* heat_data = fopen("heat_data.csv", "w");
    double *y = malloc(N*sizeof(double));
    //Stabilitätskriterium
    double stability_bound = (delta_x*delta_x)/(2.0*D);
    double t_max = 1.0;

    //Zeitschirtte
    double stable_dt = 1e-4;
    double unstable_dt = 1e-3;

    //Anfangswerte
    for(int i = 0; i<N; i++){
        y[i] = 0.0;
    }

    double t = 0.0;
    size_t j = 0;
    printf("test modulo, N ist %d, mod ist %d, j ist %d\n", N,(j-1)%N, j);
    //simul
    while(t<t_max){
        heatFTCS(t,y,unstable_dt);
        fprintf(heat_data, "%.16f", t);
        for(int i = 0; i<N; i++){
            fprintf(heat_data, ", %.16f", y[i]);
        }
        fprintf(heat_data, "\n");
        t += unstable_dt;
    }

    free(y);
    fclose(heat_data);
    return 0;
}