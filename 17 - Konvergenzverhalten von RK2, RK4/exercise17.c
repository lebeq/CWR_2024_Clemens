#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

//Konstanten für Federpendel, N_=1
const int N = 1; //Anzahl von Pendeln
const double k = 100; //Federkonstante in N/m
const double mass = 1; //masse in kg
const double base_length = 0.0; //Anfangsauslenkung in m
const double v_0 = 1.0; //Anfangsgeschwindigkeit in m/s
const double t_max = 20.0;
double delta_t = 1e-3; //Zeitschrittweite

//DGL Funktion, aus '14 - GSL Pendel/pendulums.c'
int pendulums_ODE(double t, const double y[], double f[], void *params)
{
    for(int i = 0; i<2*N; i++){
        f[i] = 0.0;
    }
    for(int i = 0; i < 2*N; i++){
        if(i < N){
            f[i] = y[N+i];
        }
        if(i == N){
            f[N] += -k*y[0]/mass;
        }
        if(i > N){
            double F = -k*(y[i-N]-y[i-N-1]-base_length)/mass;
            f[i-1] -= F; //die vorherige Masse spürt 2 Kräfte oder so
            f[i] += F;
        }
    }
    return EXIT_SUCCESS;
}

double pendulums_energy(const double y[]){
    //Kinetische energie
    double E_kin = 0.0;
    //for(int i = 0; i<N; i++){
    //    E_kin += mass*y[N+i]*y[N+i]/2;
    //}
    E_kin = mass*y[1]*y[1]/2.0; 
    //Potentielle Energie
    double E_pot = 0.0;
    //Edgecase von 1-tem Teilchen, es hat am Anfang keine Auslenkung
    E_pot = k*((y[0])*(y[0]))/2;
    //for(int i = 1; i<N; i++){
    //    E_pot += k*(y[i]-y[i-1]-base_length)*(y[i]-y[i-1]-base_length)/2;
    //}
    return E_kin + E_pot;
}

int main(){
    //Speicher allokieren
    double *y = malloc(2*sizeof(double));
    double *dy = malloc(2*sizeof(double));
    FILE *pos_file = fopen("pos_file.csv", "w");
    FILE *energy_file = fopen("energy_file.csv", "w");
    FILE *res_file = fopen("residuum_file.csv", "w");

    //Periode
    double period = 2.0*M_PI*sqrt(mass/k);
    //viertelperioden Position analytiusche LSG
    double x_t4 = v_0*sqrt(mass/k);
    printf("viertelperiodenpos iost %.16f\n", x_t4);


    //Startwerte
    y[0] = 0.0;
    y[1] = 1.0;


    //log gleichverteilte dt werte
    double t = 0.0;
    double expo_s = 0.0;
    double expo_e = -4.0;
    double dexpo = fabs(expo_e - expo_s)/100.0;
    printf("dexpo %.16e\n", dexpo);

    for(int i = 0; i<100; i++){
        delta_t = pow(10.0, expo_s - i*dexpo);
        printf("delta_t is %.16e\n", delta_t);
        //Startwerte
        y[0] = 0.0;
        y[1] = 1.0;
        dy[0] = 0.0;
        dy[1] = 0.0;
        double t = 0.0;

        //Simulationsschleif, bis T/4
        while(t < period/4.0){
            //mn_rk2_step(t,delta_t,y,dy,pendulums_ODE,2*N,NULL);
            //mn_rk4_step(t,delta_t,y,dy,pendulums_ODE,2*N,NULL);
            //mn_euler_step(t,delta_t,y,dy,pendulums_ODE,2*N,NULL);
            mn_velocity_verlet_step(t,delta_t,y,dy,pendulums_ODE,2*N,NULL);
            double E_tot = pendulums_energy(y);
            //fprintf(pos_file, "%.16f, %.16f\n", t, y[0]);
            //fprintf(energy_file, "%.16f, %.16f\n", t, E_tot);
            if(t+delta_t > period/4.0){
                delta_t = period/4.0 - t;
                t += delta_t;
            }
            else{
                t += delta_t;
            }
        }
        //double ampl = v_0*sqrt(mass/k);
        //double analytic_sol = ampl*sin(sqrt(mass/k)*(period/4.0));
        double res = fabs(y[0]-x_t4) +1e-9;
        fprintf(res_file, "%.16e, %.16e \n", delta_t, res);
    }
        

    fclose(pos_file);
    fclose(energy_file);
    fclose(res_file);
    free(y);
    free(dy);
    return 0;
}