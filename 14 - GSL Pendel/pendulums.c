/*---------------------------------------------*\
|  CWR 2024                                     |
|  Blatt 4 - Aufgabe 14                         |
|  Original Author: mail@saschalambert.de       |
\*---------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

/* TODO */
/*------------------------  PHYSIKALISCHE KONSTANTEN  -------------------------*/

const int N = 1 ;              // Anzahl der Pendel
const double k = 100.0 ;           // Federhaerte [N/m]
const double base_length = 1.0 ; // Basislaenge der Federn [m]
const double mass = 1.0 ;        // Masse der Pendel [kg]

/*-------------------------  SIMULATIONS-PARAMETER  ---------------------------*/

const double T_max = 20.0;
const double delta_t = 1e-3; // Zeitliche Schrittweite

//static const char pos_file_name[] = "position.dat";
//static const char energy_file_name[] = "energy.dat";

/*-------------------------  PHYSIKALISCHES SYSTEM  ---------------------------*/

/* TODO */
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
    return GSL_SUCCESS;
}

/* TODO */
/* Bearbeiten Sie hier Aufgabe 10.9 */
double pendulums_energy(const double y[]){
    //Kinetische energie
    double E_kin = 0.0;
    for(int i = 0; i<N; i++){
        E_kin += mass*y[N+i]*y[N+i]/2;
    }
    //Potentielle Energie
    double E_pot = 0.0;
    //Edgecase von 1-tem Teilchen, es hat am Anfang keine Auslenkung
    E_pot += k*(y[0]*y[0])/2;
    for(int i = 1; i<N; i++){
        E_pot += k*(y[i]-y[i-1]-base_length)*(y[i]-y[i-1]-base_length)/2;
    }
    return E_kin + E_pot;
}

/*-----------------------------------------------------------------------------*/

int main(void)
{
    /* ------- GSL Verwaltung -------------------------------------------------*/

    /* TODO */
    /* Berechnen Sie hier die korrekte Dimensionalitaet des DGL-Systems */
    int gsl_dimension = 2*N;

    /* Initialisieren des Systems fuer GSL
       Wir uebergeben Ihre erstellte Funktion und die Dimensionaliaet */
    gsl_odeiv2_system SYSTEM = {pendulums_ODE, NULL, gsl_dimension, NULL};

    /* Auswahl des Integrators
       Wir waehlen hier einen Runge-Kutta 4. Ordnung (rk4). */
    gsl_odeiv2_step *STEPPER =
        gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk4, gsl_dimension);

    /* ------- Anfangswerte ---------------------------------------------------*/

    /* TODO */
    /* Erstellen Sie zwei Arrays mit der notwendigen Groesse. Sie muessen
       alle Positionen und Geschwindigkeiten fassen koennen. */
    double y[gsl_dimension];
    double yerr[gsl_dimension]; /* GSL benoetigt yerr */

    /* TODO */
    /* Fuellen Sie das Array y[] mit den Startwerten */
    for(int i = 0; i < N; i++){
        y[i] = i*base_length;
    }
    y[N] = 20.0; //v_0 = 20 m/s
    for(int i = N+1; i<gsl_dimension; i++){
        y[i] = 0.0;
    }

    /* ------- Simulations-Schleife -------------------------------------------*/

    /* Ausgabe-Dateien oeffnen */
    FILE *pos_file = fopen("pos_file.csv", "w");
    FILE *energy_file = fopen("energy_file.csv", "w");

    /* Simulationszeit */
    double t = 0;
    
    //Speicher für Euler allokieren
    double *dy = malloc(gsl_dimension*sizeof(double));

    //printf("warum?\n");
    /* #####   HAUPTSCHLEIFE   ##### */
    while (t < T_max)
    {
        /* step_apply befoertdert y[] zum naechsten Zeitpunkt */
        //gsl_odeiv2_step_apply(STEPPER, t, delta_t,
        //                      y, yerr, NULL, NULL, &SYSTEM);
        //mn_euler_step(t,delta_t,y,dy,pendulums_ODE,gsl_dimension,NULL);
        //mn_rk4_step(t,delta_t,y,dy,pendulums_ODE,gsl_dimension,NULL);
        mn_velocity_verlet_step(t,delta_t,y,dy,pendulums_ODE,gsl_dimension,NULL);
        double E_tot = pendulums_energy(y);

        /* TODO */
        /* Geben Sie hier die Daten in die Dateien aus */
        //POSITIONEN
        fprintf(pos_file,"%.16f",t);
        for(int i = 0; i<N;i++){
            fprintf(pos_file, ", %.16f ", y[i]);
        }
        fprintf(pos_file, "\n");
        //Energie
        fprintf(energy_file, "%.16f, %.16f\n", t, E_tot);

        t += delta_t;
    }

    fclose(pos_file);
    fclose(energy_file);
    free(dy);
    return EXIT_SUCCESS;
}
