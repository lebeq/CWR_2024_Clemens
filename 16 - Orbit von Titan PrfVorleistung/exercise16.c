#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

typedef void RK4_STEPPER(double,double,double[],double[],ODE_FUNC,int,void*); //t,dt,y[],dy[],ode_func,dim,params

//Konstanten
const double G = 6.674*1e-11; //Gravitationskonst. in m^3/(kg*s^2)
const double mean_v_titan = 5.571*1e3; //Durchschnittl. Geschw. Titan in m/s
const double m_titan = 1.345*1e23; //Masse Titan in kg
const double P_titan = 1.379*1e6; //Periodendauer Titan in s
const double mean_r_titan = 1.222*1e9; //Durchschnittl. Abstand Titan-Saturn (mittelpunkte) in m
const double m_saturn = 5.683*1e26; //Masse Saturn in kg


double vect_betrag_2D(double x, double y){
    return sqrt(x*x+y*y);
}

double vect_betrag_4D(double v[]){
    double val = 0.0;
    for(int i = 0; i < 4; i++){
        val += v[i]*v[i];
    }
    return sqrt(val);
}

//int damit es mit dem typedef in my_numerics übereinstimmt
int ODE_gravity_2D(double t, const double y[], double dy[], void *params){
    //velocity components depend on sign of coord
    //if x is positive v_y i negative etc
    dy[0] = y[2];
    dy[1] = y[3];
    double r = vect_betrag_2D(y[0],y[1]);
    //acceleration components depend on sign of position component
    //if x is positive hte a_x component is negative and vice versa
    dy[2] = (-1.0)*G*(y[0]/r)*(m_saturn)/(r*r);
    dy[3] = (-1.0)*G*(y[1]/r)*(m_saturn)/(r*r);
    return 0;
}

//simulationsschleife, damit main mehr aufgeräumt ist
void simul(double t, double dt, double t_max, double y[], double dy[], ODE_FUNC func, RK4_STEPPER stepper, int dim, int* count, FILE* data, void* params){
    //platzhalter für count
    int dummy_count = *count;
    while(t < t_max){
        //RK4 step
        stepper(t,dt,y,dy,func,dim,NULL);
        //geschw. und abst normieren
        double v_titan_betrag_normiert = vect_betrag_2D(y[2],y[3])/mean_v_titan;
        double r_normiert = vect_betrag_2D(y[0], y[1])/mean_r_titan;
        //Daten in CSV eintragen
        fprintf(data, "%.16f, %.16f, %.16f, %.16f, %.16f\n", t, y[0], y[1], r_normiert, v_titan_betrag_normiert); //t,x,y,r,v
        t += dt;
        dummy_count++;
    }
    //wert von count überschreiben
    *count = dummy_count;
}

int main(){
    /*-------------------------- TEIL 3 ----------------------------*/
    //Speicher für Zustandvektor allokieren
    double *y = malloc(4*sizeof(double));
    double *dy = malloc(4*sizeof(double));
    FILE* data = fopen("titan_data.csv", "w");
    int dim = 4; //größe von Zustandsvektor etc

    //Anfangswerte
    y[0] = (-1.0)*mean_r_titan;
    y[1] = 0.0;
    y[2] = 0.0;
    y[3] = mean_v_titan;

    double delta_t = (1e-3)*P_titan;
    double t = 0.0;
    double t_max = 8.0*P_titan; //8 Orbitalperioden als abbruchbed. für simulation
    int counter = 0; //Schrittzähler für Simulation 

    //Simulation über 8 Orbitalperioden
    simul(t,delta_t,t_max,y,dy,ODE_gravity_2D,mn_rk4_step,dim,&counter,data,NULL);
    
    //Schrittzahl ausgeben
    printf("Simulation, Teil 3, Schrittzahl: %d\n", counter);

    fclose(data);
    free(y);
    free(dy);

    /*---------------------- TEIL 4: start im Aphelion ----------------------------*/
    double new_dt = 1e-5*P_titan;
    double new_t = 0.0;
    //Exzentrizität
    double exc = 0.998;
    //Output File
    FILE* data_e = fopen("exc_titan_data.csv", "w");

    //Speicher für Zustandvekt etc allokieren
    double *new_y = malloc(4*sizeof(double));
    double *new_dy = malloc(4*sizeof(double));
    

    //Anfangswerte - start im Aphelion (Annahme: Saturn ist im rechten Brennpunkt)
    new_y[0] = (-1)*mean_r_titan*(1+exc);
    new_y[1] = 0.0;
    new_y[2] = 0.0;
    new_y[3] = sqrt(G*m_saturn*(1-exc)/(mean_r_titan*(1+exc)));

    int new_counter = 0;
    
    //Simulation über 8 Orbitalperioden - GEHT SCHIEF (soll es auch)
    simul(new_t,new_dt,t_max,new_y,new_dy,ODE_gravity_2D,mn_rk4_step,dim,&new_counter,data_e,NULL);
    
    //Schrittzahl ausgeben
    printf("Simulation, start im Aphelion, ohne adaptiven Zeitschritt Schrittzahl: %d\n", new_counter);

    fclose(data_e);
    free(new_y);
    free(new_dy);
    
    /*--------------------------- TEIL 5 Adaptiver Zeitschritt ------------------------*/
    //CSV für Output
    FILE* adapt_data = fopen("adapt_aphelion.csv", "w");

    //Speicher für Zustände allokieren
    double *adapt_y = malloc(4*sizeof(double));
    double *adapt_dy = malloc(4*sizeof(double));

    //Speicher für Kopie für den adaptiven Zeitschritt
    double *cp_adapt_y = malloc(4*sizeof(double));

    //Anfangswerte für adaptiven Zeitschirtt
    adapt_y[0] = (-1)*mean_r_titan*(1+exc);
    adapt_y[1] = 0.0;
    adapt_y[2] = 0.0;
    adapt_y[3] = sqrt(G*m_saturn*(1-exc)/(mean_r_titan*(1+exc)));

    //counter
    int adapt_counter = 0;

    double adapt_t = 0.0;
    double adapt_delta_t = 1e-5*P_titan;

    double sigma_ok = 1e-2; //gewünschte Genauigkeit
    double epsilon = 0.0;

    while(adapt_t<t_max){
        //Zustand kopieren
        for(int i = 0; i < 4; i++){
            cp_adapt_y[i] = adapt_y[i];
        }
        
        //Original mit dt = adapt_delta_t integrieren
        mn_rk4_step(adapt_t,adapt_delta_t,adapt_y,adapt_dy,ODE_gravity_2D,dim,NULL);
        
        //2 Mal Kopie mit dt = adapt_delta_t/2.0 integrieren
        mn_rk4_step(adapt_t,adapt_delta_t/2.0,cp_adapt_y,adapt_dy,ODE_gravity_2D,dim,NULL);
        mn_rk4_step(adapt_t,adapt_delta_t/2.0,cp_adapt_y,adapt_dy,ODE_gravity_2D,dim,NULL);
        
        //Vektorielle Differenz berechnen
        double sigma = 0.0;
        for(int i = 0; i<4; i++){
            cp_adapt_y[i] = fabs(cp_adapt_y[i] - adapt_y[i]);
        }
        sigma = vect_betrag_4D(cp_adapt_y);
        
        //Korrekturfaktor
        if(sigma == 0.0){
            //dann korrekturfaktor = 1
            epsilon = 1.0;
        }
        else{
            double base = sigma_ok/sigma;
            double expo = 1.0/5.0;
            epsilon = pow(base,expo);
        }

        //Normierte Abstand und Geschw.
        double r_normed = vect_betrag_2D(adapt_y[0], adapt_y[1])/mean_r_titan;
        double v_normed = vect_betrag_2D(adapt_y[2],adapt_y[3])/mean_v_titan;
        
        //Daten inCSV schreiben
        fprintf(adapt_data, "%.16f, %.16f, %.16f, %.16f, %.16f\n", adapt_t, adapt_y[0], adapt_y[1], r_normed, v_normed); //t,x,y,r,v
        
        //dt anpassen
        adapt_delta_t = epsilon*adapt_delta_t;
        adapt_t += adapt_delta_t;
        adapt_counter++;
    }

    //Schrittzahl ausgeben
    printf("Simulation, start im aphelion, adaptiver Zeitschritt Schrittzahl: %d\n", adapt_counter);

    free(adapt_y);
    free(adapt_dy);
    free(cp_adapt_y);
    fclose(adapt_data);
    return 0;
}