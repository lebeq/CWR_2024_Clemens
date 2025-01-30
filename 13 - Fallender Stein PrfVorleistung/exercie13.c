#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef void ODE_FUNC(double, const double[], double[], int, void*);
typedef void STEPPER(double, double, double[], double[], ODE_FUNC, int, void*);

//am 10.05.24 von https://en.wikipedia.org/wiki/Earth abgelesen
//const double g = 9.80665; //surface gravity in m/s^2 
const double Mass_E = 5.9722*1e24; //Erdmasse in kg
const double Rad_E= 6.3781*1e6; //mean radius in m 

//am 10.5.24 von https://en.wikipedia.org/wiki/Gravitational_constant abgelesen
const double G = 6.67430*1e-11; //Gravitationskonstante in N*m^2*kg^{-2}

//Berechnung von konstanter Gravbeschl.
const double g = G*Mass_E/(Rad_E*Rad_E);

/*
Die Parameter t werden hier ignoriert, weil die DGL nicht von t abhängt, i.A. kann das sein, deswegen ist es mit drinnen
Das t wird an euler_step übergeben, dort nicht genutzt, vllt gibt es warning
*/

void ODE_falling_stone_const(double t, const double y[], double  dy[], int dim, void *params){
    /*
    dy/dt = (v_y,a) mit y=(y,v_y) und a = -g const
    */
    dy[0] = y[1];
    dy[1] = -g;
}

void ODE_falling_stone_varying(double t, const double y[], double dy[], int dim, void *params){
    /*
    \vect{y} = (y,dy), dann ist die Vektor DGL 1. Ordnung
    dy/dt = (dy, a); mit a = -GM_E/(y^2)
    */
   dy[0] = y[1];
   dy[1] = -G*Mass_E/(y[0]*y[0]);
}

void euler_step(double t, double dt, double y[], double dy[], ODE_FUNC ode_func, int dim, void* params){
    /*
    y_{i+1}[0] = y[0] + dy[0]*dt
    y_{i+1}[1] = y[1] + dy[1]*dt
    wo dy von ode_func berechnet wird 
    */
   ode_func(t,y,dy,2,params); 
   for(int i = 0; i<dim;i++){
    y[i] = y[i] + dy[i]*dt;
   }
}

double v_ana_const(double h){
    return sqrt(2.0*g*h);
}

double v_ana_var(double h){
    if(h > 0.0){
        return sqrt(2.0*g*h*Rad_E/(Rad_E+h));
    }
    else{
        //wenn h < Rad_E, dann sollte es nicht möglich sein, v = 0 da das Teil auf der Oberfläche liegt
        return 0.0;
    }
}

//Methode für simulation, führt euler_step aus bis Stein auf Erdoberfläche ist
void simul(double t, double dt, double y[], double dy[], ODE_FUNC ode_func, STEPPER stepper, int dim, void* params){
    while(y[0] >= Rad_E){
        stepper(t,dt,y,dy,ode_func,dim,params);
        t += dt;
    }
}

int main(){

    /*------------------- TEIL 5 -------------------------*/
    printf("----------------AUFGABENTEIL 5-------------------\n");
    //speicher für y allokieren
    double *y_const = malloc(2*sizeof(double));
    double *dy_const = malloc(2*sizeof(double));
    double *y_var = malloc(2*sizeof(double));
    double *dy_var = malloc(2*sizeof(double));

    //Fallhöhe 200km = 200*1e3 m
    double h1 = 200*1e3;

    //Anfagnswerte
    y_const[0] = Rad_E + h1;
    y_const[1] = 0;
    y_var[0] = Rad_E + h1;
    y_var[1] = 0; 

    //Die dy Arrays initialisieren mit 0,0
    dy_const[0] = 0.0;
    dy_const[1] = 0.0;
    dy_var[0] = 0.0;
    dy_var[1] = 0.0;   

    //Zeitschritt und Zeit initialisieren
    double dt = 1e-6;
    double t = 0.0;

    //Simluation für const Beschl.
    simul(t,dt,y_const,dy_const,ODE_falling_stone_const,euler_step,2,NULL);

    //Simulation für var Beschl
    simul(t,dt,y_var,dy_var,ODE_falling_stone_varying,euler_step,2,NULL);


    //Ausgabe von Gechiwndigkeiten und deren Ratio
    printf("const a; v at impact: %.16f\n",fabs(y_const[1]));
    printf("const a; v analytical at impact is: %.16f\n", v_ana_const(h1));
    printf("var a; v at impact: %.16f\n", fabs(y_var[1]));
    printf("var a; v_analytical at impact is: %.16f\n", v_ana_var(h1));
    printf("ratio von v_const/v_var: %.16f\n", y_const[1]/y_var[1]);
    

    free(y_const);
    free(dy_const);
    free(y_var);
    free(dy_var);

    /*------------------- TEIL 6 -----------------------*/
    printf("----------------AUFGABENTEIL 6-------------------\n");
    //gleichverteilte exponenenten berechnen
    double expo_s = -1.0;
    double expo_e = -7.0;
    double N = 600.0;
    double delta_expo = fabs(expo_e - expo_s)/N;

    //Fallhöhe 20m
    double h = 20.0;

    //Speicher allokieren
    double *ny_const = malloc(2*sizeof(double));
    double *ny_var = malloc(2*sizeof(double));

    //Dateien anlegen
    FILE* data_const5 = fopen("p5_euler_const_g.csv", "w");
    FILE* data_var5 = fopen("p5_euler_var_g.csv", "w");
    FILE* abs_err_const = fopen("p5_abs_err_const.csv", "w");
    FILE* abs_err_var = fopen("p5_abs_err_var.csv", "w");


    for(int i = 0; i<600; i++){

        //es ist -i*delta_expo, weil die exponenetnen alle negativ sind 
        double dt = pow(10.0,expo_s - i*delta_expo);

        //Anfangswerte
        ny_const[0] = Rad_E + h; //Fallen ua Höhe h=20m
        ny_const[1] = 0.0;
        ny_var[0] = Rad_E + h;
        ny_var[1] = 0.0;
        double t = 0.0;
        double *dy_const = malloc(2*sizeof(double));
        double *dy_var = malloc(2*sizeof(double));
        
        //Anfangswerte für dy
        dy_const[0] = 0.0;
        dy_const[1] = 0.0;
        dy_var[0] = 0.0;
        dy_var[1] = 0.0;

        //Simulationen laufen
        simul(t,dt,ny_const,dy_const,ODE_falling_stone_const,euler_step,2,NULL);
        simul(t,dt,ny_var,dy_var,ODE_falling_stone_varying,euler_step,2,NULL);
        
        //Betrag von Fehler berechnen
        double err_const = fabs(fabs(ny_const[1])-v_ana_const(h));
        double err_var = fabs(fabs(ny_var[1])-v_ana_var(h));
        
        //CSV Dateien beschriften
        fprintf(abs_err_const, "%.16f, %.16f\n", dt, err_const);
        fprintf(abs_err_var, "%.16f, %.16f\n", dt, err_var);
        
        free(dy_const);
        free(dy_var);
    }
    printf("DONE\n");
    fclose(data_const5);
    fclose(data_var5);
    fclose(abs_err_const);
    fclose(abs_err_var);
    free(ny_const);
    free(ny_var);
    return EXIT_SUCCESS;
}