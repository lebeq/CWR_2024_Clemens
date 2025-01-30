#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

/*--------------- Konstanten -----------------*/
//Gravitationskonstante
const double G = 6.674*1e-11; //m^3/(kg*s^2)
//die Massen in kg
const double mass_1 = 1.345*1e24;
const double mass_2 = 1.345*1e24;
const double mass_3 = 5.683*1e26;

const double t_max = 5000; //in s
const double delta_t = 0.5; //in s

const int dimension = 12;

typedef struct{
    double x;
    double y;
}s_vec2D;

double vect_betrag2D(double x, double y){
    return sqrt(x*x + y*y);
}

double system_energy(double y[], void* params){
    //Kinetische energie
    double e_kin_m1 = (mass_1/2.0)*(y[6]*y[6] + y[7]*y[7]);
    double e_kin_m2 = (mass_2/2.0)*(y[8]*y[8] + y[9]*y[9]);
    double e_kin_m3 = (mass_3/2.0)*(y[10]*y[10] + y[11]*y[11]);
    double e_kin_tot = e_kin_m1 + e_kin_m2 + e_kin_m3;
    //Potenzielle energie
    //für masse k: pot von masse i + pot von masse j, i,j \neq k
    double e_pot_m1 = (-1)*G*mass_1*mass_2/vect_betrag2D(y[0]-y[2],y[1]-y[3]) +(-1)* G*mass_1*mass_3/vect_betrag2D(y[0]-y[4],y[1]-y[5]);
    double e_pot_m2 = (-1)*G*mass_2*mass_1/vect_betrag2D(y[0]-y[2],y[1]-y[3]) + (-1)*G*mass_2*mass_3/vect_betrag2D(y[2]-y[4],y[3]-y[5]);
    double e_pot_m3 = (-1)*G*mass_3*mass_2/vect_betrag2D(y[4]-y[2],y[5]-y[3]) + (-1)*G*mass_1*mass_3/vect_betrag2D(y[0]-y[4],y[1]-y[5]);
    //gesamt potentialenergie ist 0.5*summe, weil wir jede wirkung doppel zählen: e_pot_12 = e_pot_21
    double e_pot_tot = (e_pot_m1 + e_pot_m2 + e_pot_m3)/2.0;

    return e_kin_tot + e_pot_tot;
}

//Gravitationsbeschleunigung von einer masse auf eine andere
s_vec2D gravity_accel(double target_x, double target_y, double target_m, double actor_x, double actor_y, double actor_m){
    s_vec2D result;
    double dist_actor_target = vect_betrag2D(actor_x - target_x, actor_y - target_y);
    result.x = (-1)*G*actor_m*target_m*(target_x - actor_x)/(dist_actor_target*dist_actor_target*dist_actor_target);
    result.y = (-1)*G*actor_m*target_m*(target_y - actor_y)/(dist_actor_target*dist_actor_target*dist_actor_target);
    return result;
}

int ODE_gravity_2D(double t, const double y[], double dy[], void* params){
    //y = {positionen von allen drei Massen: 1,2,3, geschw. von allen drei massen: 1,2,3}
    
    //Geschwindigkeiten in dy[] schrieben 
    for(int i = 0; i<6; i++){
        dy[i] = y[6+i];
    }

    //Beschl von Masse 1
    s_vec2D grav12 = gravity_accel(y[0],y[1],mass_1,y[2],y[3],mass_2); //in richtung Masse2
    s_vec2D grav13 = gravity_accel(y[0],y[1],mass_1,y[4],y[5],mass_3); //in richtung Masse3
    s_vec2D m1_accel = {m1_accel.x = grav12.x/mass_1 + grav13.x/mass_1, m1_accel.y = grav12.y/mass_1 + grav13.y/mass_1}; //Resultierende Beschl
    //in dy Schrieben
    dy[6] = m1_accel.x;
    dy[7] = m1_accel.y;

    //Beschleunigung von Masse 2
    s_vec2D grav21 = gravity_accel(y[2],y[3],mass_2,y[0],y[1],mass_1); //in richtung Masse2
    s_vec2D grav23 = gravity_accel(y[2],y[3],mass_2,y[4],y[5],mass_3); //in richtung Masse3
    s_vec2D m2_accel = {m1_accel.x = grav21.x/mass_2 + grav23.x/mass_2, m1_accel.y = grav21.y/mass_2 + grav23.y/mass_2}; //Resultierende Beschl
    //in dy schrieben
    dy[8] = m2_accel.x;
    dy[9] = m2_accel.y;

    //Beschleunigung von Masse 3
    s_vec2D grav31 = gravity_accel(y[4],y[5],mass_3,y[0],y[1],mass_1); //in richtung Masse2
    s_vec2D grav32 = gravity_accel(y[4],y[5],mass_3,y[2],y[3],mass_2); //in richtung Masse3
    s_vec2D m3_accel = {m1_accel.x = grav31.x/mass_3 + grav32.x/mass_3, m1_accel.y = grav31.y/mass_3 + grav32.y/mass_3}; //Resultierende Beschl
    //in dy schrieben
    dy[10] = m3_accel.x;
    dy[11] = m3_accel.y;
    return 0;
}

//Anfangswertemethode für reset von Zustandsvektor zwischen Simuls
void init_state(double y[], int dim){
    //pos masse 1 in m
    y[0] = 1.2098*1e7;
    y[1] = 0.0;
    //pos masse 2 in m
    y[2] = 1.2342*1e7;
    y[3] = 0.0;
    //pos masse 3 in m
    y[4] = 0.0;
    y[5] = 0.0;
    //geschw masse 1 in m/s
    y[6] = 0.0;
    y[7] = 6.9262*1e4;
    //geschw masse 2 in m/s
    y[8] = 0.0;
    y[9] = 4.2158*1e4;
    //geschw masse 3 in m/s
    y[10] = 0.0;
    y[11] = 0.0;
}

int main(){
    FILE* euler_data = fopen("euler_data.csv", "w");
    FILE* euler_energy = fopen("euler_energy.csv", "w");
    FILE* RK2_data = fopen("RK2_data.csv", "w");
    FILE* RK2_energy = fopen("RK2_energy.csv", "w");
    FILE* RK4_data = fopen("RK4_data.csv", "w");
    FILE* RK4_energy = fopen("RK4_energy.csv", "w");
    FILE* VV_data = fopen("VV_data.csv", "w");
    FILE* VV_energy = fopen("VV_energy.csv", "w");

    double *y_euler = malloc(12*sizeof(double));
    double *y_rk2 = malloc(12*sizeof(double));
    double *y_rk4 = malloc(12*sizeof(double));
    double *y_vv = malloc(12*sizeof(double));
    double *dy_euler = malloc(12*sizeof(double));
    double *dy_rk2 = malloc(12*sizeof(double));
    double *dy_rk4 = malloc(12*sizeof(double));
    double *dy_vv = malloc(12*sizeof(double));

    /*----------- Simulationen -----------------*/

    //Euler
    init_state(y_euler,dimension);
    init_state(y_rk2,dimension);
    init_state(y_rk4,dimension);
    init_state(y_vv,dimension);

    double t = 0.0;
    while(t < t_max){
        //Euler
        mn_euler_step(t,delta_t,y_euler,dy_euler,ODE_gravity_2D,dimension,NULL);
        double energy_euler = system_energy(y_euler,NULL);
        fprintf(euler_energy, "%.16f, %.16e\n", t, energy_euler);
        fprintf(euler_data, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f\n", y_euler[0], y_euler[1], y_euler[2], y_euler[3], y_euler[4], y_euler[5]); //x,y Komp. von von Masse 1, Masse 2 und Masse 3
        
        //RK2
        mn_rk2_step(t,delta_t,y_rk2,dy_rk2,ODE_gravity_2D,dimension,NULL);
        double energy_rk2 = system_energy(y_rk2,NULL);
        fprintf(RK2_energy, "%.16f, %.16e\n", t, energy_rk2);
        fprintf(RK2_data, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f\n", y_rk2[0], y_rk2[1], y_rk2[2], y_rk2[3], y_rk2[4], y_rk2[5]);
        
        //RK4
        mn_rk4_step(t,delta_t,y_rk4,dy_rk4,ODE_gravity_2D,dimension,NULL);
        double energy_rk4 = system_energy(y_rk4,NULL);
        fprintf(RK4_energy, "%.16f, %.16e\n", t, energy_rk4);
        fprintf(RK4_data, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f\n", y_rk4[0], y_rk4[1], y_rk4[2], y_rk4[3], y_rk4[4], y_rk4[5]);
        
        //VV
        mn_velocity_verlet_step(t,delta_t,y_vv,dy_vv,ODE_gravity_2D,dimension,NULL);
        double energy_vv = system_energy(y_vv,NULL);
        fprintf(VV_energy, "%.16f, %.16e\n", t, energy_vv);
        fprintf(VV_data, "%.16f, %.16f, %.16f, %.16f, %.16f, %.16f\n", y_vv[0], y_vv[1], y_vv[2], y_vv[3], y_vv[4], y_vv[5]);
        
        t += delta_t;
    }

    free(y_euler);
    free(y_rk2);
    free(y_rk4);
    free(y_vv);
    free(dy_euler);
    free(dy_rk2);
    free(dy_rk4);
    free(dy_vv);
    fclose(euler_data);
    fclose(RK2_data);
    fclose(RK4_data);
    fclose(VV_data);
    fclose(euler_energy);
    fclose(RK2_energy);
    fclose(RK4_energy);
    fclose(VV_energy);
    return 0;
}