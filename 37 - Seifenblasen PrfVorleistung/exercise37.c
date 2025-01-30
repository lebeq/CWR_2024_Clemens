#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#include "../MyNumerics/my_numerics.h"

/*typedef struct{
    double x;
    double y;
    double z;
} point_3D;

double euclidean_norm_3D(point_3D P){
    return sqrt(P.x*P.x + P.y*P.y + P.z*P.z);
}*/
/*------------------------ gsl_rng setup and co --------------------------*/
gsl_rng* get_gsl_gen(){
    gsl_rng* generator = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(generator, time(NULL));
    return generator;
}

void free_gsl(gsl_rng* generator){
    if(generator != NULL){
        gsl_rng_free(generator);
    }
}

double random_uniform(gsl_rng* generator, double min, double max){
    double extent = max - min;
    double r = gsl_rng_uniform(generator)*extent;
    return r + min;
}



/*------------------ Density funcs ------------------------*/
double point_in_sphere(double *p, double *c){
    double diff[3] = {0.0, 0.0, 0.0};
    for(int j = 0; j<3; j++){
        diff[j] = c[j] - p[j];
    }
    if(mn_euclid_norm(3,diff) <= 1.0){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

double point_in_sphere_union(double *c1, double *c2, double *p, double func(double *, double *)){
    if(func(p,c1) == 1.0 || func(p,c2) == 1.0){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

/*------------------- MC Domain setup -------------------*/

void init_domain(double *c1, double *c2, double *R){
    R[0] = fmin(c1[0],c2[0]) - 1.0;
    R[1] = fmax(c1[0],c2[0]) + 1.0;
    for(int i = 2; i<6; i+=2){
        R[i] = -1.0;
        R[i+1] = 1.0;
    }
}

double domain_volume(double *R){
    return (R[1] - R[0])*(R[3] - R[2])*(R[5] - R[4]);
}

/*-------------- Analytic solution ---------------------- */
double ana_volume(double d){
    //für radius = 1
    if(fabs(d) <= 2.0){
        return (4.0/3.0)*M_PI + M_PI*fabs(d)*(1.0 - (d*d)/12.0);
    }
    else if(fabs(d) > 2.0){
        return 8.0*M_PI/3;
    }
}
/*--------------------- Monte Carlo integration -------------*/
double mc_integrate(gsl_rng *generator, double *R, double *c1, double *c2, double func(double *, double *), double integrand(double *, double *, double *, double *), int N){
    double factor = domain_volume(R)/N;
    double result = 0.0;
    double *x_rnd = malloc(3*sizeof(double));
    //MC integration
    for(int l = 0; l<N; l++){
        //rand vect innerhalb domain erstellen
        for(int j = 0; j<3; j++){
            x_rnd[j] = random_uniform(generator,R[2*j], R[2*j+1]);
        }
        result += integrand(c1,c2,x_rnd,func);
    }
    result = factor*result;

    free(x_rnd);
    return result;
}

int main(){
    FILE* vol_data = fopen("volume_data.csv", "w");
    FILE* mean_data = fopen("mean_data.csv", "w");
    FILE* ana_data = fopen("analytic_data.csv", "w");
    //Die Mittelpunkte der beiden Blasen
    double c1[3] = {0.0,0.0,0.0};
    double c2[3] = {0.0,0.0,0.0};


    //Domain init
    double R[6];
    for(int i = 0; i<6; i++){
        R[i] = 0.0;
    }
    init_domain(c1,c2,R);
    printf("domain volume is %f\n", domain_volume(R));

    gsl_rng* r = get_gsl_gen();

    double d = -2.5;
    double delta_d = 0.05;
    int N = 1e5;
    int M = 100;
    
    while(d <= 2.5){
        double mean = 0.0;

        //update center of second sphere
        c2[0] = d;

        //init domain
        init_domain(c1,c2,R);

        for(int i = 0; i<500; i++){
            double *I_k = malloc(M*sizeof(double));
            double res = 0.0;
            //I_k berechnen
            for(int j = 0; j<M; j++){
                I_k[j] = mc_integrate(r,R,c1,c2,point_in_sphere,point_in_sphere_union,N/M);
                res += I_k[j];
            }
            //Mittel der I_k ist das integral
            double I = (1.0/M)*res;

            //mean aktualisieren
            mean += I;

            //Variance
            double var = 0.0;
            for(int l = 0; l<M; l++){
                var += (I_k[l] - I)*(I_k[l] - I)/((M-1)*(M-1));
            }
            var = sqrt(var);
            
            //CSV Beschriften: d, I, var
            fprintf(vol_data, "%f, %.16f, %.16f\n", d,I,var);
        }
        
        //Mittelwert der ergebnisse für geg. d
        mean = (1.0/500.0)*mean;
        fprintf(mean_data, "%f, %.16f\n", d, mean);
        
        //analytische LSG
        double an_vol = ana_volume(d);
        fprintf(ana_data, "%f, %.16f\n", d, an_vol);
        
        d += delta_d;
    }


    free_gsl(r);
    fclose(mean_data);
    fclose(ana_data);
    fclose(vol_data);
    return 0;
}