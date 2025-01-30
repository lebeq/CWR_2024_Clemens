#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include "../MyNumerics/my_numerics.h"

/*----------------------- Parameter -----------------------------*/
const int dimension = 5;

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

/*------------------ auxilliary functions ------------------------------*/

//konstante dichte = 1
double density_func(int D, double *x){
    if(mn_euclid_norm(D,x) <= 1.0){
        return 1.0;
    }
    else{
        return 0.0;
    }
}

void init_hypercube(int D, double *R){
    //Grenzen von Würfel
    for(int i = 0; i<2*D - 1; i+=2){
        R[i] = -1.0;
        R[i+1] = 1.0;
    }
}
/*
NOTE: anscheinend ist es genug die gamma funktion mit max_iter = 10 laufen zu lassen.
Ab 1e2 dauert es viel zu lange und die werte von 1e1 sind schon sehr genau.
Wahrscheinlisch ist meine gamma funktion overkill und einfach ne simpson integration 
mit großer oberer Grenze (1e3 sogar schon) genügt vollkommen.
*/
double volume_analytical(int D, double r){
    double factor = pow(r,D)*pow(M_PI,D/2.0);
    return factor/mn_gamma_function_real(D/2.0 + 1,10,1e-2);
}

double domain_volume(int D, double *R){
    //Volume of hypercube with sidelenght L = 2
    double side_length = fabs(R[0]-R[1]);
    return pow(side_length,D);
}

/*--------------------- MC Integration ----------------------------*/

double mc_integrate(gsl_rng* generator, int D, double *R, double integrand(int, double *), int N){
    //Vorfaktor von integral
    double faktor = domain_volume(D,R)/N;
    double result = 0.0;
    double *x_rnd = malloc(D*sizeof(double));
    //Integral berechnen
    //N stützsdtellen aufaddieren
    for(int l = 0; l<N; l++){
    //generate random Vector
        for(int i = 0; i<D; i++){
            x_rnd[i] = random_uniform(generator, -1.0, 1.0);
            //printf("mc_integrate rand vect loop, x_rnd[%d] is %.16f\n", i, x_rnd[i]);
        }
        //the integrand is only 1 if we are inside the sphere
        result += integrand(D,x_rnd);
    }
    result = (faktor)*result;
    
    free(x_rnd);
    return result;
}

int main(){
    /*-------------------- TEIL 4 -------------------------*/
    //Dimension etc
    int D[7] = {1,2,3,4,5,6,7};
    int N = 100000;
    int M = 100;
    FILE* mc_data = fopen("mc_data.csv", "w");
    FILE* analytical_data = fopen("ana_data.csv", "w");
    /*
    for(int j = 0; j<7; j++){
        int dim = D[j];
        double *R = malloc(2*dim*sizeof(double));
        gsl_rng* r = get_gsl_gen();
        //Grenzen von Würfel
        init_hypercube(dim,R);

        for(int t = 0; t<M; t++){
            double *I_k = malloc(M*sizeof(double));
            double mc_output = 0.0;
            //die I_k in tegrale berechnen und in arreay speichern für variance
            for(int i = 0; i<M; i++){
                I_k[i] = mc_integrate(r,dim,R,density_func,N/M);
                mc_output += I_k[i];
            }
            //Das Mittel der I_k berechenen
            double I = (1.0/M)*mc_output;

            //Variance
            double var = 0.0;
            for(int i = 0; i<M; i++){
                var += (I_k[i] - I)*(I_k[i] - I)/((M-1)*(M-1));
            }
            var = sqrt(var);
            //integration in CSV eintragen: dim, I, variation
            fprintf(mc_data, "%d, %.16f, %16f\n",dim,I,var);
            //zugehörige variance
            //fprintf(mc_data, ", %.16f", var);
            free(I_k);
        }
        //analytische LSg in datei schreiben
        fprintf(analytical_data, "%.d, %.16f\n", dim, volume_analytical(dim,1.0));
        free(R);
        free_gsl(r);
    }
    */

    
    fclose(analytical_data);
    fclose(mc_data);

    printf("############### PART 4 DONE ##############\n");

    /*----------------------------- TEIL % -------------------------*/
    FILE* hypersphere_data = fopen("hypersphere_data.csv", "w");
    FILE* analytical_hypersphere = fopen("analytical_hypersphere.csv", "w");

    double expo_s = 2.0;
    double expo_e = 8.0;
    double delta_expo = (expo_e - expo_s)/100.0;

    //analytische LSG
    double ana_sol = volume_analytical(dimension,1.0);

    for(int i = 0; i<100; i++){
        double expo = expo_s + i*delta_expo;
        printf("expo is %f\n", expo);
        double N = pow(10,expo);
        printf("N is %f\n",N);
        double *R = malloc(2*dimension*sizeof(double));
        //Hypercube Grenzen setzen
        init_hypercube(dimension,R);
        //rng holen
        gsl_rng *r = get_gsl_gen();

        //Numerisch das Volumen rechnen
        for(int j = 0; j<M; j++){
            double *J_k = malloc(M*sizeof(double));
            double sum = 0.0;
            //100x kleniere integrationen
            for(int t = 0; t<M; t++){
                J_k[t] = mc_integrate(r,dimension,R,density_func,N/M);
                //printf("J_K[%d] is %f\n",t,J_k[t]);
                sum += J_k[t];
                //printf("sum is %f\n", sum);
            }
            //Mittel
            double J = (1.0/M)*sum;

            //Variance
            double var = 0.0;
            for(int l = 0; l<M; l++){
                var += (J_k[i] - J)*(J_k[i] - J)/((M-1)*(M-1));
            }
            var = sqrt(var);
            //Abweichung zu analytischen LSG
            double diff_num_ana = fabs(J - ana_sol);
            //Ergebnisse in CSV schreiben: N, J, var
            fprintf(hypersphere_data, "%.16f, %.16f, %.16f, %.16f\n", N,J,var,diff_num_ana);
            
            free(J_k);
        }
        printf("### expo %f done.\n", expo);
        //analytische LSG in CSV schreiben
        fprintf(analytical_hypersphere, "%.16f, %.16f\n", N, ana_sol);

        free(R);
        free_gsl(r);
    }

    fclose(analytical_hypersphere);
    fclose(hypersphere_data);
    return 0;
}