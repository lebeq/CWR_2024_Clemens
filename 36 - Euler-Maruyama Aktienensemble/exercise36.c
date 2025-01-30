#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"

typedef struct{
    Tuple constants;
    int dimension;
    gsl_rng* generator;
} s_params;

/*------------------------------ Simulationsparameter -------------------------*/
double rendite = 0.1;
double voli = 0.1;
double K_0 = 10.0;


/*------------------------ Die SDE ---------------------------------------*/
int stockSDE(double t, const double y[], double f[], double g[], void* params){
    s_params *p = (s_params*)params;
    double mu = p->constants.x;
    double sigma = p->constants.y;
    int dim = p->dimension;
    for(int i = 0; i<dim; i++){
        f[i] = mu*y[i];
        g[i] = sigma*y[i];
    }
    return 0;
}


int main(){
    double t = 0.0;
    double dt = 1e-3; //in Jahren
    double t_max = 20.0; //in jahren
    FILE* stonks = fopen("stock_data.csv", "w");
    double *y = malloc(200*sizeof(double));
    gsl_rng* gen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gen, time(NULL));
    //Anfangswerte
    for(size_t i = 0; i<200; i++){
        y[i] = K_0;
    }
    Tuple consts = {rendite, voli};
    s_params params = {consts, 200, gen};

    while(t<t_max){
        
        mn_euler_maruyama_step(t,dt,y,stockSDE,params.dimension,gen,&params);
        
        //CSV beschriften
        fprintf(stonks, "%f",t);
        for(int l = 0; l<params.dimension; l++){
            fprintf(stonks, ", %.16e", y[l]);
        }
        fprintf(stonks, "\n");
        t+=dt;
    }

    free(y);
    fclose(stonks);
    gsl_rng_free(gen);
    return 0;   
}