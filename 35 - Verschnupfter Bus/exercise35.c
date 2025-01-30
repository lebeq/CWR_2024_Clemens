#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>

#include  "../MyNumerics/my_numerics.h"

typedef enum{
    INFECTED,
    HEALTHY
} Node_Type;

typedef struct{
    Node_Type type;
} node;

/*-------------------------- GSL Generator setup ----------------------------*/
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

/*------------------------ Simulationsparameter --------------------------*/
int N = 40; //Sitzplätze
double p = 0.25; //Wahrscheinlichkeit der Ansteckung in Prozent
double initial_infected = 0.2; //Prozent der initial angesteckten


/*---------------------------- Iinit Domain -----------------------------*/
void InitDomain(int N, node y[N]){
    for(size_t i = 0; i<N; i++){
        y[i].type = HEALTHY;
    }
}

void initBus(int N, node y[N], double init_infect){
    int counter = 0;
    gsl_rng* r = get_gsl_gen();
    while(counter < init_infect*N){
        int seat = random_uniform(r,0,N);
        if(y[seat].type != INFECTED){
            y[seat].type = INFECTED;
            counter ++;
        }
    }
    free_gsl(r);
}

/*void simul(int N, node y[N], double probability, double dt, int t_max){
    initBus(N,y,initial_infected);
    int t = 0;
    while(t < t_max){

    }
}*/

int main(){
    FILE* infect_data = fopen("infect_data.csv", "w");
    node *y = malloc(sizeof(node[N]));
    InitDomain(N,y);
    initBus(N,y,initial_infected);
    int t = 0;
    int t_max = 10000;
    int infected_counter = initial_infected*N;
    gsl_rng* gen = get_gsl_gen();
    while(t < t_max){
        int seat = random_uniform(gen, 0, N);
        if(y[seat].type == INFECTED){
            //mit Wkeit p genesung
            double rnd = gsl_rng_uniform(gen);
            if(rnd < p){
                y[seat].type = HEALTHY;
            }
            //mit Wkeit 1-p nachbarn anstecken
            else{
                //Sitze 0 und N-1 stecken nur eine Person an
                if(seat == 0){
                    y[seat + 1].type = INFECTED;
                }
                else if(seat == N - 1){
                    y[seat - 1].type = INFECTED;
                }
                else{
                    y[seat - 1].type = INFECTED;
                    y[seat + 1].type = INFECTED;
                }
            }
        }
        //Zähle alle infizierten zum Zeitpunkt t
        infected_counter = 0;
        for(size_t i = 0; i<N; i++){
            if(y[i].type == INFECTED){
                infected_counter ++;
            }
        }
        fprintf(infect_data, "%d, %d\n", t, infected_counter);
        t ++;
    }
    fclose(infect_data);
    free_gsl(gen);
    free(y);
    return 0;
}