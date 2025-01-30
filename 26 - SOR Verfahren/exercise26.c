#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include "../MyNumerics/my_numerics.h"


typedef enum{
    INSIDE,
    DIRICHLET,
    NEUMANN,
    DISABLED
} Node_Type;

typedef struct{
    Node_Type type;
    double val;
} node;

double dist(double x, double y){
    return sqrt(x*x + y*y);
}

void InitDomain(int N, node F[][N], double N_val, double E_val, double S_val, double W_val){
    //the double params are the cardinal directions, not the same N as the lattice point number!!
    //boundary gets type DIRICHLET and according _val
    //points on inside get type INSIDE and val = 0.0
    for(int i = 0; i < N; i++){
        for(int j = 0; j<N; j++){
            //bottom boundary
            if(i == 0){
                F[i][j].type = DIRICHLET;
                F[i][j].val = S_val;
            }
            //top boundary
            else if(i == N-1){
                F[i][j].type = DIRICHLET;
                F[i][j].val = N_val;
            }
            //left boundary
            else if(j == 0){
                F[i][j].type = DIRICHLET;
                F[i][j].val = W_val;
            }
            //right boundary
            else if(j == N-1){
                F[i][j].type = DIRICHLET;
                F[i][j].val = E_val;
            }
            //points inside
            else{
                F[i][j].type = INSIDE;
                F[i][j].val = 0.0;
            }

        }
    }
}

void InsertCircle(int N, double cx, double cy, double R, node F[][N], double val){
    //inserts circle/hole in potential grid
    //whole disc gets DIRICHLET and val
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            if(dist(cy - i, cx - j) < R){
                F[i][j].type = DIRICHLET;
                F[i][j].val = val;
            }
        }
    }
}

int PoissonGaussSeidel(int N, node F[][N], int iter_max, double tolerance, double alpha){
    //alpha is the relacation parameter
    //to calc. the k+1 dstep values we need the k-th step values
    node (*F_old)[N] = malloc(sizeof(node[N][N]));
    //variable for quadratic relative change
    double rel_chng_sq = 100.0;
    int iteration = 0;
    while(iteration < iter_max && sqrt(rel_chng_sq) > tolerance){
        for(int i = 0; i<N;i++){
            for(int j = 0; j<N; j++){
                //save k-th step values in f_old
                F_old[i][j].type = F[i][j].type;
                F_old[i][j].val = F[i][j].val;
            }
        }
        rel_chng_sq = 0.0;
        for(int i = 0; i<N; i++){
            for(int j = 0; j<N; j++){
                switch(F[i][j].type){
                    case INSIDE:
                        //inner lattice point
                        //North + South + West + East
                        //MIT RELAXATION PARAMETER ALPHA
                        F[i][j].val = (1.0 - alpha)*F_old[i][j].val + alpha*(F[i-1][j].val + F_old[i+1][j].val + F[i][j-1].val + F_old[i][j+1].val)/4.0;
                        break;
                    case DIRICHLET:
                        break;
                    default:
                        break;
                }
                if(F_old[i][j].val == 0.0){
                    rel_chng_sq += (F[i][j].val)*(F[i][j].val);
                }
                else{
                    rel_chng_sq += ((F[i][j].val - F_old[i][j].val)/F_old[i][j].val)*((F[i][j].val - F_old[i][j].val)/F_old[i][j].val);
                }
            }
        }
        iteration ++;
    }
    free(F_old);
    return iteration;
}

int main(){
    int N = 256; //Anzalh gitterpunkte in eine Richtung -> N^2 Gitterpunkte insg.
    node (*F)[N] = malloc(sizeof(node[N][N]));
    FILE* potential_data = fopen("potential_data.csv", "w");
    FILE* iteration_data = fopen("iteration_data.csv", "w");
    double delta_x = 1.0/N; // in m
    double r_electrode = 0.1; // in m
    double wand = 0.0; //Spannung der äußeren Wände

    //Iterationsparamter
    double max_iter = 1e4;
    double err_tol = 1e-3;
    
    //Spannung der Elektorden
    double electrode_left = -1.0;
    double electorde_right = 1.0;  

    //F initialisieren
    InitDomain(N,F,wand,wand,wand,wand);

    //Elektroden reinsätzen
    InsertCircle(N,0.25/delta_x,0.5/delta_x,r_electrode/delta_x,F,electrode_left);
    InsertCircle(N,0.75/delta_x,0.5/delta_x,r_electrode/delta_x,F,electorde_right);

    //Laplacegleichung Lösen
    double alpha = 1.0;
    double delta_alpha = 0.01;
    //int iter = PoissonGaussSeidel(N,F,max_iter,err_tol,alpha);

    while(alpha < 2){
        //F initialisieren
        InitDomain(N,F,wand,wand,wand,wand);

        //Elektroden reinsätzen
        InsertCircle(N,0.25/delta_x,0.5/delta_x,r_electrode/delta_x,F,electrode_left);
        InsertCircle(N,0.75/delta_x,0.5/delta_x,r_electrode/delta_x,F,electorde_right);

        int iter = PoissonGaussSeidel(N,F,max_iter,err_tol,alpha);
        fprintf(iteration_data, "%.16f, %d\n", alpha, iter);
        alpha += delta_alpha;
    }

    //CSV-Datei beschriften
    for(int i = 0; i<N; i++){
        fprintf(potential_data, "%.16f", F[i][0].val);
        for(int j = 1; j<N; j++){
            //printf("check, value of F[%d][%d] is %.16f\n",i,j,F[i][j].val);
            fprintf(potential_data, ", %.16f", F[i][j].val);
        }
        fprintf(potential_data, "\n");
    }

    fclose(potential_data);
    fclose(iteration_data);
    free(F);
    return 0;
}