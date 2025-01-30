#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>

#include "../MyNumerics/my_numerics.h"

/*---------------- Parameter ----------------*/
const double L = 10.0;
const double N = 10.0;
const double delta_x = L/N;
const double delta_x_sq = delta_x*delta_x;

/*--------------- Potenziale ---------------*/
double null_pot(){
    return 0.0;
}

double quad_pot(double x){
    //Fedekonstante k=1
    return -0.5*x*x;
}

/*------------------ Matrix Print func --------------*/
void print_matrix_rowmajor(char *desc, lapack_int m, double *mat, lapack_int ldm){
    lapack_int i,j;
    printf("\n %s \n", desc);
    for(i=0; i<ldm; i++){
        for(j=0; j<m; j++){
            printf(" %6.2f", mat[i*ldm+j]);
        }
        printf("\n");
    }
}

/*------------- Hamiltonian -------------*/

void init_hamiltonian(double *H, int N, double dx, double V(double)){
    //initialisiert die NxN Hamiltonian Matrix
    for(int i = 0; i<N; i++){
        for(int j = 0; j<N; j++){
            if(j==i){
                H[i*N+j] = V(i) + delta_x_sq;
            }
            else if(j == i+1 || j== i-1){
                H[i*N+j] = -0.5*delta_x_sq;
            }
            else{
                H[i*N+j] = 0.0;
            }
        }
    }
}

/*--------------- LAPACKE stuff ---------------------------*/
lapack_int lapacke_diagonalize(double *A, double *W, int N){
    double wkopt;
    double *work;
    int n = N, lda = N, info, lwork;
    /*--------- query and allocate the optimal workspace ------*/
    info = LAPACKE_dsyev( LAPACK_COL_MAJOR,'V', 'U', n, A, lda, W);
    if(info > 0){
        printf("Failed to compute iegnevals.\n");
        exit(1);
    }
    /*------- Print eigenvals and eigenvects in matrix -------------*/
    char *EV = "Eigenvectors in columns";
    char *EW = "Eigenvalues on diagonal";
    print_matrix_rowmajor(EW ,n,W,1);
    print_matrix_rowmajor(EV, n, A, lda);

    /*------ Free workspace ------*/
    free((void*)work);
    return 0;
}

int main(){
    double *H = malloc(N*N*sizeof(double));
    double *W = malloc(N*sizeof(double));
    init_hamiltonian(H,N,delta_x,null_pot);
    char *descrp = " Hamiltonian mit Potenzial = 0.0";
    print_matrix_rowmajor(descrp,N,H,N);
    lapacke_diagonalize(H,W,N);
    free(H);
    free(W);
    return 0;
}