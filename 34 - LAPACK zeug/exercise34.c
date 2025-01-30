#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <math.h>
#include <time.h>

#include "../MyNumerics/my_numerics.h"

double solve_linear_system(lapack_int n){
    FILE* file = fopen("matrices.dat", "r");
    /*------- LAPACKE  Params ------*/
    lapack_int info;
    lapack_int m = n;
    lapack_int lda = m;
    lapack_int ldb = 1;
    lapack_int nrhs = 1;
    /*------- MAtrix und Vektor erstellen -----*/
    FILE *file_A = fopen("A.dat", "r");
    FILE *file_b = fopen("b.dat", "r");
    double *A = (double *)malloc(m*m*sizeof(double));
    double *b = (double *)malloc(m*sizeof(double));

    /*-------- array das änderungen beim pivoting speichert ----------*/
    lapack_int *ipiv = (lapack_int *)malloc(m*sizeof(lapack_int));

    /*-------- Matrix und Vektor beschrieben -----------*/
    for(int i = 0; i<m; i++){
        for(int j = 0; j<m; j++){
            fscanf(file, "%lf", &A[i*m+j]);
        }
        fscanf(file, "%lf", &b[i]);
    }

    clock_t time_1 = clock();

    info = LAPACKE_dsysv(LAPACK_ROW_MAJOR, 'U', m, nrhs, A, lda, ipiv, b, ldb);

    clock_t time_2 = clock();

    double exec_time = (double)(time_2 - time_1)/CLOCKS_PER_SEC;
    //printf("Time for solving %dx%d is %.16f seconds\n", m, m, exec_time);

    fclose(file);
    free(ipiv);
    free(A);
    free(b);
    fclose(file_A);
    fclose(file_b);
    return exec_time;
}


int main(){
    FILE* timing = fopen("time_data.csv", "w");

    lapack_int n = 500;
    lapack_int delta_n = 500;
    lapack_int n_max = 5000;

    while(n <= n_max){
        double time = solve_linear_system(n);
        fprintf(timing, "%d, %.16f\n", n, time);
        n += delta_n;
    }


    fclose(timing);
    return 0;
}