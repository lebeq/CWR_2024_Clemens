#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_numerics.h"


double f(double x){
    return x*x*(x-1);
}



int main(){
    FILE* data = fopen("diff_num.csv", "w");
    double wert = f(2.0);
    double ablt = mn_diff(1.0,1e-5,f);
    printf("Wert von f an 2: %f\n", wert);
    printf("Ableitung numerisch an 1: %f\n",ablt);

    int expo_s = 0;
    int expo_e = -16;
    int N = 100;
    double expo_delta = (expo_e - expo_s)/N;
    for(int i = 0; i<100; i++){
        double delta = pow(10,expo_s + i*expo_delta);
        double err = fabs(ablt-mn_diff(1.0,delta,f));
        fprintf(data, "%.16e, %.16e\n",delta, err);
    }

    return EXIT_SUCCESS;
}