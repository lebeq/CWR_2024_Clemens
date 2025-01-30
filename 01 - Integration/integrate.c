#include <stdio.h>
#include <stdlib.h>
#include <math.h>


//berechnet integranden
double p1(double x){
    return x*x*x - x/2;
}

//Integration mit Rechtecken
double integrate_midpoint(double left, double right, int N, double integrand(double y)){
    double delta_y = (left - right)/N;
    double sum = 0.0;

    for(int k = 0; k < N; k++){
        double y = left + (k-0.5)*delta_y;
        double f = integrand(y);
        double A = f*delta_y;
        sum += A;
    }
    return sum;
}

//Integration mit Trapezen
double integrate_trap(double left, double right, int N, double integrand(double y)){
    double delta_y = (left - right)/N;
    double sum = 0.0;
    for(int i = 0; i < N; i++){
        double y1 = left + i*delta_y;
        double y2 = left + (i+1)*delta_y;
        double f1 = integrand(y1);
        double f2 = integrand(y2);
        double A = (f1+f2)*delta_y/2;
        sum += A;
    }
    return sum;
}

int err_plot(int max_exp, double left, double right, double integrand(double y), double integrator(double left, double right, int N, double integrand(double y)), FILE* data){
    int N = 10;
    for(int i = 1; i < max_exp; i++){
        int exp = (int) pow(N,i); 
        double num_int = integrator(left,right,exp,integrand);
        double exct_int = 3.0;
        double abs_diff = fabs(num_int - exct_int);
        fprintf(data, "%i, %.16f\n", exp, abs_diff);
    }
    return 0;
}

int main(){
    double left = 0.0; //untere integrationsgrenze
    double right = 2.0; //obere int-grenze
    int N = 10000; //intervallschritt
    FILE* midpt_data = fopen("data_midpoint.csv", "w");
    FILE* trpz_data = fopen("data_trap.csv", "w");

    //berechnung von integral
    double I_midpoint = integrate_midpoint(left, right, N, p1);
    double I_trapez = integrate_trap(left,right,N,p1);

    //Fehler in CSV speichern
    err_plot(10,left,right,p1,integrate_midpoint, midpt_data);
    err_plot(10,left,right,p1,integrate_trap,trpz_data);

    //ausgabe von integralwert
    printf("Integralwert midpoint: %.16f\n", I_midpoint);
    printf("Integralwert trapeze: %.16f\n", I_trapez);

    fclose(midpt_data);
    fclose(trpz_data);

    return EXIT_SUCCESS;
}