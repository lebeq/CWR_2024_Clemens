#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../MyNumerics/my_numerics.h"


//Integrand von Gl.2 mit \omega = 2.0
double fourier_integrand_sinc_real(double x, double k){
    if(x==0){
        return cos(k*x);
    }
    else{
        return (sin(2.0*x)/(2.0*x))*cos(k*x);
    }
}

//Integrand von Imaginärteil, Teilaufgabe 5
double fourier_integrand_sinc_imag(double x, double k){
    if(x==0){
        return sin(-k*x);
    }
    else{
        return (sin(2.0*x)/(2.0*x))*sin(-k*x);
    }
}

//Analytische Lösung, weiterhin \omega = 2.0
double analytic_sol(double k){
    return (M_PI/2.0)*mn_rect(k/4.0);
}

//integrator schleifen auswertung
int data_dump(int M_size, double M[], double step, double prec, double integrator(double,double,double,double(double,double),double), double integrand(double,double), FILE* data[]){
    /*
    M ist array mit den Intgrenzen von größe M_size, prec ist delta_x, eig fest = 1e3, data ist array mit pointern zu csv dateien
    step ist Schrittweite
    */
    for(int idx = 0; idx < M_size;idx++){
        double i = -5.0;
        printf("M ist %f, -M ist %f\n, prec is %f", M[idx], -1*M[idx], prec);
        while(i<=5){
            double val = integrator(-1*M[idx],M[idx],prec,integrand,i);
            fprintf(data[idx], "%f, %.16f\n", i, val);
            i += step;
        };
    };
    return EXIT_SUCCESS;
}

int main(){
    FILE* data10 = fopen("fourier_int_10.csv", "w"); //M = 10
    FILE* data30 = fopen("fourier_int_30.csv", "w"); // M = 30
    FILE* data100 = fopen("fourier_int_100.csv", "w"); // M = 100
    FILE* data500 = fopen("fourier_int_500.csv", "w"); // M = 500
    FILE* anadata = fopen("ana_sol.csv", "w"); //analytische Lösung
    FILE* imagdata = fopen("imag_fourier_int_100.csv", "w"); //imaginätreil integral für M = 100

    double delta_k = 0.05; //Schrittweite für Auswertungsschleife - 200 Gleichvrtlte Werte in [-5,5]
    
    double M_arr[] = {10,30,100,500}; //Werte von M
    FILE** outputdata = malloc(sizeof(FILE*) * 4);
    FILE* data[] = {data10,data30,data100,data500};

    //auswertung approx Lsg
    data_dump(4,M_arr,delta_k,1e-3,mn_integrate_simpson_two,fourier_integrand_sinc_real,data);

    double t = -5;
    //auswertun analytische Lsg
    while(t<=5){
        double val = analytic_sol(t);
        fprintf(anadata, "%f, %.16f\n", t, val);
        t += delta_k;
    }

    //zweck einheitlichkeit imaginärteil auch mit data_dump ausgewertet
    double M_imag_arr[] = {100};
    FILE* imagdata_arr[] = {imagdata};
    data_dump(1,M_imag_arr,delta_k,1e-3,mn_integrate_simpson_two,fourier_integrand_sinc_imag,imagdata_arr);
    /*double r = -5.0;
    //auswertung Imaginärteil für M=100
    while(r<=5){
        double val = mn_integrate_simpson_two(-100,100,1e-3,fourier_integrand_sinc_imag,r);
        fprintf(imagdata, "%f, %.16f\n", r, val);
        r += delta_k;
    }*/
    free(outputdata);
    fclose(data10);
    fclose(data30);
    fclose(data100);
    fclose(data500);
    fclose(anadata);
    fclose(imagdata);
    return EXIT_SUCCESS;
}