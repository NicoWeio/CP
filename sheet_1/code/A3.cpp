#include "integrate.h"
#include <iostream>
#include <fstream>
#include <math.h>


double f1(double x){
    return exp(-x)/x;
}

double f2(double x){
    return x*sin(1/x);
}

// Ich habe es leider nicht geschafft die drei unten stehenden Funktionen zu einer zusammenzufassen mit einem
// zusätzlichen Argument: Zeiger auf Funktion (trapezregel, mittelpunkt, simpson)
// da diese funktionen selbst bereits als argument einen zeiger auf eine funktion besitzen

void integrate_trapez(double (*f)(double), double a, double b, std::string path_output, double eps){
    std::ofstream myfile(path_output);

    double err = 1.0;
    double h = abs(b-a)/2;

    double I;
    double I_old;
    
    I = trapezregel(f1, a, b, h);
    while(err>eps){
        myfile << I << std::endl;
        h = h/2;
        I_old = I;
        I = trapezregel(f1, a, b, h);

        err = abs((I-I_old)/I_old);
    }
    myfile.close();
}

void integrate_mittelpunkt(double (*f)(double), double a, double b, std::string path_output, double eps){
    std::ofstream myfile(path_output);

    double err = 1.0;
    double h = abs(b-a)/2;


    double I;
    double I_old;
    
    I = mittelpunktsregel(f1, a, b, h);
    while(err>eps){
        myfile << I << std::endl;
        h = h/2;
        I_old = I;
        I = mittelpunktsregel(f1, a, b, h);

        err = abs((I-I_old)/I_old);
    }
    myfile.close();
}

void integrate_simpson(double (*f)(double), double a, double b, std::string path_output, double eps){
    std::ofstream myfile(path_output);

    double err = 1.0;
    double h = abs(b-a)/2;


    double I;
    double I_old;
    
    I = simpsonregel(f1, a, b, h);
    while(err>eps){
        myfile << I << std::endl;
        h = h/2;
        I_old = I;
        I = simpsonregel(f1, a, b, h);

        err = abs((I-I_old)/I_old);
    }
    myfile.close();
}


int main() {
    const double eps = pow(10,-4);
    const double I1_exact = 0.219384;
    const double I2_exact =  0.378530;

    // a)
    integrate_trapez(f1, 1, 100, "build/I1_trapez.csv", eps);
    integrate_mittelpunkt(f1, 1, 100, "build/I1_mittelpunkt.csv", eps);
    integrate_simpson(f1, 1, 100, "build/I1_simpson.csv", eps);


    // b)
    //integrate_trapez(f2, 0, 1, "build/I2_trapez.csv", eps);
    //integrate_mittelpunkt(f2, 0, 1, "build/I2_mittelpunkt.csv", eps);
    //integrate_simpson(f2, 0, 1, "build/I2_simpson.csv", eps);

    return 0;
}