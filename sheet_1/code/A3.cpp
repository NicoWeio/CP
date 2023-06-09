#include "integrate.h"
#include <fstream>
#include <iostream>
#include <math.h>

double f1(double x) {
    return exp(-x) / x;
}

double f2(double x) {
    if (x == 0) {
        return 0;
    }
    return x * sin(1 / x);
}

typedef double (*integrator_fn)(double (*f)(double), double a, double b, double h);

void integrate_general(integrator_fn integrator, double (*f)(double), double a, double b, std::string path_output, double eps) {
    std::ofstream myfile(path_output);

    double err = 1.0;
    double h = abs(b - a) / 2; // start value for h(N=1): h=(b-a)/N

    double I;
    double I_old;

    I = integrator(f1, a, b, h);
    while (err > eps) {
        myfile << h << ", " << I << std::endl; // write h and result of integral I to csv file for each iteration
        h = h / 2;
        I_old = I;
        I = integrator(f1, a, b, h);

        err = abs((I - I_old) / I_old);
    }
    myfile.close();
}

int main() {
    const double eps = pow(10, -4);
    const double I1_exact = 0.219384;
    const double I2_exact = 0.378530;

    // a)
    integrate_general(trapezregel, f1, 1, 100, "build/I1_trapez.csv", eps);
    integrate_general(mittelpunktsregel, f1, 1, 100, "build/I1_mittelpunkt.csv", eps);
    integrate_general(simpsonregel, f1, 1, 100, "build/I1_simpson.csv", eps);

    // b)
    // integrate_general(trapezregel, f2, 0, 1, "build/I2_trapez.csv", eps);
    // integrate_general(mittelpunktsregel, f2, 0, 1, "build/I2_mittelpunkt.csv", eps);
    // integrate_general(simpsonregel, f2, 0, 1, "build/I2_simpson.csv", eps);

    return 0;
}
