#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>

double trapezregel(double (*f)(double), double a, double b, double h) {
    double sum = 0;
    sum += f(a) / 2;
    for (double i = a + h; i < (b - h); i += h) {
        sum += f(i);
    }
    sum += f(b) / 2;
    return sum * h;
}

double mittelpunktsregel(double (*f)(double), double a, double b, double h) {
    double sum = 0;
    for (double i = (a + h / 2); i < (b - h / 2); i += h) {
        sum += f(i);
    }
    return sum * h;
}

double simpsonregel(double (*f)(double), double a, double b, double h) {
    // check if N is even
    int N = int((b - a) / h);
    if (N % 2 == 1) {
        std::cerr << "N is not even! This can cause problems." << std::endl;
    }

    double sum = 0;
    sum += f(a) + f(b);
    for (double i = a + h; i < (b - h); i += h) {
        sum += 2 * f(i);
    }
    for (double i = a + h / 2; i < (b - h / 2); i += h) {
        sum += 4 * f(i);
    }
    return sum * h / 6;
}

#endif
