#include <iostream>
#include <fstream>
#include <functional>
#include <math.h>

double f2(double x) {
    // TODO: Richtige Funktion(en) verwenden
    return x*sin(1/x);
}

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

int main() {
    std::cout << trapezregel(f2, 0.1, 1, 0.00000001) << std::endl;
    std::cout << mittelpunktsregel(f2, 0.1, 1, 0.00000001) << std::endl;
    std::cout << simpsonregel(f2, 0.1, 1, 0.00000001) << std::endl;
    return 0;
}
