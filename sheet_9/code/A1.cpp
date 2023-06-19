#include <iostream>
#include "optimize.h"

double f(double x){
    return x*x - 2;
}
int main(){
    std::cout << "\nExercise 1 >>>" << std::endl;
    double eps = 1E-9;

    // a) interval biscetion method
    std::cout << "a) Bisection method\n";
    double x0 = -0.5;
    double y0 = -0.1;
    double z0 = 2.0;

    Bisection bisection(f, x0, y0, z0, eps);
    double min_bisection = bisection.min("build/A1_bisection.csv");
    std::cout << "Minimum at x0 = " << min_bisection << "\n\n";
    
    // b) Newtons method
    std::cout << "b) Newtons method\n";
    double x00 = 1.0;

    Newton newton(f, x00, eps);
    double min_newton = newton.min("build/A1_newton.csv");
    std::cout << "Minimum at x0 = " << min_newton << "\n\n";

    return 0;
}