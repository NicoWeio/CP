#include <iostream>
#include <math.h>
using namespace std;

// █ numerically unstable functions >>>>
// x ≫ 1
double f1(double x) {
    return 1.0 / sqrt(x) - 1.0 / sqrt(x + 1);
}

// x ≪ 1
double f2(double x) {
    return (1.0 - cos(x)) / sin(x);
}

// delta ≪ 1
double f3(double x, double delta) {
    return sin(x + delta) - sin(x);
}
// <<<<

// █ numerically stable functions >>>>
//  x ≫ 1
double g1(double x) {
    return (sqrt(x + 1) - sqrt(x)) / (sqrt(x + 1) * sqrt(x));
}

// x ≪ 1
double g2(double x) {
    return 1 / sin(x) - cos(x) / sin(x);
}

// delta ≪ 1
double g3(double x, double delta) {
    return sin(x) * cos(delta) + cos(x) * sin(delta) - sin(x);
}
// <<<<

// function error calculation
double rel_error(double a, double b) {
    return (a - b) / a;
}

int main() {
    // values
    const double x_low = pow(10, -6);
    const double x_high = pow(10, +10);
    const double delta = pow(10, -6);

    // console output
    cout << "Relative error between the nummerically unstable function f and "
            "stable function g:"
         << endl;
    cout << "a) x = " << x_high << ", f(x) = " << f1(x_high)
         << " , g(x) = " << g1(x_high)
         << ", rel. error: " << rel_error(g1(x_high), f1(x_high)) << endl;
    cout << "b) x = " << x_low << ", f(x) = " << f2(x_low)
         << " , g(x) = " << g2(x_low)
         << ", rel. error: " << rel_error(g2(x_low), f2(x_low)) << endl;
    cout << "c) x = " << x_high << ", delta = " << delta
         << ", f(x) = " << f3(x_high, delta)
         << " , g(x) = " << g3(x_high, delta)
         << ", rel. error: " << rel_error(g3(x_high, delta), f3(x_high, delta))
         << endl;

    return 0;
}
