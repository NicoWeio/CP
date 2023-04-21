#include "integrate.h"
#include <iostream>
#include <math.h>

double f2(double x) {
    return x * sin(1 / x);
}

int main() {
    // test integrate.h
    std::cout << trapezregel(f2, 0.1, 1, 0.0001) << std::endl;
    std::cout << mittelpunktsregel(f2, 0.1, 1, 0.0001) << std::endl;
    std::cout << simpsonregel(f2, 0.1, 1, 0.0001) << std::endl;
    return 0;
}
