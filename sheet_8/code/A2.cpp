#include "LCG.h"
#include <fstream>
#include <iostream>
#include <math.h>

// print array
void printarr(double* r, int size){
    std:: cout << "[";
    for (int i = 0; i<size; i++){
        std::cout << r[i];
        if (i!=size-1){
            std::cout << ", ";
        }
    }
    std:: cout << "]\n";
}

// write array to file
void writearr(double* r, int size, std::string filename){
    std::ofstream file;
    file.open(filename);
    for (int i = 0; i<size; i++){
        file << r[i] << "\n";
    }
    file.close();
}

// function for c)
double f_sin(double x){
    return sin(x)/2;
}

// inverse function for d)
// f = 3*x^2
double f_inv_dev(double x){
    //p(x)*dx/dy
    return pow(x, 1.0/3.0);
}

int main(){
    // set seed (Matlab)
    double r0 = r0 = 1234;
    long long a = 16807;
    long long c = 0;
    long long m = 2147483647;
    LCG rng(r0, a, c, m);

    int N = 1E5;

    // a) Box-Mueller
    double* r_norm = rng.normal(N, "boxmueller");
    writearr(r_norm, N, "build/A2_boxmueller.csv");

    // b) Central limit theorem
    double* r_central = rng.normal(N, "central");
    writearr(r_central, N, "build/A2_central.csv");
    
    // c) Neumann rejection
    double* r_neumann = rng.neumann(N, 0.0, M_PI, f_sin);
    writearr(r_neumann, N, "build/A2_neumann.csv");

    // d) Inverse method
    double* r_inv = rng.inverse(N, 0.0, 3.0, f_inv_dev);
    writearr(r_inv, N, "build/A2_inverse.csv");

    return 0;
}