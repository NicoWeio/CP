#ifndef LCG_H
#define LCG_H

#include <fstream>
#include <iostream>
#include <functional>
#include <math.h>

class LCG
{
private:
    double r_0;
    long long a;
    long long c;
    long long m;

    // random numbers
    double *r;

    // helper methods
    void scale(int n, double lower, double upper);

public:
//a, c and, m is 64bit integer
    LCG(double r_0, long long a, long long c, long long m);
    double* uniform(int n, double lower, double upper);
    double* normal(int n, std::string method);
    double* neumann(int n, double lower, double upper, std::function<double(double)> f);
    double* inverse(int n, double lower, double upper, std::function<double(double)> f);
    double* boxmueller(int n);
    double* centrallimit(int n);
    ~LCG();
};

LCG::LCG(double r_0, long long a, long long c, long long m){
    this->r_0 = r_0;
    this->a = a;
    this->c = c;
    this->m = m;

    this->r = NULL;
}
LCG::~LCG(){
    if (this->r != NULL){
        delete[] this->r;
    }
}

void LCG::scale(int n, double lower, double upper){
    // max value currently: m-1
    // scale r to [lower, upper]
    for (int i = 0; i<n; i++){
        r[i] = lower + (upper-lower)*r[i]/(m-1);       
    }
}

double* LCG::uniform(int n, double lower, double upper){
    r = new double[n]; // create array with size n

    // LCG algorithm
    for (int ind = 0; ind<n; ind++){  
        r[ind] = fmod((a*r_0 + c), m); // calc r_n+1
        r_0 = r[ind]; // update r_0 (seed)
    }
    scale(n, lower, upper);
    return r;
}

double* LCG::normal(int n, std::string method="boxmueller"){
    if (method == "boxmueller"){
        return boxmueller(n);
    }
    else if (method == "central"){
        return centrallimit(n);
    }
    else {
        std::cerr << "Error: method not found" << std::endl;
        return NULL;
    }
}

double* LCG::boxmueller(int n){
    double* r = new double[n];
    double v1, v2;
    double R2, R;
    double y1, y2;
 
    // polar method
    for (int ind = 0; ind<n;){
        // get pair of random numbers
        double v1 = uniform(1, -1, 1)[0];
        double v2 = uniform(1, -1, 1)[0];
        
        // check if inside unit circle
        R2 = v1*v1 + v2*v2;
        if (R2 < 1){
            R = sqrt(R2);
            y1 = sqrt(-2*log(R2)) * v1/R;
            y2 = sqrt(-2*log(R2)) * v2/R;

            // store in array
            r[ind] = y1;
            if (ind+1 < n){
                r[ind+1] = y2;
            }
            ind += 2;
        }
    }

    return r;
}

double* LCG::centrallimit(int n){
    int const N = 100; //number of unitnumbers per normal distributed number
    double* r = new double[n];

    double* r_unit;
    double y;

    for (int i = 0; i<n; i++){ // iterate over normal distributed numbers
        r_unit = uniform(N, 0, 1);

        // sum of unit drawn numbers - N/2
        y = -N/2;
        for (int j = 0; j<N; j++){
            y += r_unit[j];
        }
        r[i] = y;
    }

    return r;
}

double* LCG::neumann(int n, double lower, double upper, std::function<double(double)> f){
    double* r = new double[n];
    double x, y; // random numbers x, y
    
    for (int i = 0; i<n; i++){
        // take x as random number if y>f(x)
        do {
            x = uniform(1, lower, upper)[0];
            y = uniform(1, 0.0, 1.0)[0];

        } while(y > f(x));
        r[i] = x;
    }
    return r;
}

double* LCG::inverse(int n, double lower, double upper, std::function<double(double)> f){
    double* r = new double[n];
    double x, y; // random numbers x, y
    
    for (int i = 0; i<n; i++){
        // take x as random number if y>f(x)
        do {
            y = uniform(1, 0.0, 1.0)[0];
            x = f(y);

        } while(x < lower || x > upper);
        r[i] = x;
    }
    return r;
}


#endif