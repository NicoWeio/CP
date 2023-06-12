#ifndef LCG_H
#define LCG_H

#include <fstream>
#include <iostream>
#include <math.h>

class LCG
{
private:
    double r_0;
    long long a;
    long long c;
    long long m;

    // number of random numbers
    int n;
    // random numbers
    double *r;

    // helper methods
    void scale(double lower, double upper);

public:
//a, c and, m is 64bit integer
    LCG(double r_0, long long a, long long c, long long m);
    double* uniform(int n, double lower, double upper);
    ~LCG();
};

LCG::LCG(double r_0, long long a, long long c, long long m){
    this->r_0 = r_0;
    this->a = a;
    this->c = c;
    this->m = m;

    this->n = 0;
    this->r = NULL;
}
LCG::~LCG(){
    if (this->r != NULL){
        delete[] this->r;
    }
}

void LCG::scale(double lower, double upper){
    // max value currently: m-1
    // scale r to [lower, upper]
    for (int i = 0; i<n; i++){
        r[i] = lower + (upper-lower)*r[i]/(m-1);       
    }
}

double* LCG::uniform(int n, double lower, double upper){
    this->n = n;
    r = new double[n]; // create array with size n

    // LCG algorithm
    double r0temp = r_0;
    for (int ind = 0; ind<n; ind++){  
        r[ind] = fmod((a*r0temp + c), m); // calc r_n+1
        r0temp = r[ind]; // update r_n
    }
    scale(lower, upper);
    return r;
}

#endif