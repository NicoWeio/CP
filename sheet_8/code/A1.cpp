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

int main() {
    // number of parameters
    int N = 1E5;
    double min = 0.0;
    double max = 1.0;

    // (a) r0 = 1234, a = 20, c = 120, m = 6075
    double r0 = 1234; 
    long long a = 20;
    long long c = 120;
    long long m = 6075;

    LCG rng_a(r0, a, c, m);
    double* r_a = rng_a.uniform(N, min, max);
    writearr(r_a, N, "build/A1_a.csv");
    
    
    // (b) r0 = 1234, a = 137, c = 187, m = 256
    r0 = 1234;
    a = 137;
    c = 187;
    m = 256;

    LCG rng_b(r0, a, c, m);
    double* r_b = rng_b.uniform(N, min, max);
    writearr(r_b, N, "build/A1_b.csv");

    // (c) r0 = 123456789, a = 65539, c = 0, m = 2^31 = 2147483648 (RANDU Generator von IBM)
    r0 = 123456789;
    a = 65539;
    c = 0;
    m = 2147483648;
    
    LCG rng_c(r0, a, c, m);
    double* r_c = rng_c.uniform(N, min, max);
    writearr(r_c, N, "build/A1_c.csv");

    // (d) r0 = 1234, a = 7^5 = 16807, c = 0, m = 2^31 − 1 (ran1() aus Numerical Recipes, 2. Ausgabe, bzw. Matlab bis Version 4)
    r0 = 1234;
    a = 16807;
    c = 0;
    m = 2147483647;

    LCG rng_d(r0, a, c, m);
    double* r_d = rng_d.uniform(N, min, max);
    writearr(r_d, N, "build/A1_d.csv");
    
    return 0;
}
