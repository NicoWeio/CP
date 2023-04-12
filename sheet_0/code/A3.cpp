#include <iostream>
#include <fstream>
#include <functional>
#include <math.h>
using namespace std;

void euler(float y0, float dt, string filename) {
    ofstream file(filename);
    float y = y0;
    file << "# t , y "<< endl;
    for(float t = 0 ; t<=10; t += dt ){
        file << t << "," << y << endl;
        y = y*(1 -dt );
    }
    file.close();
}

void symeuler(float y0, float y1, float dt, string filename) {
    ofstream file(filename);
    float yn_1 = y0;
    float yn = y1;
    float y_zwischen;
    file << "# t , y "<< endl;
    file <<  0 << "," << y0 << endl;
    file << dt << "," << y1<< endl;
    for(float t = dt*2 ; t<=10; t += dt ){
        y_zwischen=yn;
        yn = yn_1 - 2*dt*yn; 
        file << t << "," << yn<< endl;
        yn_1 = y_zwischen;
    }
    file.close();
}

int main() {
    float dt = 0.01;
    euler(1, dt, "build/A3_euler1.csv");
    symeuler(1, exp(-dt), dt, "build/A3_symeuler1.csv");
    euler(1-dt, dt, "build/A3_euler2.csv"); 
    symeuler(1, 1-dt, dt, "build/A3_symeuler2.csv");
    return 0;
}