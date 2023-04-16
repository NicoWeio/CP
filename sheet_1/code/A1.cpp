#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
using namespace std;

float df_zwei(std::function<float(float)> f, float x, float h) {
    return (f(x+h) - f(x-h))/(2*h);
}


float df2_zwei(std::function<float(float)> f, float x, float h) {
    return (df_zwei(f, x+h, 0.01) - df_zwei(f, x-h, 0.01))/(2*h);
}

float df_vier(std::function<float(float)> f, float x, float h) {
    return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/(12*h);
}

float f1(float x) {
    return sin(x);
}


float f2(float x) {
    if (x >= 0) {
        return 2 * floor(x/M_PI) - cos( fmod(x, M_PI)) +1;       
    } else {
        return 2 * floor(x/M_PI) + cos( fmod(x, M_PI)) +1;
    }
}


int main() {
//a)
ofstream File1("build/A1_a_h.csv");
    for (float h = 0.000000001; h <= 1.0; h *=10) {
        File1 << h << "," << df_zwei(f1, 1.5, h) << endl;
    }
    File1.close();

ofstream File2("build/A1_a_x.csv");
    for (float x = -M_PI; x <= M_PI; x +=0.1) {
        File2 << x << "," << df_zwei(f1, x, 0.01) << endl;
    }
    File2.close();

//b)
ofstream File3("build/A1_b_h.csv");
    for (float h =  0.000000001; h <= 1.0; h *=10) {
        File3 << h << "," << df2_zwei(f1, 1.5, h) << endl;
    }
    File3.close();

ofstream File4("build/A1_b_x.csv");
    for (float x = -M_PI; x <= M_PI; x +=0.1) {
        File4 << x << "," << df2_zwei(f1, x, 0.01) << endl;
    }
    File4.close();

//c)
ofstream File5("build/A1_c_h.csv");
    for (float h = 0.000000001; h <= 1.0; h *=10) {
        File5 << h << "," << df_vier(f1, 1.5, h) << endl;
    }
    File5.close();

ofstream File6("build/A1_c_x.csv");
    for (float x = -M_PI; x <= M_PI; x +=0.1) {
        File6 << x << "," << df_vier(f1, x, 0.01) << endl;
    }
    File6.close();


return 0;
}

