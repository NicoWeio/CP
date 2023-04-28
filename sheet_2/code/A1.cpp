#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/Dense>
using namespace std;

Eigen::Vector3d df(Eigen::Vector3d y){
    Eigen::Vector3d dy;
    if(y(0)==0){
        dy(0) = cos(y(2));
        dy(1) = sin(y(2));
        dy(2) = (2 - 0.1* y(1))/2;
    }
    else{
        dy(0) = cos(y(2));
        dy(1) = sin(y(2));
        dy(2) = 2 - 0.1* y(1)- sin(y(2))/y(0);
    }
    return dy;
}

Eigen::Vector3d RungeKutta(Eigen::Vector3d y0, double h){
    Eigen::Vector3d y = y0;
    Eigen::Vector3d k1= df(y);
    Eigen::Vector3d k2 = df ( y + h/2*k1);
    Eigen::Vector3d k3 = df ( y + h/2*k2);
    Eigen::Vector3d k4 = df ( y + h*k3);
    return y + h * (k1 + 2*k2+2*k3 + k4)/6;
}

int main() {
    double h =0.01;
    Eigen::Vector3d y(0,0,0);
    Eigen::Vector3d yn(0,0,0);
    ofstream File("build/A1.csv");
    for(double i=0; yn(0)< 1.5; i+=h){
    yn=RungeKutta(y, h);
    File << i*h << "," << y(0)<< ","<< yn(1)<< "," << yn(2)<< endl;
    y = yn;
    }
    File.close();
  return 0;
}
