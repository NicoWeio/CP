#include <fstream>
#include <functional>
#include <iostream>
#include <math.h>
#include <eigen3/Eigen/Dense>
using namespace std;

Eigen::Vector3d RungeKutter(double r0, double z0, double phi0, double h){
    Eigen::Vector3d y(r0,z0,phi0);
    Eigen::Vector3d k1= df(y);
    Eigen::Vector3d k2 = df ( y + h/2*k1);
    Eigen::Vector3d k3 = df ( y + h/2*k2);
    Eigen::Vector3d k3 = df ( y + h*k3);
    return y + h * ()




    return y1;
}


Eigen::Vector3d df(Eigen::Vector3d y){
    Eigen::Vector3d dy(cos(y(2)), sin(y(2)), 2 - 0.1* y(1)- sin(y(2))/y(0));
    return dy;
}


int main() {
    Eigen::Vector3d y(1,1,0);
   cout<< df(y);
  return 0;
}
