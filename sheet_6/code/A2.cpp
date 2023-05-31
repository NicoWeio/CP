#include <iostream>
#include <fstream>
#include <Eigen/Dense>

class LorenzModel{
private:
    // constants
    double sigma;
    double b;
    double r;

    // state
    Eigen::Vector3d x;

    // methods
    Eigen::Vector3d f(Eigen::Vector3d x);

public:
    LorenzModel(double sigma, double b, double r, Eigen::Vector3d x0);
    void rk4(double dt, int nsteps, std::string path_out);
};

LorenzModel::LorenzModel(double sigma, double b, double r, Eigen::Vector3d x0){
    this->sigma = sigma;
    this->b = b;
    this->r = r;
    this->x = x0;
}

Eigen::Vector3d LorenzModel::f(Eigen::Vector3d x){
    // f contains the derivatives
    Eigen::Vector3d f;
    f(0) = -sigma*x(0) + sigma*x(1);
    f(1) = -x(0)*x(2) + r*x(0) - x(1);
    f(2) = x(0)*x(1) - b*x(2);

    return f;
}


void LorenzModel::rk4(double dt, int nsteps, std::string path_out){
    Eigen::Vector3d k1, k2, k3, k4;

    // write x to csv
    std::ofstream file(path_out);
    file << "x,y,z\n";
    file << x(0) << "," << x(1) << "," << x(2) << "\n";

    // rk4
    std::cout << "Run Runge-Kutta method (4th order) ..." << std::endl;
    for (int n = 0; n<=nsteps; n++) {
        std::cout << n << " / " << nsteps << "\r" << std::flush;

        k1 = dt*f(x);
        k2 = dt*f(x + 0.5*k1);
        k3 = dt*f(x + 0.5*k2);
        k4 = dt*f(x + k3);

        x = x + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4);


        file << x(0) << "," << x(1) << "," << x(2) << "\n"; //write to csv
    }
    std::cout << std::endl << "Done." << std::endl;
}

int main() {
    // parameters
    double sigma = 10.0;
    double b = 8.0/3.0;
    double r1 = 20.0;
    double r2 = 28.0;
    
    double dt = 0.01;
    int nsteps = 10000;

    // initial conditions
    Eigen::Vector3d x0;
    Eigen::Vector3d x1;
    x0 << 1.0, 1.0, 1.0;
    x1 << -5.0, -1.0, 1.0;

    // create model
    LorenzModel model1(sigma, b, r1, x0);
    model1.rk4(dt, nsteps, "build/A2_r20.csv");

    LorenzModel model2(sigma, b, r2, x0);
    model2.rk4(dt, nsteps, "build/A2_r28.csv");

    return 0;
}
