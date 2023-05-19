#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <math.h>

class Poisson2d{
private:
    int xsize;
    int ysize;
    double** phi; // 2d array: shape(xsize, ysize)
    double** phi_old;
    double** rho;

    double calc_error();

public:
    Poisson2d(int xsize, int ysize);
    ~Poisson2d();
    void set_initial_conditions();
    void gauss_seidel(double delta, double kappa);
};
    
Poisson2d::Poisson2d(int xsize, int ysize){
    // set values
    this->xsize = xsize;
    this->ysize = ysize;

    // create 2d arrays: phi_n+1, phi_n, rho
    phi = new double*[xsize];
    phi_old = new double*[xsize];
    rho = new double*[xsize];

    for (int i = 0; i<xsize; i++){
        phi[i] = new double[ysize];
        phi_old[i] = new double[ysize];
        rho = new double*[xsize];
    }

    // TODO: set_initial_conditions for phi and rho
    // use function set_initial_conditions()
}

Poisson2d::~Poisson2d(){
    // clear memory
    for (int i = 0; i<xsize; i++){
        delete[] phi[i];
        delete[] phi_old[i];
        delete[] rho[i];
    }

    delete[] phi;
    delete[] phi_old;
    delete[] rho;
}

void Poisson2d::set_initial_conditions(){
    // TODO: initial/boundary conditions
    // calculate rho[i][j]
}

void Poisson2d::gauss_seidel(double delta, double kappa){
    // kappa: stop criteria
    // relative error between to time steps
    double error;
    
    // copy phi to phi_old
    for (int i = 0; i<xsize; i++){
        for (int j = 0; j<ysize; j++){
            phi_old[i][j] = phi[i][j];
        }
    }

    do{
        for (int i = 1; i<xsize-1; i++){
            for (int j = 1; j<ysize-1; j++){
                phi[i][j] = 0.25 * (phi[i+1][j] + phi[i][j+1] + phi[i][j-1] + phi[i-1][j]) + 0.25 * delta*delta*rho[i][j];
            }
        }

    // calc relative change of euclidean distance between phi and phi_old
    error = calc_error();
    std::cout << "Error: " << error << std::endl;
    } while (error>kappa);
}

double Poisson2d::calc_error(){
    // calc rel error between phi and phi_old
    // distance measure: euclidean distance

    double norm_phi = 0.0;
    double norm_phiold = 0.0;

    for (int i = 0; i<xsize; i++){
        for (int j = 0; j<ysize; j++){
            norm_phi += phi[i][j]*phi[i][j];
            norm_phiold += phi_old[i][j]*phi_old[i][j];
        }
    }
    norm_phi = std::sqrt(norm_phi);
    norm_phiold = std::sqrt(norm_phiold);

    return std::abs(norm_phi-norm_phiold) / norm_phi;
}





    
int main(){
    // parameters
    int xsize, ysize; //square
    xsize = 1;
    ysize = 1;

    double delta = 0.05;
    double kappa = 0.00001;

    Poisson2d poisson(xsize, ysize);
    poisson.gauss_seidel(delta, kappa);
    return 0;
}