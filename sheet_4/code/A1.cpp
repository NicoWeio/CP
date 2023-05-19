#include <iostream>
#include <fstream>
#include <Eigen/Dense>


class Poisson2d{
private:
    int xsize;
    int ysize;
    double** phi; // 2d array: shape(xsize, ysize)
    double** phi_old;
    double** rho;

    void set_initial_conditions();
    double calc_error();

public:
    Poisson2d(int xsize, int ysize);
    ~Poisson2d();
    void gauss_seidel(double delta, double kappa);
};
    
Poisson2d::Poisson2d(int xsize, int ysize){
    // set values
    this->xsize = xsize;
    this->ysize = ysize;

    // create 2d array
    phi = new double*[xsize];
    phi_old = new double*[xsize];


    for (int i = 0; i<xsize; i++){
        phi[i] = new double[ysize];
        phi_old[i] = new double[ysize];
    }

    // TODO: set_initial_conditions for phi and rho
    // use function set_initial_conditions()
}

Poisson2d::~Poisson2d(){
    // clear memory
    for (int i = 0; i<xsize; i++){
        delete[] phi[i];
        delete[] phi_old[i];

    }
    delete[] phi;
    delete[] phi_old;
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

    // TODO: calc relative error
    error = calc_error();

    } while (error>kappa);
    

}

double Poisson2d::calc_error(){
    //TODO: calc error
    return 0.0;
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