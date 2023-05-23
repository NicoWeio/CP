#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <fstream>

using Cmatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;
using Cdouble = std::complex<double>;

class Schroedinger1D
{
private:
    // grid
    double dx;
    double x_min;
    double x_max;
    int xsize;

    // time
    int N;
    double dt;

    // constants
    double dx2;

    // operators and observables: Complex matrix
    Cmatrix H; //Hamilton
    Cmatrix eye; //unit
    Cmatrix S; //time evolution

    // wave function psi_n, psi_n+1
    Cmatrix psi_n;
    Cmatrix psi_np1;

    // helper functions
    void init_H();
    void init_S();
    void write_to_csv(std::ofstream& file);
     

public:
    Schroedinger1D(double xmin, double xmax, double dx);
    void init_psi_gauss(double x0, double sigma);
    void run(double dt, int N, std::string path);
};

Schroedinger1D::Schroedinger1D(double xmin, double xmax, double dx){
    this->dx = dx;
    this->dx2 = dx*dx;

    this->x_min = xmin;
    this->x_max = xmax;

    xsize = int((x_max-x_min)/dx);
    std::cout << "\nA1) 1D Schrödinger Equation" << std::endl;
    std::cout << "Grid size: " << xsize << std::endl;


    // init unit matrix
    eye.resize(xsize, xsize);
    eye = Cmatrix::Identity(xsize, xsize);

    // initialize Hamilton
    H.resize(xsize, xsize);
    init_H(); 

    // init time evaluation operator (const for each time step?)
    S.resize(xsize, xsize);
    init_S();
    std::cout << "Number of entries in Hamilton operator: " << H.size() << std::endl;

    // init wave function
    psi_n.resize(xsize, 1);
    psi_np1.resize(xsize, 1);
}

void Schroedinger1D::init_H(){
    // not so bad if slow, using one time
    for (int n = 0; n<xsize; n++){
        for (int m = 0; m<xsize; m++){
            if(n == m-1 || n == m+1){ //side-diagonal
                H(n,m) = Cdouble(-1/dx2, 0.0);
            }
            else if (n == m){ //diagonal
                H(n,m) = Cdouble(2/dx2 + dx2 * n*n, 0.0); 
            }
            else{
                H(n,m) = Cdouble(0.0, 0.0);
            }
        }
    }
}

void Schroedinger1D::init_S(){
    // init time evolution operator
    Cmatrix term = Cdouble(0.0, 1.0) / 2.0 * H*dt;
    
    S = (eye + term).inverse() * (eye - term);
}

void Schroedinger1D::init_psi_gauss(double x0, double sigma){
    // init wave function
    // SET BY USER
    // normalized gaussian packet with position x0 and standard deviation sigma
    
    // check if sigma is positive
    if (sigma<=0){
        std::cerr << "sigma not positive" << std::endl;
        exit(1);
    }
    // check if x0 is in grid
    if (x0<x_min || x0>x_max){
        std::cerr << "x0 not in grid" << std::endl;
        exit(1);
    }

    // init wave function
    for (int n = 0; n<xsize; n++){
        psi_n(n,0) = Cdouble(pow(1.0/(2*M_PI*sigma), 0.25), 0.0) * exp(Cdouble(-pow(x_min + n*dx - x0, 2.0) / (4.0*sigma), 0.0));
    }
}

void Schroedinger1D::write_to_csv(std::ofstream& file){
    // write wave function to csv file
    for (int n = 0; n<xsize; n++){
        file << psi_n(n,0).real() << ",";
    }
    file << "\n";
}

void Schroedinger1D::run(double dt, int N, std::string path){
    // run simulation: time development of psi
    // SET BY USER
    this->dt = dt;
    this->N = N;

    std::cout << "Running simulation..." << std::endl;

    std::ofstream file(path);
    write_to_csv(file);
    for (int n = 0; n<=N; n++){
        std::cout << "Time step: " << n << "\r" << std::flush;
        psi_np1 = S * psi_n;
        psi_n = psi_np1;
        write_to_csv(file);
    }
    file.close();

    std::cout << "Simulation finished" << std::endl;
}

int main() {
    // parameters
    double xmin = -10.0;
    double xmax = 10.0;
    double dx = 1.0;
    double dt = 0.02;
    double x0 = 1.0;
    double sigma = 1.0;

    // run simulation
    Schroedinger1D simulation(xmin, xmax, dx);
    simulation.init_psi_gauss(x0, sigma);
    simulation.run(dt, 1E5, "build/A1.csv");

    return 0;
}
