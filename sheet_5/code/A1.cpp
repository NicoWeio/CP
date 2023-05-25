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
    void prob_to_csv(std::ofstream& file);
     

public:
    Schroedinger1D(double xmin, double xmax, double dx, double dt);
    void init_psi_gauss(double x0, double sigma);
    void run(int N, std::string path);
    double calc_rho();
};

Schroedinger1D::Schroedinger1D(double xmin, double xmax, double dx, double dt){
    this->dt = dt;
    this->dx = dx;
    this->dx2 = dx*dx;

    this->x_min = xmin;
    this->x_max = xmax;

    this->xsize = int((x_max-x_min)/dx);
    std::cout << "\nA1) 1D Schrödinger Equation" << std::endl;
    std::cout << "Grid size: " << xsize << std::endl;


    // init unit matrix
    eye.resize(xsize, xsize);
    eye = Cmatrix::Identity(xsize, xsize);

    // init Hamilton
    H.resize(xsize, xsize);
    init_H(); 
    std::cout << "Number of entries in Hamilton operator: " << H.size() << std::endl;

    // init time evaluation operator (const for each time step)
    S.resize(xsize, xsize);
    init_S();
    
    // init wave function
    psi_n.resize(xsize, 1);
    psi_np1.resize(xsize, 1);
}

void Schroedinger1D::init_H(){
    // init Hamilton operator for 1D harmonic oscillator
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

    // print H to csv
    std::ofstream file("build/A1_H.csv");
    for (int n = 0; n<xsize; n++){
        for (int m = 0; m<xsize; m++){
            file << H(n,m).real() << ",";
        }
        file << "\n";
    }
}

void Schroedinger1D::init_S(){
    // init time evolution operator
    Cmatrix term = dt/2 *Cdouble(0.0, 1.0) * H;

    S = (eye + term).inverse() * (eye - term);

    // print S to csv
    std::ofstream file("build/A1_S_real.csv");
    for (int n = 0; n<xsize; n++){
        for (int m = 0; m<xsize; m++){
            file << S(n,m).real() << ",";
        }
        file << "\n";
    }
    file.close();

    file.open("build/A1_S_imag.csv");
    for (int n = 0; n<xsize; n++){
        for (int m = 0; m<xsize; m++){
            file << S(n,m).imag() << ",";
        }
        file << "\n";
    }
    file.close();
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

    // init wave function: gaussian packet with position x0 and standard deviation sigma normalized to 1
    for (int n = 0; n<xsize; n++){
        psi_n(n,0) = Cdouble(pow(1.0/(2*M_PI*sigma), 0.25), 0.0) * exp(Cdouble(-pow(x_min + n*dx - x0, 2.0) / (4.0*sigma), 0.0));
    }

    // normalize
    psi_n = psi_n / sqrt(calc_rho());
}

double Schroedinger1D::calc_rho(){
    // calc probability density: rho = |psi|^2
    double sum = 0.0;
    for (int n = 0; n<xsize; n++){
        sum += psi_n(n,0).real() * psi_n(n,0).real() + psi_n(n,0).imag() * psi_n(n,0).imag();
    }
    return sum;
}

void Schroedinger1D::prob_to_csv(std::ofstream& file){
    // write propability density to csv file
    double rho = 0.0;

    for (int n = 0; n<xsize; n++){
        rho = psi_n(n,0).real() * psi_n(n,0).real() + psi_n(n,0).imag() * psi_n(n,0).imag();
        file << rho << ",";
    }
    file << "\n";
}

void Schroedinger1D::run(int N, std::string path){
    // run simulation: time development of psi
    // SET BY USER
    this->N = N;

    // sum of propability density, should be 1
    double rho_max;
    double rho_min;
    double rho;

    std::cout << "Running simulation..." << std::endl;
    std::ofstream file(path);

    // initial state
    prob_to_csv(file);
    rho = calc_rho();
    rho_max = rho;
    rho_min = rho;

    // simulation
    for (int n = 0; n<=N; n++){
        std::cout << "Time step: " << n << "\r" << std::flush;

        // time development
        psi_np1 = S * psi_n;
        psi_n = psi_np1;
        
        prob_to_csv(file);

        // check if sum of rho = |psi|^2 is equal to 1
        rho = calc_rho();
        if (rho > rho_max){
            rho_max = rho;
        }
        if (rho < rho_min){
            rho_min = rho;
        }
    }
    file.close();

    std::cout << "Sum of probability density rho should be 1 for each time step." << std::endl;
    std::cout << "Max: " << rho_max << std::endl;
    std::cout << "Min:" << rho_min << std::endl;

    std::cout << "Simulation finished" << std::endl;
}

int main() {
    // parameters
    double xmin = -10.0;
    double xmax = 10.0;
    double dx = 0.1;
    double dt = 0.02;
    double x0 = 1.0;
    double sigma = 1.0;

    // run simulation
    Schroedinger1D simulation(xmin, xmax, dx, dt);
    simulation.init_psi_gauss(x0, sigma);
    simulation.run(1E3, "build/A1_psi.csv");

    return 0;
}
