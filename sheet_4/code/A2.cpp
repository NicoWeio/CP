#include <fstream>
#include <iostream>

using namespace std;

// class for the diffusion equation
class Diffusion {
    public:
        // parameters
        double D;
        double L;
        double dx;
        double dt;
        uint t_steps;
        uint x_steps;
        string initial;
        // initialize array u
        double *u;

        // methods
        Diffusion(double D, double L, string initial);
        void solve(double dt, double dx, uint t_steps, string path_write);
        void set_inital_conditions();
};

// constructor
Diffusion::Diffusion(double D, double L, string initial) {
    this->D = D;
    this->L = L;
    
}

// set inital conditions
void Diffusion::set_inital_conditions(){
    // initial condition
    u = new double[x_steps+1];

    for (uint i=0; i<=x_steps; i++) {
        u[i] = 1;
    }
}

// solve the diffusion equation
void Diffusion::solve(double dt, double dx, uint t_steps, string path_write) {
    // set parameters
    this->dx = dx;
    this->dt = dt;
    this->t_steps = t_steps;
    this->x_steps = uint(L/dx);

    set_inital_conditions();

    // constant
    double C = D*dt/(dx*dx);

    // edges are isolating, no flow through them
    u[0] = 0;
    u[x_steps] = 0;

    ofstream outfile;
    outfile.open(path_write);

    // FTCS scheme
    for (uint j=0; j<t_steps; j++) {
        for (uint i=1; i<x_steps; i++) {
            u[i] = C*u[i-1] + (1-2*C)*u[i] + C*u[i+1];
        }
        
        cout << "Progress: " << j << "/" << t_steps << "\r";
        cout.flush();

        outfile << j*dt;
        // write to csv file
        for (int i=0; i<=x_steps; i++) {
             outfile << "," << u[i];
        }
        outfile << endl;

    }
    
    outfile.close();
}

int main() {
    double D = 1;
    const double L = 1.0;
    double dx = 0.01;
    string path = "build/A1.csv";
    
    double dt = 0.00001;
    uint t_steps = 40000;
    

    Diffusion dif(D, L, "step");
    dif.solve(dt, dx, t_steps, path);

    
    return 0;
}
