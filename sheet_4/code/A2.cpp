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
        Diffusion(double D, double L, double dx);
        void solve(double dt, uint t_steps, string path_write);

        // set initial conditions: must be called before solve()
        void set_initial_const(double a); // a) set initial condition to constant 
        void set_initial_delta(double x0, double a); // b) set initial condition to delta function with height "a"
        void set_initial_heaviside(double x0, double a); //c) Heaviside function with height "a"
};

// constructor
Diffusion::Diffusion(double D, double L, double dx) {
    // set parameters
    this->D = D;
    this->L = L;
    this->dx = dx;
    this->x_steps = uint(L/dx);
}

// set inital conditions
void Diffusion::set_initial_const(double a){
    // u_0 = a 
    u = new double[x_steps+1];

    for (uint i=0; i<=x_steps; i++) {
        u[i] = a;
    }
}

void Diffusion::set_initial_delta(double x0, double a){
    // u_0 = delta(x-a) with height a
    u = new double[x_steps+1];
    
    // get bin position of x0
    uint pos = uint(x0/dx);

    for (uint i=0; i<=x_steps; i++){
        if (i==pos){
            u[i] = a;
        }
        else{
            u[i] = 0.0;
        }
    }
}

// solve the diffusion equation
void Diffusion::solve(double dt, uint t_steps, string path_write) {
    // set parameters
    this->dt = dt;
    this->t_steps = t_steps;

    // constant
    double C = D*dt/(dx*dx);


    ofstream outfile;
    outfile.open(path_write);

    // FTCS scheme
    for (uint j=0; j<t_steps; j++) {
        // boundary conditions: reflective

        for (uint i=1; i<x_steps-1; i++) {
            u[i] = u[i] + C*(u[i+1] - 2*u[i] + u[i-1]);
        }

        // boundary conditions: reflective
        u[0] = u[1];
        u[x_steps-1] = u[x_steps-2];

        // print progress
        cout << "Progress: " << j << "/" << t_steps << "\r";
        cout.flush();

        // write to csv file
        outfile << j*dt;
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
    
    // a)
    double dt = 0.00001;
    uint t_steps = 40000;
    
    Diffusion dif_a(D, L, dx);
    dif_a.set_initial_const(1);
    dif_a.solve(dt, t_steps, "build/A2_a).csv");
    

    // b)
    // dt fullfills stability condition -> stable
    dt = 0.00001;
    t_steps = 1000;

    Diffusion dif_b_slow(D, L, dx);
    dif_b_slow.set_initial_delta(0.5, 1);
    dif_b_slow.solve(dt, t_steps, "build/A2_b)slow.csv");

    // dt too large -> unstable
    dt = 0.0001;
    t_steps = 80;

    Diffusion dif_b_fast(D, L, dx);
    dif_b_fast.set_initial_delta(0.5, 1);
    dif_b_fast.solve(dt, t_steps, "build/A2_b)fast.csv");

    
    return 0;
}
