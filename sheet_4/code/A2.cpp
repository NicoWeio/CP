#include <fstream>
#include <iostream>

using namespace std;

// class for the diffusion equation
class Diffusion {
    private:
        void set_delta(double x0, double a);
    public:
        // parameters
        double D;
        double L;
        double dx;
        double dt;
        uint t_steps;
        uint x_steps;
        double *u; // declare 1d array u

        // methods
        Diffusion(double D, double L, double dx);
        ~Diffusion();
        void solve(double dt, uint t_steps, string path_write);
        
        // set initial conditions: must be called before solve()
        void set_initial_const(double a); // a) set initial condition to constant 
        void set_initial_delta(double x0, double a); // b)c) set initial condition to delta function with height "a"
        void set_initial_heaviside(double x0, double a); //c) Heaviside function with height "a"
        void set_initial_dirac_ridge(double x0, uint N); // c) Dirac ridge with N peaks at x0*n
};

// constructor
Diffusion::Diffusion(double D, double L, double dx) {
    // set parameters
    this->D = D;
    this->L = L;
    this->dx = dx;
    this->x_steps = uint(L/dx);
}

Diffusion::~Diffusion(){
    // clear memory
    delete[] u;
}

// set inital conditions
void Diffusion::set_initial_const(double a){
    // u_0 = a 
    u = new double[x_steps+1];

    for (uint i=0; i<=x_steps; i++) {
        u[i] = a;
    }
}

// helper function to set delta peaks
void Diffusion::set_delta(double x0, double a){ 
    // get bin position of x0
    uint pos = uint(x0/dx); 
    
    // set delta peak u_0 = delta(x-a) with height a
    u[pos] = a;
}

void Diffusion::set_initial_delta(double x0, double a){
    if (x0>L || x0<0){
        cerr << "x0 musst be smaller than L!" << endl;
        exit(1);
        }
    
    // initialize empty array
    u = new double[x_steps+1];
    for (uint i=0; i<=x_steps; i++){
        u[i] = 0.0;
        }

    // set delta peak u_0 = delta(x-a) with height a
    set_delta(x0, a);
}

void Diffusion::set_initial_heaviside(double x0, double a){
    if (x0>L || x0<0){
        cerr << "x0 musst be smaller than L!" << endl;
        exit(1);
        }

    // theta(x-x0) * a, a:height
    u = new double[x_steps+1];
    uint pos = uint(x0/dx);

    for (uint i=0; i<=pos; i++){
        u[i] = 0;
    }
    
    for (uint i=pos+1; i<=x_steps; i++){
        u[i] = 1;
    }
}

void Diffusion::set_initial_dirac_ridge(double x0, uint N){
    if (x0>L || x0<0){
        cerr << "x0 musst be smaller than L!" << endl;
        exit(1);
        }

    // initialize empty array
    u = new double[x_steps+1];
    for (uint i=0; i<=x_steps; i++){
        u[i] = 0;
    }

    // set delta peaks with height a=1/N
    for (uint n=1; n<=N; n++){
        set_delta(x0*n, 1.0/N);
    }
}

// solve the diffusion equation
void Diffusion::solve(double dt, uint t_steps, string path_write) {
    // check if initial conditions are set
    if (u == NULL) {
        cerr << "Initial conditions not set!" << endl;
        exit(1);
    }

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

        for (uint i=1; i<=x_steps; i++) {
            u[i] = u[i] + C*(u[i+1] - 2*u[i] + u[i-1]);
        }

        // boundary conditions: reflective
        u[0] = u[1];
        u[x_steps] = u[x_steps-1];

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

    // c)
    // Delta-distribution u1 = delta(x-0.5)
    dt = 0.00001;
    t_steps = 10000;

    Diffusion dif_u1(D, L, dx);
    dif_u1.set_initial_delta(0.5, 1);
    dif_u1.solve(dt, t_steps, "build/A2_c)u1.csv");

    // Heaviside u2 = theta(x-0.5)
    dt = 0.00001;
    t_steps = 100000;

    Diffusion dif_u2(D, L, dx);
    dif_u2.set_initial_heaviside(0.5, 1);
    dif_u2.solve(dt, t_steps, "build/A2_c)u2.csv");

    // Dirac ridge
    dt = 0.000001;
    t_steps = 20000;
    Diffusion dif_u3(D, L, dx);
    dif_u3.set_initial_dirac_ridge(0.1, 9);
    dif_u3.solve(dt, t_steps, "build/A2_c)u3.csv");
    return 0;
}
