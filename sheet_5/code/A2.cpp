#include <fstream>
#include <iostream>
#include <math.h>

class Wave2D
{
private:
    // wave speed
    const double c = 1.0;

    // field variables: u_n+1, u_n, u_n-1
    double** u_np1;
    double** u_n;
    double** u_nm1;

    // constants to speed up calculations
    double dx2;
    double dy2;
    double c2dt2;

    // helper functions
    void update_u_np1();
    void swap_u();
    void write_to_csv(std::ofstream& file);

public:
    // size of the membrane axb
    int a;
    int b;

    // space steps
    double dx;
    double dy;
    int xsize;
    int ysize;

    // time steps
    int n;
    double dt;

    Wave2D(double a, double b, double dx, double dy);
    void init_u();
    ~Wave2D();

    void run(int n, double dt, std::string path_out, int print_step=1);
};

Wave2D::Wave2D(double a, double b, double dx, double dy){
    // set parameters
    this->a = a;
    this->b = b;
    this->dx = dx;
    this->dy = dy;
    dx2 = dx*dx;
    dy2 = dy*dy;

    // grid size
    xsize = int(a / dx);
    ysize = int(b / dy);
    std::cout << "\nA2) 2D Wave Equation" << std::endl;
    std::cout << "Rectangle with a = " << a << ", b = " << b << "\n(Resolution: ["<< xsize << ", " << ysize << "])" <<"\n";

    init_u();
}

void Wave2D::init_u(){
    // initialize u_n+1, u_n, u_n-1
    u_np1 = new double*[xsize];
    u_n = new double*[xsize];
    u_nm1 = new double*[xsize];

    for (int i = 0; i < xsize; i++){
        u_np1[i] = new double[ysize];
        u_n[i] = new double[ysize];
        u_nm1[i] = new double[ysize];
    }

    // set initial conditions
    for (int i = 0; i < xsize; i++){
        for (int j = 0; j < ysize; j++){
            // x = dx * i
            // y = dy * j

            // calc u^1 with given u(x,y,t=0) = sin(pi*x/a)*sin(2*pi*y/b)
            u_nm1[i][j] = sin(M_PI*dx*i/a)*sin(2*M_PI*dy*j/b);

            // calc u^0 with given du/dt(x,y,t=0) = 0
            if (i == 0 || i == xsize-1 || j == 0 || j == ysize-1){
                u_n[i][j] = u_nm1[i][j]; // fixed boundaries
            }
            else{
                u_n[i][j] = c2dt2/2 * ( (u_nm1[i+1][j] - 2*u_nm1[i][j] + u_nm1[i-1][j])/dx2 + (u_nm1[i][j+1] - 2*u_nm1[i][j] + u_nm1[i][j-1])/dy2)
                            + u_nm1[i][j] + 2*dt*0; // du/dt|t=0 = 0
            }
        }
    }
}

Wave2D::~Wave2D(){
    // clean up storage
    for (int i = 0; i<xsize; i++){
        delete[] u_np1[i];
        delete[] u_n[i];
        delete[] u_nm1[i];
    }
    delete[] u_np1;
    delete[] u_n;
    delete[] u_nm1;
}

void Wave2D::update_u_np1(){
    // wave equation discretized and solved for u_n+1
    for (int i = 1; i < xsize-1; i++){
        for (int j = 1; j < ysize-1; j++){
            u_np1[i][j] = c2dt2* ( (u_n[i+1][j] - 2*u_n[i][j] + u_n[i-1][j])/dx2 + (u_n[i][j+1] - 2*u_n[i][j] + u_n[i][j-1])/dy2) + 2*u_n[i][j] - u_nm1[i][j];    
        }
    }
    // boundaries: fixed
    for (int i = 0; i < xsize; i++){
        u_np1[i][0] = u_n[i][0];
        u_np1[i][ysize-1] = u_n[i][ysize-1];
    }
    for (int j = 0; j < ysize; j++){
        u_np1[0][j] = u_n[0][j];
        u_np1[xsize-1][j] = u_n[xsize-1][j];
    }
}

void Wave2D::swap_u(){
    // u_n+1 -> u_n    u_n -> u_n-1
    for (int i = 0; i<xsize; i++){
        for (int j = 0; j<ysize; j++){
            u_nm1[i][j] = u_n[i][j];
            u_n[i][j] = u_np1[i][j];
        }
    }
}

void Wave2D::run(int n, double dt, std::string path_out, int print_step){
    // set time steps
    this->n = n;
    this->dt = dt;
    c2dt2 = c*c*dt*dt;

    std::ofstream file(path_out);

    // run simulation
    std::cout << "Running simulation...\n";

    for (int t_step = 0; t_step < n; t_step++){
        // print time step and write u_n to csv
        if (t_step % print_step == 0){
            std::cout << "t = " << t_step*dt << " / " << n*dt << "\r";
            std::cout.flush();

            write_to_csv(file);
        }
        // calc u_n+1
        update_u_np1();

        // u_n+1 -> u_n
        // u_n   -> u_n-1
        swap_u();
    }
    write_to_csv(file); // write last time step

    file.close();
    std::cout << "Done.\n" << "Data written to " << path_out << "\n";
}

void Wave2D::write_to_csv(std::ofstream& file){
    // write u_n to csv file
    for (int i = 0; i < xsize; i++){
        for (int j = 0; j < ysize; j++){
            file << u_n[i][j] << ",";
        }
        file << "|";
    }
    file << "\n";
}

int main(){
    // set parameters
    double a = 1.0; //1920.0;
    double b = 1.5; //1080.0;
    double dx = 0.01; // 1.0;
    double dy = dx;

    double dt = 0.000003;
    int n = 800000;

    // create Wave2D object
    Wave2D wave(a, b, dx, dy);
    wave.run(n, dt, "build/u_n.csv", 1000);

    return 0;
}
