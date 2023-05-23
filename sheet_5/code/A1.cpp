#include <Eigen/Dense>
using namespace Eigen;

class Schroedinger1D
{
private:
    // grid
    double dx;
    double x_min;
    double x_max;
    int xsize;

    // time
    int n;
    double dt;

    // operators and observables
    Eigen::Matrix2d H; //Hamilton TODO: symmetric??

    // helper functions
    void init_H();

public:
    Schroedinger1D(double xmin, double xmax, double dx);
    ~Schroedinger1D();
};

Schroedinger1D::Schroedinger1D(double xmin, double xmax, double dx){
    this->dx = dx;
    this->x_min = x_min;
    this->x_max = x_max;

    xsize = int(dx/(x_max-x_min));
    
    init_H();
}

void Schroedinger1D::init_H(){

}

Schroedinger1D::~Schroedinger1D(){

}

int main() {
    // parameters
    double xmin = -10.0;
    double xmax = 10.0;
    double dx = 0.1;
    double dt = 0.01;

    return 0;
}
