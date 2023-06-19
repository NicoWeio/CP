#include <iostream>
#include <functional>

class Bisection
{
private:
    // variables
    double x, y, z, eps;
    std::function<double(double)> f;

    // functions
    double min(std::function<double(double)> f0);

public:
    Bisection(std::function<double(double)> f, double x0, double y0, double z0, double eps);
};

Bisection::Bisection(std::function<double(double)> f, double x0, double y0, double z0, double eps){
    if(x0>y0 || y0>z0 || eps<=0){
        std::cerr << "Invalid arguments!" << std::endl;
        return;
    }
    this->x = x0;
    this->y = y0;
    this->z = z0;
    this->f = f;
}

double Bisection::min(std::function<double(double)> f){
    double u_left, u_right;


}

int main(){

    return 0;
}