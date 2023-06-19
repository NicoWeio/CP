#ifndef OPTIMIZE_H
#define OPTIMIZE_H

#include <iostream>
#include <fstream>
#include <functional>

class Bisection{
private:
    // variables
    double x, y, z, eps;
    std::function<double(double)> f;

public:
    Bisection(std::function<double(double)> f, double x0, double y0, double z0, double eps);
    double min(std::string path_out);
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
    this->eps = eps;
}

double Bisection::min(std::string path_out){
    double u_np1;
    double error = z-x;
    int counter = 0;

    // write to csv file
    std::ofstream file(path_out);
    file << "step,x,y,z,error\n";
    file << counter << "," << x << "," << y << "," << z << "," << error << "\n";

    while(error > eps){
        counter++;
        if ((y-x) > (z-y)){ // half left side
            u_np1 = (y+x)/2.0;
            if (f(u_np1) < f(y)){
                z = y;
                y = u_np1;    
            }
            else{
                x = u_np1;
            }
        }
        else{ // half right side
            u_np1 = (y+z)/2.0;
            if (f(u_np1) < f(y)){
                x = y;
                y = u_np1;
            }
            else{
                z = u_np1;
            }

        }
        // calc error
        error = z-x;
        std::cout << "Error: " << error << " / " << eps << "\r";
    
        // write to csv file
        file << counter << "," << x << "," << y << "," << z << "," << error << "\n";
    }
    file.close();
    std::cout << "\nNumber of iterations: " << counter << std::endl;
   
    return y;
}

class Newton{
private:
    // variables
    double x, eps;
    std::function<double(double)> f;

public:
    Newton(std::function<double(double)> f, double x0, double eps);
    double min(std::string path_out);
};

Newton::Newton(std::function<double(double)> f, double x0, double eps){
    this->f = f;
    this->x = x0;
    this->eps = eps;
}

double Newton::min(std::string path_out){
    double h = 1E-5; // step size for derivative
    double df, df2;

    double error = 1.0;
    double x_nm1;
    int counter = 0;

    // write to csv file
    std::ofstream file(path_out);
    file << "step,x,error\n";
    file << counter << "," << x << "," << error << "\n";

    while(error > eps){
        counter++;
        // calc derivatives
        df = (f(x+h) - f(x-h)) / (2.0*h);
        df2 = (f(x+h) - 2.0*f(x) + f(x-h)) / (h*h);

        // calc next step
        x_nm1 = x;
        x = x - df/df2;

        // calc error
        if (x>x_nm1){
            error = x - x_nm1;
        }
        else{
            error = x_nm1 - x;
        }

        // write to csv file
        file << counter << "," << x << "," << error << "\n";

        std::cout << "Error: " << error << " / " << eps << "\r";
    }
    file.close();
    std::cout << "\nNumber of iterations: " << counter << std::endl;

    return x;
}

#endif