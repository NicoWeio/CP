#include <iostream>
#include <Eigen/Dense>

class EVpower // determines ev for 4x4 matrix
{
private:
    Eigen::Matrix4f A;
    Eigen::Vector4f v_n;
    Eigen::Vector4f w_n;
    int maxiter;

    void powerMethod();

public:
    EVpower(Eigen::Matrix4f A, Eigen::Vector4f x0, int maxiter);
    double getEigenvalue();
    Eigen::Vector4f getEigenvector();
};

EVpower::EVpower(Eigen::Matrix4f A, Eigen::Vector4f x0, int maxiter){
    this->A = A;
    this->v_n = x0;
    this->maxiter = maxiter;

    powerMethod();
}

void EVpower::powerMethod(){
    for (int i = 0; i < maxiter; i++){
        w_n = A*v_n;
        v_n = w_n/w_n.norm();
        std::cout << "Step: " << i+1 << " / " << maxiter << "\r";    
        //std::cout << "v_n:\n" << v_n << std::endl;
    }
}

double EVpower::getEigenvalue(){
    return (v_n.transpose() * A * v_n)(0);
}

Eigen::Vector4f EVpower::getEigenvector(){
    return v_n;
}


int main(){
    // define given matrix A
    Eigen::Matrix4f A;
    A << 1, -2, -3, 4,
        -2, 2, -1, 7,
        -3, -1, 3, 6,
         4, 7, 6, 4;

    // calc ev with Eigen::Eigenvalues()
    auto ev_EIGEN = A.eigenvalues();
    std::cout << "Exercise 1: Matrix diagonalization – Power method" << std::endl;
    std::cout << "Eigenvalues determined by Eigen:\n" << ev_EIGEN << std::endl;
    
    // calc ev with power method
    Eigen::Vector4f x0;
    x0 << 1, 1, 1, 1;

    EVpower evpower(A, x0, 1000);
    double eval_power = evpower.getEigenvalue();
    Eigen::Vector4f evec_power = evpower.getEigenvector();
    std::cout << "Eigenvalues determined by power method:\n" << eval_power << std::endl;
    std::cout << "Eigenvector determined by power method:\n" << evec_power << std::endl;
    return 0;
}