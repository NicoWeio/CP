#include <iostream>
#include <Eigen/Dense>

class EVpower // determines eigenvalues and -vectors for 4x4 matrix
{
private:
    Eigen::Matrix4f A;
    Eigen::Vector4f x0;
    int maxiter;

    Eigen::Vector4f v_n;
    Eigen::Vector4f w_n;

    //results
    Eigen::Matrix4f eigenvector;
    Eigen::Vector4f eigenvalue;

    void powerMethod();

public:
    EVpower(Eigen::Matrix4f A, Eigen::Vector4f x0, int maxiter);
    Eigen::Vector4f getEigenvalue();
    Eigen::Matrix4f getEigenvector();
};

EVpower::EVpower(Eigen::Matrix4f A, Eigen::Vector4f x0, int maxiter){
    this->A = A;
    this->x0 = x0;
    this->maxiter = maxiter;

    powerMethod();
}

void EVpower::powerMethod(){
    for (int ind_ev = 0; ind_ev<4; ind_ev++){ // loop over eigenvalues
        v_n = x0;
        for (int i = 0; i < maxiter; i++){ // loop over iterations
            w_n = A*v_n;
            v_n = w_n/w_n.norm();
        }
        // save results
        eigenvector.col(ind_ev) = v_n;
        eigenvalue(ind_ev) = (v_n.transpose() * A * v_n)(0);

        // adjust A for next iteration (set current eigenvalue to 0)
        A = A - eigenvalue(ind_ev) * v_n * v_n.transpose();
    }
}

Eigen::Vector4f EVpower::getEigenvalue(){
    return eigenvalue;
}

Eigen::Matrix4f EVpower::getEigenvector(){
    return eigenvector;
}


int main(){
    // define given matrix A
    Eigen::Matrix4f A;
    A << 1, -2, -3, 4,
        -2, 2, -1, 7,
        -3, -1, 3, 6,
         4, 7, 6, 4;

    // calc eigenvalues and -vectors with Eigen
    auto eval_EIGEN = A.eigenvalues();
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix4f> eigensolver(A);
    Eigen::Matrix4f evec_EIGEN = eigensolver.eigenvectors();
    std::cout << "Exercise 1: Matrix diagonalization – Power method" << std::endl;
    std::cout << "Eigenvalues determined by Eigen:\n" << eval_EIGEN << std::endl;
    std::cout << "Eigenvectors determined by Eigen:\n" << evec_EIGEN << "\n\n";

    // define initial vector x0
    Eigen::Vector4f x0;
    x0 << 1, 1, 1, 1;

    // calc eigenvalues and -vectors with power method
    EVpower evpower(A, x0, 100);
    Eigen::Vector4f eval_power = evpower.getEigenvalue();
    Eigen::Matrix4f evec_power = evpower.getEigenvector();

    std::cout << "Eigenvalues determined by power method:\n" << eval_power << "\n";
    std::cout << "Eigenvectors determined by power method:\n" << evec_power << "\n\n";
    return 0;
}