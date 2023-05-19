#include <fstream>
#include <iostream>
#include <math.h>

class Poisson2d {
  private:
    int xsize;
    int ysize;
    double delta;
    double **phi; // 2d array: shape(xsize, ysize)
    double **phi_old;
    double **rho;
    double **E_x;
    double **E_y;

    double calc_error();
    void calc_E();

  public:
    Poisson2d(double L_x, double L_y, double delta);
    ~Poisson2d();
    void set_phi(double inner, double left, double right, double bottom, double top);
    void add_charge(double x, double y, double q);
    void gauss_seidel(double delta, double kappa, std::string output_path_base);
};

Poisson2d::Poisson2d(double L_x, double L_y, double delta)
    : delta(delta),
      xsize(int(L_x / delta)),
      ysize(int(L_y / delta)) {
    // create 2d arrays: phi_{n+1}, phi_n, rho
    phi = new double *[xsize];
    phi_old = new double *[xsize];
    rho = new double *[xsize];
    E_x = new double *[xsize];
    E_y = new double *[ysize];

    for (int i = 0; i < xsize; i++) {
        // NOTE: These are suppposed to be initialized using set_phi
        phi[i] = new double[ysize];
        phi_old[i] = new double[ysize];

        rho[i] = new double[xsize];
        // fill rho with zeroes
        for (int j = 0; j < ysize; j++) {
            rho[i][j] = 0.0;
        }

        E_x[i] = new double[xsize];
        E_y[i] = new double[ysize];
        // fill E_x, E_y with zeroes
        for (int j = 0; j < ysize; j++) {
            E_x[i][j] = 0.0;
            E_y[i][j] = 0.0;
        }
    }
}

Poisson2d::~Poisson2d() {
    // clear memory
    for (int i = 0; i < xsize; i++) {
        delete[] phi[i];
        delete[] phi_old[i];
        delete[] rho[i];
    }

    delete[] phi;
    delete[] phi_old;
    delete[] rho;
}

void Poisson2d::set_phi(double inner, double left, double right, double bottom, double top) {
    for (int x = 0; x < xsize; x++) {
        phi[x][0] = bottom;
        phi[x][ysize - 1] = top;
    }
    for (int y = 0; y < ysize; y++) {
        phi[0][y] = left;
        phi[xsize - 1][y] = right;
    }
    for (int x = 1; x < xsize - 1; x++) {
        for (int y = 1; y < ysize - 1; y++) {
            phi[x][y] = inner;
        }
    }
}

void Poisson2d::add_charge(double x, double y, double q) {
    // x and y are in units of L
    int x_index = int(x / delta);
    int y_index = int(y / delta);

    rho[x_index][y_index] += q;
}

void Poisson2d::gauss_seidel(double delta, double kappa, std::string output_path_base) {
    // kappa: stop criteria
    // relative error between to time steps
    double error;

    std::ofstream output_file_phi(output_path_base + "_phi.txt");
    std::ofstream output_file_E(output_path_base + "_E.txt");

    unsigned int iter = 0;
    const unsigned int max_iter = 1000;
    do {
        // copy phi to phi_old
        for (int i = 0; i < xsize; i++) {
            for (int j = 0; j < ysize; j++) {
                phi_old[i][j] = phi[i][j];
            }
        }

        // update phi
        for (int i = 1; i < xsize - 1; i++) {
            for (int j = 1; j < ysize - 1; j++) {
                phi[i][j] = 0.25 * (phi[i + 1][j] + phi[i][j + 1] + phi[i][j - 1] + phi[i - 1][j]) + 0.25 * delta * delta * rho[i][j];
            }
        }

        // update E
        calc_E();

        // calc relative change of euclidean distance between phi and phi_old
        error = calc_error();
        std::cout << "Error: " << error << std::endl;

        // phi write a line to our csv file, e.g.
        // 1,2,3,|4,5,6,|7,8,9,|
        for (int i = 0; i < xsize; i++) {
            for (int j = 0; j < ysize; j++) {
                output_file_phi << phi[i][j] << ",";
                output_file_E << E_x[i][j] << ";" << E_y[i][j] << ",";
            }
            output_file_phi << "|";
            output_file_E << "|";
        }
        output_file_phi << std::endl;
        output_file_E << std::endl;

        if (iter++ > max_iter) {
            std::cout << "Max iterations reached" << std::endl;
            break;
        }

    } while (error > kappa);

    output_file_phi.close();
    output_file_E.close();
}

void Poisson2d::calc_E() {
    for (int i = 1; i < xsize - 1; i++) {
        for (int j = 1; j < ysize - 1; j++) {
            E_x[i][j] = (phi[i + 1][j] - phi[i - 1][j]) / (2 * delta);
            E_y[i][j] = (phi[i][j + 1] - phi[i][j - 1]) / (2 * delta);
        }
    }
}

double Poisson2d::calc_error() {
    // calc rel error between phi and phi_old
    // distance measure: euclidean distance

    double norm2_dif = 0.0;

    for (int i = 0; i < xsize; i++) {
        for (int j = 0; j < ysize; j++) {
            norm2_dif += std::pow((phi[i][j] - phi_old[i][j]), 2);
        }
    }

    return std::sqrt(norm2_dif);
}

int main() {
    // given constants
    double delta = 0.05;
    double kappa = 0.00001;
    double L = 1.0;

    Poisson2d poisson(L, L, delta);
    poisson.set_phi(1.0, 0.0, 0.0, 0.0, 0.0);
    poisson.add_charge(0.5, 0.5, 1.0);
    poisson.gauss_seidel(delta, kappa, "build/A1_b");
    return 0;
}
