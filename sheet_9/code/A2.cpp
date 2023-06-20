#include <iostream>
#include <fstream>
#include <functional>
#include <Eigen/Dense>
#include <random>
#include <math.h>

class ParticleSwarm
{
private:
    // variables
    std::function<double(double, double)> f;
    double r_min, r_max;
    double w, c1, c2;
    double r1,r2;
    int n_particles, n_iterations;

    // position and velocity
    Eigen::MatrixXd r, v, r_best;
    Eigen::VectorXd r_global_best;

    // init random number generator for r1,r2
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

    // functions
    void init();
    void global_best();
    void personal_best();
    void update();

public:
    ParticleSwarm(std::function<double(double, double)> f, double r_min, double r_max, int n_particles, int n_iterations, double w, double c1, double c2);
    Eigen::VectorXd run(std::string path_out);
};

ParticleSwarm::ParticleSwarm(std::function<double(double, double)> f, double r_min, double r_max, int n_particles, int n_iterations, double w, double c1, double c2){
    this->f = f;
    this->r_min = r_min;
    this->r_max = r_max;
    this->n_particles = n_particles;
    this->n_iterations = n_iterations;
    this->w = w;
    this->c1 = c1;
    this->c2 = c2;

    // init position and velocity
    init();

    // init random number generator for r1,r2
    gen = std::mt19937(rd());
    dis = std::uniform_real_distribution<>(0, 1);
}

void ParticleSwarm::init(){
    // init position: scale random matrix [0,1] to [r_min, r_max]
    r = Eigen::MatrixXd::Random(n_particles, 2);
    r = r_min*Eigen::MatrixXd::Constant(n_particles, 2, 1) + (r_max - r_min) * r;

    // init velocity [-1, 1]
    v = Eigen::MatrixXd::Random(n_particles, 2);
    v = -1*Eigen::MatrixXd::Constant(n_particles, 2, 1) + 2 * v;

    // init best position for each particle and global best
    r_best = r;
    r_global_best = r.row(0);
    global_best();
}

void ParticleSwarm::global_best(){
    for (int i = 0; i < n_particles; i++){
        if (f(r(i,0), r(i,1)) < f(r_global_best(0), r_global_best(1))){
            r_global_best = r.row(i);
        }
    }
}

void ParticleSwarm::personal_best(){
    for (int i = 0; i < n_particles; i++){
        if (f(r(i,0), r(i,1)) < f(r_best(i,0), r_best(i,1))){
            r_best.row(i) = r.row(i);
        }
    }
}

void ParticleSwarm::update(){
    // draw random numbers r1, r2
    r1 = dis(gen);
    r2 = dis(gen);

    // update position
    r = r + v;

    // update best position
    personal_best();
    global_best();

    // update velocity
    for (int i = 0; i < n_particles; i++){
        v.row(i) = w * v.row(i) + c1 * r1 * (r_best.row(i) - r.row(i)) + c2 * r2 * (r_global_best.transpose() - r.row(i));
    }
}

Eigen::VectorXd ParticleSwarm::run(std::string path_out){
    // open file
    std::ofstream file(path_out);
    file << "step, x, y, f(x,y)" << std::endl;
    file << 0 << ", " << r_global_best(0) << ", " << r_global_best(1) << ", " << f(r_global_best(0), r_global_best(1)) << std::endl;
    // run algorithm
    for (int i = 1; i < n_iterations; i++){
        update();
        file << i << ", " << r_global_best(0) << ", " << r_global_best(1) << ", " << f(r_global_best(0), r_global_best(1)) << std::endl;
    }

    // close file
    file.close();

    // return global best
    return r_global_best;
}

double f(double x, double y){
    return pow(x-1.9, 2) + pow(y-2.1, 2) + 2*cos(4*x+2.8) + 3*sin(2*y+0.6);
}

int main(){
    // variables
    double r_min = -5.0;
    double r_max = 5.0;
    double w = 0.8;
    double c1 = 0.1;
    double c2 = 0.1;

    int Niter = 1000;
    int Nparticles = 10;

    // minimize f(x,y)
    ParticleSwarm PS(f, r_min, r_max, Nparticles, Niter, w, c1, c2);
    Eigen::VectorXd r = PS.run("build/A2.csv");
    std::cout << "x = " << r(0) << ", y = " << r(1) << ", f(x,y) = " << f(r(0), r(1)) << std::endl;
    return 0;
}