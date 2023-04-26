#include <iostream>
#include <vector>
#include "vector_functions.h"
#include <fstream>
#include <cstdlib> //for exit(1) to catch errors
#include <string>
#include <algorithm>

// calc acceleration vector a
std::vector<double> get_a(const std::vector<double>& r1, const std::vector<double>& r2, const double M){
    // a = F(r1, r2)/m
    std::vector<double> r_dif = vec_sub(r1, r2);
    
    //check if r1-r2 == 0
    if(vec_abs(r_dif)==0){
        std::cerr << "Division by zero: r1 and r2 must be different!";
        std::exit(1);
    }

    double factor = -M / pow(vec_abs(r_dif), 3);
    return vec_scalar(r_dif, factor);
}

// calc velocity vector v
std::vector<double> get_v(const std::vector<double>& v_old, const std::vector<double>& a, const double h){
    // v + a*h
    return vec_add(v_old, vec_scalar(a, h));
}

// calc space vector r
std::vector<double> get_r(const std::vector<double>& r_old, const std::vector<double>& v, const double h){
    // r + v*h
    return vec_add(r_old, vec_scalar(v, h));
}

// Euler algorithm
//results (space vectors r1, r2) in csv file
void euler(std::vector<double> r1, std::vector<double> r2, std::vector<double> v1, std::vector<double> v2, double m1, double m2, double h, double T_max){
    // var m1
    std::vector<double> r1_old{0.0, 0.0};
    std::vector<double> v1_old{0.0, 0.0};
    std::vector<double> a1{0.0, 0.0};
    std::vector<double> a1_old{0.0, 0.0};

    // var m2
    std::vector<double> r2_old{0.0, 0.0};
    std::vector<double> v2_old{0.0, 0.0};
    std::vector<double> a2{0.0, 0.0};
    std::vector<double> a2_old{0.0, 0.0};


    // open csv files to save the space vectors r1, r2
    // use stepsize h in filename
    std::string h_string = std::to_string(h);
    std::replace(h_string.begin(), h_string.end(), '.', '_');

    std::string path_r1 = std::string("build/euler_r1_h") + h_string + std::string(".csv");
    std::string path_r2 = std::string("build/euler_r2_h") + h_string + std::string(".csv");

    std::ofstream file_r1(path_r1);
    std::ofstream file_r2(path_r2);

    // algorithm
    for(int i = 1; i*h<T_max; i++){ //t = i*h, i=1,2,..
        //m1
        a1_old = a1;
        a1 = get_a(r1, r2, m2); //update acceleration a
        
        v1_old = v1;
        v1 = get_v(v1_old, a1, h); //calc new velocity
        
        r1_old = r1;
        r1 = get_r(r1_old, v1, h); //calc new position

        file_r1 << r1[0] << ", " << r1[1] << std::endl; //write to csv

        //m2
        a2_old = a2;
        a2 = get_a(r2, r1_old, m1); //update acceleration a
        
        v2_old = v2;
        v2 = get_v(v2_old, a2, h); //calc new velocity
        
        r2_old = r2;
        r2 = get_r(r2_old, v2, h); //calc new position
        
        file_r2 << r2[0] << ", " << r2[1] << std::endl; //write to csv        
    }

    // close csv files
    file_r1.close();
    file_r2.close();
}

// Verlet algorithm
void verlet(std::vector<double> r1, std::vector<double> r2, std::vector<double> v1, std::vector<double> v2, double m1, double m2, double h, double T_max){
    // var m1
    std::vector<double> r1_old{0.0, 0.0}; // r_n
    std::vector<double> r1_oold{0.0, 0.0}; // r_n-1

    std::vector<double> a1{0.0, 0.0};
    std::vector<double> a1_old{0.0, 0.0};

    // var m2
    std::vector<double> r2_old{0.0, 0.0};
    std::vector<double> r2_oold{0.0, 0.0};

    std::vector<double> a2{0.0, 0.0};
    std::vector<double> a2_old{0.0, 0.0}; 

    // calc initial value for r_-1
    // r_-1 = r_0 - v_0*h + 0.5*a_0*h^2
    r1_old = vec_add(vec_sub(r1, vec_scalar(v1,h)), vec_scalar(get_a(r1, r2, m2), 0.5*h*h));
    r2_old = vec_add(vec_sub(r2, vec_scalar(v2,h)), vec_scalar(get_a(r2, r1, m1), 0.5*h*h));

    // open csv files to save the space vectors r1, r2
    // use stepsize h in filename
    std::string h_string = std::to_string(h);
    std::replace(h_string.begin(), h_string.end(), '.', '_');

    std::string path_r1 = std::string("build/verlet_r1_h") + h_string + std::string(".csv");
    std::string path_r2 = std::string("build/verlet_r2_h") + h_string + std::string(".csv");

    std::ofstream file_r1(path_r1);
    std::ofstream file_r2(path_r2);

    // algorithm
    for(int i = 1; i*h<T_max; i++){ //t = i*h, i=1,2,..
        //m1
        a1_old = a1;
        a1 = get_a(r1, r2, m2); // update a

        r1_oold = r1_old;
        r1_old = r1;

        // r_n+1 = 2*r_n - r_n-1 + a_n*h^2
        r1 = vec_add(vec_sub(vec_scalar(r1_old, 2), r1_oold), vec_scalar(a1, h*h));

        file_r1 << r1[0] << ", " << r1[1] << std::endl; //write to csv

        //m2
        a2_old = a2;
        a2 = get_a(r2, r1_old, m1); // update a

        r2_oold = r2_old;
        r2_old = r2;

        // r_n+1 = 2*r_n - r_n-1 + a_n*h^2
        r2 = vec_add(vec_sub(vec_scalar(r2_old, 2), r2_oold), vec_scalar(a2, h*h));

        file_r2 << r2[0] << ", " << r2[1] << std::endl; //write to csv
    }

    // close csv files
    file_r1.close();
    file_r2.close();
}

int main(){
    // Var mass 1: initial conditions
    const double m1 = 1.0;
    const std::vector<double> r1{0.0, 1.0};
    const std::vector<double> v1{0.8, 0.0};

    // Var mass 2: initial conditions
    const double m2 = 2.0;
    const std::vector<double> r2{0.0, -0.5};
    const std::vector<double> v2{-0.4, 0.0};

    // Var time
    // t = i*h = 1*h, 2*h, ..., T_max
    const double T_max = 100.0;

    // Euler algorithm
    euler(r1, r2, v1, v2, m1, m2, 1.0, T_max);
    euler(r1, r2, v1, v2, m1, m2, 0.1, T_max);
    euler(r1, r2, v1, v2, m1, m2, 0.01, T_max);

    // Verlet algorithm
    verlet(r1, r2, v1, v2, m1, m2, 1.0, T_max);
    verlet(r1, r2, v1, v2, m1, m2, 0.1, T_max);
    verlet(r1, r2, v1, v2, m1, m2, 0.01, T_max);


    return 0;
}