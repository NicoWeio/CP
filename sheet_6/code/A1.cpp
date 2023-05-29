#include <iostream>
#include <fstream>

class Logmap{
private:
    double x; //current x0
    double r;   

    // helper functions
    void map(int N);
    
public:
    Logmap(double x0, double r);

    void equilibrate(int N);
    double measure(int N, std::string path_out); // same as equi but with write_to_csv()
};

Logmap::Logmap(double x0, double r){
    this-> x = x0;
    this-> r = r;
}

void Logmap::equilibrate(int N){
    std::cout << "Equilibrate..." << std::endl;

    // logistic map
    for (int n_iter = 0; n_iter<N; n_iter++){
        std::cout << n_iter << " / " << N << "\r";
        x = r*x*(1-x);
    }
}

double Logmap::measure(int N, std::string path_out){
    std::cout << "Measure..." << std::endl;
    
    std::ofstream file(path_out);

    // logistic map
    for (int n_iter = 0; n_iter<N; n_iter++){
        file << n_iter << "," << x << "\n";

        std::cout << n_iter+1 << " / " << N << "\r";
        x = r*x*(1-x);
        write_to_csv(file);
    }

    file.close();

    return x; //return last x_n+1
}


int main(){

    return 0;
}