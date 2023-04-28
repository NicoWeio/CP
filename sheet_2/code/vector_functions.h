#ifndef VECTOR_FUNCTIONS
#define VECTOR_FUNCTIONS

#include <iostream>
#include <vector>
#include <math.h>

// print vector of any type
template <typename T>
void vec_print(const std::vector<T>& vec) {
    std::cout << "Vector elements: ";
    for (const auto& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

// vector addition
std::vector<double> vec_add(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> result(a.size());

    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] + b[i];
    }
    return result;
}

// vector subtraction
std::vector<double> vec_sub(const std::vector<double>& a, const std::vector<double>& b){
    std::vector<double> result(a.size());

    for(int i = 0; i < a.size(); i++){
        result[i] = a[i] - b[i];
    }
    return result;
}

// absolute of a vector
double vec_abs(const std::vector<double>& a){
    double value = 0;
    for(int i = 0; i<a.size(); i++){
        value += pow(a[i], 2);
    }
    return sqrt(value);
}

// vector times scalar
std::vector<double> vec_scalar(const std::vector<double>& vec, const double a){
    std::vector<double> result(vec.size(), 0);
    for(int i = 0; i < vec.size(); i++){
        result[i] = vec[i]*a;
    }
    return result;
}

#endif // VECTOR_FUNCTIONS
