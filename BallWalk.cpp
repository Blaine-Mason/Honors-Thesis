#include <iostream>
#include "matplotlibcpp.h"
#include <Eigen/Dense> 

using namespace Eigen;

namespace plt = matplotlibcpp;

int main(){
    VectorXf* all_pos;
    all_pos = new VectorXf[1];
    
    VectorXf X(2);
    X << 0, 0;
    int calls = 0;

    Matrix2Xf pos;
    pos.resize(2,1);
    pos << 0, 0;

    VectorXf Test(2);
    Test << 1, 2;

    pos.resize(2,2);
    all_pos[0] = X;
    all_pos[1] = Test;
    all_pos[2] = Test;
    pos << X, Test;
    for(int i = 0; i < 3; i++){
        std::cout << all_pos[i];
    }
    std::cout << X << std::endl;
    std::cout << pos << std::endl;

    delete[] all_pos;
    
}