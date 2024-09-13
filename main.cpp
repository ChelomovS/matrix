#include <iostream>

#include "matrix.hpp"

int main() {
    long long matrix_degree;
    std::cin >> matrix_degree;
    if ((!std::cin.good()) || (matrix_degree < 0)) {
        std::cerr << "Error input" << std::endl;
        std::abort();
    }

    std::vector<double> elems = {};
    double elem;
    for (auto i = 0; i < matrix_degree * matrix_degree; i++) {
        std::cin >> elem;
        elems.push_back(elem);
    }
    
    Matrix::Matrix<double> matrix(matrix_degree, matrix_degree, elems.begin(), elems.end());
    std::cin >> matrix;
    #if 1
    std::cout << matrix;
    #endif
    std::cout << Matrix::det(matrix) << std::endl;

    return 0;
}