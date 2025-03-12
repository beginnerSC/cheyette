#include <pybind11/pybind11.h>
#include "piecewisefunction.hpp"
#include "example1.hpp"
#include "example2.hpp"

PYBIND11_MODULE(_ppc, m) {
    m.doc() = "This is _cheyette's docstring.";
    m.def("add", &add, "Add up two numbers.");
    m.def("sub", &sub, "Find difference of two numbers.");
}

void foo(){
    int c, d;
    
    c = add(1, 2);
    d = sub(1, 2);
}

int main()
{
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30}; 
    std::vector<double> zeros = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
    PiecewiseFunction f1(times, zeros), f2(times, zeros), g1(times, zeros), g2(times, zeros), zeroFunction(times, zeros); 

    PiecewiseFunction h = f1*f2 + g1*g2;
    // std::cout << "(f1*f2)(3.5) = " << (f1*f2)(3.5) << std::endl;
    // std::cout << "(f1*f2 + g1*g2)(3.5) = " << (f1*f2 + g1*g2)(3.5) << std::endl;
    std::cout << "h(3.5) = " << h(3.5) << std::endl;

    PiecewiseFunction q(times, zeros);
    q = f1*f2 + g1*g2;
    std::cout << "q(3.5) = " << q(3.5) << std::endl;

    MatrixPiecewiseFunction A(zeroFunction), B(zeroFunction); 
    (A * B).printEvaluated(3.5);
}