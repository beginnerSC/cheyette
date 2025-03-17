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
    
    std::vector<double> k1 = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    std::vector<double> k2 = {0.15, 0.2, 0.25, 0.27, 0.28, 0.29, 0.3, 0.3, 0.3, 0.3}; 
    std::vector<double> k3 = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.6, 0.6, 0.6};
    
    std::vector<double> lambda1 = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    std::vector<double> lambda2 = {0.15, 0.2, 0.25, 0.27, 0.28, 0.29, 0.3, 0.3, 0.3, 0.3}; 
    std::vector<double> lambda3 = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.6, 0.6, 0.6};
    
    std::vector<double> a1 = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    std::vector<double> a2 = {0.15, 0.2, 0.25, 0.27, 0.28, 0.29, 0.3, 0.3, 0.3, 0.3}; 
    std::vector<double> a3 = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.6, 0.6, 0.6};
    
    std::vector<double> b1 = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    std::vector<double> b2 = {0.15, 0.2, 0.25, 0.27, 0.28, 0.29, 0.3, 0.3, 0.3, 0.3}; 
    std::vector<double> b3 = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.6, 0.6, 0.6};
    
    std::vector<double> f1 = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    std::vector<double> f2 = {0.15, 0.2, 0.25, 0.27, 0.28, 0.29, 0.3, 0.3, 0.3, 0.3}; 
    std::vector<double> f3 = {0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.6, 0.6, 0.6};
    
    PiecewiseFunction zeroFunction(times, zeros);
    
    PiecewiseFunction k1_int(times, k1, true);
    PiecewiseFunction k2_int(times, k1, true);
    PiecewiseFunction k3_int(times, k1, true);
    
    PiecewiseFunction h1 = exp(-k1_int);
    PiecewiseFunction h2 = exp(-k2_int);
    PiecewiseFunction h3 = exp(-k3_int);
    
    PiecewiseFunction h1_inv = exp(k1_int);
    PiecewiseFunction h2_inv = exp(k2_int);
    PiecewiseFunction h3_inv = exp(k3_int);
    
    PiecewiseFunction lambda1_(times, lambda1);
    PiecewiseFunction lambda2_(times, lambda2);
    PiecewiseFunction lambda3_(times, lambda3);
    
    PiecewiseFunction a1_(times, a1);
    PiecewiseFunction a2_(times, a2);
    PiecewiseFunction a3_(times, a3);
    
    PiecewiseFunction b1_(times, b1);
    PiecewiseFunction b2_(times, b2);
    PiecewiseFunction b3_(times, b3);
    
    PiecewiseFunction f1_(times, f1);
    PiecewiseFunction f2_(times, f2);
    PiecewiseFunction f3_(times, f3);

    MatrixPiecewiseFunction H(zeroFunction), Lambda(zeroFunction), A(zeroFunction), B(zeroFunction), F(zeroFunction);

    std::vector<std::vector<double>> DDT(3, std::vector<double>(3, 0));

    DDT[0][0] = 1.00; DDT[0][1] = 0.85; DDT[0][2] = 0.75;
    DDT[1][0] = 0.85; DDT[1][1] = 1.00; DDT[1][2] = 0.65;
    DDT[2][0] = 0.75; DDT[2][1] = 0.65; DDT[2][2] = 1.00;
    
    H[0][0] = h1;
    H[1][1] = h2;
    H[2][2] = h3;

    A[0][0] = a1_;
    A[1][1] = a2_;
    A[2][2] = a3_;

    B[0][0] = b1_;
    B[1][1] = b2_;
    B[2][2] = b3_;

    F[0][0] = f1_;
    F[1][1] = f2_;
    F[2][2] = f3_;

    Lambda[0][0] = lambda1_;
    Lambda[1][1] = lambda2_;
    Lambda[2][2] = lambda3_;

    MatrixPiecewiseFunction sigma_f_0 = Lambda*(A + B*F);
    MatrixPiecewiseFunction H_f = H;

    (H*( H.inverse() * H_f.transpose().inverse() * sigma_f_0 * DDT * sigma_f_0 * H_f.inverse() * H.inverse() ).integral() * H).printEvaluated(3.5);

    return 0;
}