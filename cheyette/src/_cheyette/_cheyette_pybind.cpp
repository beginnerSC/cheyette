#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "piecewisefunction.hpp"
#include "example1.hpp"
#include "example2.hpp"
#include <array>
#include <chrono>

double h(double t)
{
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> k = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    PiecewiseFunction k_int(times, k, true);
    PiecewiseFunction h_ = exp(-k_int);

    return h_(t);
}

double h_shifted(double t, double shift)
{
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> k = {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};
    PiecewiseFunction k_int(times, k, true);
    PiecewiseFunction h_ = exp(-k_int);

    return h_.shift(shift)(t);
}

std::vector<std::vector<double>> y_bar(double t);

PYBIND11_MODULE(_cheyette, m) {
    m.doc() = "This is _cheyette's docstring.";
    m.def("add", &add, "Add up two numbers.");
    m.def("sub", &sub, "Find difference of two numbers.");

    m.def("h",          &h);
    m.def("h_shifted",  &h_shifted);
    m.def("y_bar",      &y_bar);
}

int main0(){
    int c, d;
    
    c = add(1, 2);
    d = sub(1, 2);

    return 0;
}

int main1()
{
    std::vector<double> times1 = {5, 13, 24};
    std::vector<double> values1 = {10, 23, 25};

    std::vector<double> times2 = {8, 15, 25};
    std::vector<double> values2 = {12, 24, 26};

    PiecewiseFunction pf1(times1, values1), pf2(times2, values2);

    std::vector<double> points = {4, 7, 12, 14, 23, 24.5};
    for (double p : points) {
        std::cout << (pf1 / pf2)(p) << std::endl;
    }
    return 0;
}

int main3()
{
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> values =   {0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03};

    PiecewiseFunction f(times, values);

    double shift = 0.33; 
    double time = 0.2;
    std::cout << f(time+shift) << ", " << f.shift(shift)(time) << std::endl;

    return 0;
}

int main_pcf_integrate_will_crash()
{
    PiecewiseFunction f, g, h;
    f = PiecewiseFunction({1, 2, 3}, {1, 2, 3}, true);
    g = PiecewiseFunction({1, 2, 3}, {1, 2, 3}).integral();
    h = PiecewiseFunction({1, 2, 3}, {1, 2, 3}) + PiecewiseFunction({1, 2, 3}, {1, 2, 3});

    std::cout << f(3.5) << std::endl;   // crash
    std::cout << g(3.5) << std::endl;   // fine
    std::cout << h(3.5) << std::endl;   // fine 
    
    return 0;
}

std::vector<std::vector<double>> y_bar(double t)
{
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    
    std::vector<std::vector<double>> k =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.016, 0.021, 0.026, 0.028, 0.029, 0.030, 0.031, 0.031, 0.031, 0.031}, 
                                            {0.018, 0.023, 0.028, 0.03, 0.031, 0.032, 0.033, 0.033, 0.033, 0.033}};
    
    std::vector<std::vector<double>> lambda =  {{0.00025   , 0.00033333, 0.00041667, 0.00045   , 0.00046667, 0.00048333, 0.0005    , 0.0005    , 0.0005    , 0.0005    },
                                                {0.0025    , 0.00333333, 0.00416667, 0.0045    , 0.00466667, 0.00483333, 0.005     , 0.005     , 0.005     , 0.005     },
                                                {0.005     , 0.00583333, 0.00666667, 0.0075    , 0.00833333, 0.00916667, 0.01      , 0.01      , 0.01      , 0.01      }};

    std::vector<std::vector<double>> a =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
        
    std::vector<std::vector<double>> b =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
        
    std::vector<std::vector<double>> f =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
    

    std::vector<double> delta = {2, 5, 10};
    std::vector<PiecewiseFunction> k_int(3), h(3), h_inv(3), lambda_(3), a_(3), b_(3), f_(3);
    MatrixPiecewiseFunction H, Lambda, A, B, F, H_f;

    for (size_t i=0 ; i<3 ; ++i)
    {
        // k_int[i] = PiecewiseFunction(times, k[i], true);   // will crash and cannot fix with shared_ptr of this in ctor

        k_int[i] = PiecewiseFunction(times, k[i]).integral();
        h[i] = exp(-k_int[i]);
        h_inv[i] = exp(k_int[i]);
        lambda_[i] = PiecewiseFunction(times, lambda[i]);
        a_[i] = PiecewiseFunction(times, a[i]);
        b_[i] = PiecewiseFunction(times, b[i]);
        f_[i] = PiecewiseFunction(times, f[i]);

        H[i][i] = h[i];
        Lambda[i][i] = lambda_[i];
        A[i][i] = a_[i];
        B[i][i] = b_[i];
        F[i][i] = f_[i];
    }

    MatrixPiecewiseFunction sigma_f_0 = Lambda*(A + B*F);
    std::vector<std::vector<double>> DDT(3, std::vector<double>(3, 0));

    DDT[0][0] = 1.00; DDT[0][1] = 0.85; DDT[0][2] = 0.75;
    DDT[1][0] = 0.85; DDT[1][1] = 1.00; DDT[1][2] = 0.65;
    DDT[2][0] = 0.75; DDT[2][1] = 0.65; DDT[2][2] = 1.00;

    for (size_t i=0 ; i<3 ; ++i){
        for (size_t j=0 ; j<3 ; ++j){
            H_f[i][j] = h[i].shift(delta[j])/h[i];
        }
    }
    PiecewiseFunction k_int_(times, k[0], true);
    PiecewiseFunction h0 = exp(-k_int_);

    return (H*( H.inverse() * H_f.transpose().inverse() * sigma_f_0 * DDT * sigma_f_0 * H_f.inverse() * H.inverse() ).integral() * H).evaluate(t);
}

int main()
{
    return 0;
}

int main_Eigen() {

    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    
    std::vector<std::vector<double>> k =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.016, 0.021, 0.026, 0.028, 0.029, 0.030, 0.031, 0.031, 0.031, 0.031}, 
                                            {0.018, 0.023, 0.028, 0.03, 0.031, 0.032, 0.033, 0.033, 0.033, 0.033}};

    std::vector<std::vector<double>> lambda =  {{0.00025   , 0.00033333, 0.00041667, 0.00045   , 0.00046667, 0.00048333, 0.0005    , 0.0005    , 0.0005    , 0.0005    },
                                                {0.0025    , 0.00333333, 0.00416667, 0.0045    , 0.00466667, 0.00483333, 0.005     , 0.005     , 0.005     , 0.005     },
                                                {0.005     , 0.00583333, 0.00666667, 0.0075    , 0.00833333, 0.00916667, 0.01      , 0.01      , 0.01      , 0.01      }};
    
    std::vector<std::vector<double>> a =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
    
    std::vector<std::vector<double>> b =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
    
    std::vector<std::vector<double>> f =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};

    std::vector<double> delta = {2, 5, 10};
    std::vector<PiecewiseFunction> k_int(3), h(3), h_inv(3), lambda_(3), a_(3), b_(3), f_(3);

    for (size_t i=0 ; i<3 ; ++i)
    {
        // k_int[i] = PiecewiseFunction(times, k[i], true);   // will crash and cannot fix with shared_ptr of this in ctor

        k_int[i] = PiecewiseFunction(times, k[i]).integral();
        h[i] = exp(-k_int[i]);
        h_inv[i] = exp(k_int[i]);
        lambda_[i] = PiecewiseFunction(times, lambda[i]);
        a_[i] = PiecewiseFunction(times, a[i]);
        b_[i] = PiecewiseFunction(times, b[i]);
        f_[i] = PiecewiseFunction(times, f[i]);
    }

    DiagonalMatrixPF Lambda(lambda_), A(a_), B(b_), F(f_), H(h), H_inv(h_inv);
    DiagonalMatrixPF sigma_f_0 = Lambda*(A + B*F);
    EigenMatrixPF H_f;
    
    Eigen::Matrix3d DDT;

    DDT(0, 0) = 1.00; DDT(0, 1) = 0.85; DDT(0, 2) = 0.75;
    DDT(1, 0) = 0.85; DDT(1, 1) = 1.00; DDT(1, 2) = 0.65;
    DDT(2, 0) = 0.75; DDT(2, 1) = 0.65; DDT(2, 2) = 1.00;

    for (size_t i=0 ; i<3 ; ++i){
        for (size_t j=0 ; j<3 ; ++j){
            H_f(i, j) = h[i].shift(delta[j])/h[i];
        }
    }

    auto start = std::chrono::high_resolution_clock::now();

    (H*( H_inv * H_f.transpose().inverse() * sigma_f_0 * DDT * sigma_f_0 * H_f.inverse() * H_inv ).integral() * H).printEvaluated(25);     // integral is slow, otherwise it's fast

    // Code to be timed
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}

int main_vec_vec() {
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    
    std::vector<std::vector<double>> k =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.016, 0.021, 0.026, 0.028, 0.029, 0.030, 0.031, 0.031, 0.031, 0.031}, 
                                            {0.018, 0.023, 0.028, 0.03, 0.031, 0.032, 0.033, 0.033, 0.033, 0.033}};

    std::vector<std::vector<double>> lambda =  {{0.00025   , 0.00033333, 0.00041667, 0.00045   , 0.00046667, 0.00048333, 0.0005    , 0.0005    , 0.0005    , 0.0005    },
                                                {0.0025    , 0.00333333, 0.00416667, 0.0045    , 0.00466667, 0.00483333, 0.005     , 0.005     , 0.005     , 0.005     },
                                                {0.005     , 0.00583333, 0.00666667, 0.0075    , 0.00833333, 0.00916667, 0.01      , 0.01      , 0.01      , 0.01      }};
    
    std::vector<std::vector<double>> a =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
    
    std::vector<std::vector<double>> b =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};
    
    std::vector<std::vector<double>> f =   {{0.015, 0.02, 0.025, 0.027, 0.028, 0.029, 0.03, 0.03, 0.03, 0.03}, 
                                            {0.02, 0.03, 0.035, 0.037, 0.038, 0.039, 0.04, 0.04, 0.04, 0.04}, 
                                            {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.06, 0.06, 0.06}};

    std::vector<double> delta = {2, 5, 10};
    std::vector<PiecewiseFunction> k_int(3), h(3), h_inv(3), lambda_(3), a_(3), b_(3), f_(3);
    MatrixPiecewiseFunction H, H_inv, Lambda, A, B, F, H_f;

    for (size_t i=0 ; i<3 ; ++i)
    {
        // k_int[i] = PiecewiseFunction(times, k[i], true);   // will crash and cannot fix with shared_ptr of this in ctor

        k_int[i] = PiecewiseFunction(times, k[i]).integral();
        h[i] = exp(-k_int[i]);
        h_inv[i] = exp(k_int[i]);
        lambda_[i] = PiecewiseFunction(times, lambda[i]);
        a_[i] = PiecewiseFunction(times, a[i]);
        b_[i] = PiecewiseFunction(times, b[i]);
        f_[i] = PiecewiseFunction(times, f[i]);

        H[i][i] = h[i];
        H_inv[i][i] = h_inv[i];
        Lambda[i][i] = lambda_[i];
        A[i][i] = a_[i];
        B[i][i] = b_[i];
        F[i][i] = f_[i];
    }

    MatrixPiecewiseFunction sigma_f_0 = Lambda*(A + B*F);
    std::vector<std::vector<double>> DDT(3, std::vector<double>(3, 0));

    DDT[0][0] = 1.00; DDT[0][1] = 0.85; DDT[0][2] = 0.75;
    DDT[1][0] = 0.85; DDT[1][1] = 1.00; DDT[1][2] = 0.65;
    DDT[2][0] = 0.75; DDT[2][1] = 0.65; DDT[2][2] = 1.00;

    for (size_t i=0 ; i<3 ; ++i){
        for (size_t j=0 ; j<3 ; ++j){
            H_f[i][j] = h[i].shift(delta[j])/h[i];
        }
    }

    auto start = std::chrono::high_resolution_clock::now();

    (H*( H_inv * H_f.transpose().inverse() * sigma_f_0 * DDT * sigma_f_0 * H_f.inverse() * H_inv ).integral() * H).printEvaluated(25);     // integral is slow, otherwise it's fast

    // Code to be timed
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Execution time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}