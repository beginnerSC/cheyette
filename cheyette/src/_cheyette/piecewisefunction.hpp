#include <functional>
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <utility>
#include <limits>


class Subfunction {
public:
    Subfunction() : m_func([](double t) { return 0; }) {}
    Subfunction(std::function<double(double)> func) : m_func(func) {}
    Subfunction(const Subfunction& other) : m_func(other.m_func) {}
    double operator()(double t) const {
        return m_func(t);
    }
    Subfunction operator+(const Subfunction& other) const {
        auto shared_this = std::make_shared<Subfunction>(*this);
        auto shared_other = std::make_shared<Subfunction>(other);
        return Subfunction([=](double t) { return (shared_this->m_func)(t) + (shared_other->m_func)(t); });
    }
    Subfunction operator-(const Subfunction& other) const {
        auto shared_this = std::make_shared<Subfunction>(*this);
        auto shared_other = std::make_shared<Subfunction>(other);
        return Subfunction([=](double t) { return (shared_this->m_func)(t) - (shared_other->m_func)(t); });
    }
    Subfunction operator*(const Subfunction& other) const {
        auto shared_this = std::make_shared<Subfunction>(*this);
        auto shared_other = std::make_shared<Subfunction>(other);
        return Subfunction([=](double t) { return (shared_this->m_func)(t) * (shared_other->m_func)(t); });
    }
    Subfunction operator/(const Subfunction& other) const {
        auto shared_this = std::make_shared<Subfunction>(*this);
        auto shared_other = std::make_shared<Subfunction>(other);
        return Subfunction([=](double t) { return (shared_this->m_func)(t) / (shared_other->m_func)(t); });
    }
    Subfunction operator-() const {
        auto shared_this = std::make_shared<Subfunction>(*this);
        return Subfunction([=](double t) { return -(shared_this->m_func)(t); });
    }
    Subfunction& operator=(const Subfunction& other) {
        if (this != &other) {
            m_func = other.m_func;
        }
        return *this;
    }
    Subfunction& operator=(Subfunction&& other) noexcept {
        if (this != &other) {
            m_func = std::move(other.m_func);
        }
        return *this;
    }
    friend Subfunction operator*(double scalar, const Subfunction& func) {
        auto shared_func = std::make_shared<Subfunction>(func);
        return Subfunction([=](double t) { return scalar * (shared_func->m_func)(t); });
    }
    friend Subfunction operator*(const Subfunction& func, double scalar) {
        auto shared_func = std::make_shared<Subfunction>(func);
        return Subfunction([=](double t) { return (shared_func->m_func)(t) * scalar; });
    }
    Subfunction shift(double shift) const {   // f.shift(shift) returns the shifted function f(x + shift)
        auto shared_func = std::make_shared<Subfunction>(*this);
        return Subfunction([=](double t) { return (shared_func->m_func)(t + shift); });
    }
    // Integral function using 10-point Gaussian quadrature
    Subfunction integral(double a) const {
        auto shared_func = std::make_shared<Subfunction>(*this);
        return Subfunction([=](double t) {
            // 10-point Gaussian quadrature weights and abscissae
            static const double weights[10] = {
                0.29552422471475287, 0.29552422471475287,
                0.26926671930999635, 0.26926671930999635,
                0.21908636251598204, 0.21908636251598204,
                0.14945134915058059, 0.14945134915058059,
                0.06667134430868814, 0.06667134430868814
            };
            static const double abscissae[10] = {
                -0.14887433898163122, 0.14887433898163122,
                -0.4333953941292472, 0.4333953941292472,
                -0.6794095682990244, 0.6794095682990244,
                -0.8650633666889845, 0.8650633666889845,
                -0.9739065285171717, 0.9739065285171717
            };

            double integral_value = 0.0;
            double midpoint = 0.5 * (a + t);
            double half_length = 0.5 * (t - a);

            for (int i = 0; i < 10; ++i) {
                double x = midpoint + half_length * abscissae[i];
                integral_value += weights[i] * (shared_func->m_func)(x);
            }

            integral_value *= half_length;
            return integral_value;
        });
    }

private:
    std::function<double(double)> m_func;
};


class PiecewiseFunction {
public:
    PiecewiseFunction() {  // zero function by default
        m_times = {200};
        m_functions.resize(m_times.size(), Subfunction([](double) { return 0.0; }));
    }
    PiecewiseFunction(const std::vector<double>& times, const std::vector<Subfunction>& functions)
        : m_times(times), m_functions(functions) {
            m_times.back() = 200;  // ensure the last time is 200y
        }
    PiecewiseFunction(const PiecewiseFunction& other) 
        : m_times(other.m_times), m_functions(other.m_functions) {
            m_times.back() = 200;  // ensure the last time is 200y
        }
    PiecewiseFunction(PiecewiseFunction&& other) 
        : m_times(std::move(other.m_times)), m_functions(std::move(other.m_functions)) {
            m_times.back() = 200;  // ensure the last time is 200y
        }

    // convient constructor for piecewise constant function or integral of a piecewise constant function
    PiecewiseFunction(const std::vector<double>& times, const std::vector<double>& constants, bool integrate = false) {
        m_times = times;
        m_times.back() = 200;
    
        if (integrate) {
            // Compute the integral of the piecewise constant function
            std::vector<Subfunction> integral_functions;

            for (size_t i = 0; i < m_times.size(); ++i) {
                Subfunction integral_func([=](double t) {
                    size_t q_t = std::upper_bound(m_times.begin(), m_times.end(), t) - m_times.begin() + 1;
                    double integral_value = 0.0;

                    // Add the sum of the intervals [T_j, T_{j+1}] for j from 0 to q(t)-1
                    for (size_t j = 0; j < q_t ; ++j) {
                        integral_value += constants[j] * (m_times[j] - (j==0 ? 0 : m_times[j-1]));
                    }
                    // Add the remaining term -g_{q(t)-1} * (T_{q(t)} - t)
                    integral_value -= constants[q_t - 1] * (m_times[q_t - 1] - t);
                    
                    return integral_value;
                });

                integral_functions.push_back(integral_func);
            }
            m_functions = integral_functions;
        } else {
            // Initialize the piecewise constant functions
            for (const auto& constant : constants) {
                m_functions.emplace_back(Subfunction(std::function<double(double)>([constant](double) { return constant; })));
            }
        }
    }

    double operator()(double t) const {
        for (size_t i = 0; i < m_times.size(); ++i) {
            if (t < m_times[i]) {
                return m_functions[i](t);
            }
        }
        return m_functions.back()(t);  // use the last subfunction if t > T_N
    }

    PiecewiseFunction operator+(const PiecewiseFunction& other) const {
        std::vector<double> new_times = {};
        std::vector<Subfunction> new_functions = {};
        
        int i=0, j=0;

        while (i<m_times.size() && j<other.m_times.size())
        {
            if (m_times[i] < other.m_times[j]) {
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] + other.m_functions[j]);
                i++;
            } else if (m_times[i] > other.m_times[j]) {
                new_times.push_back(other.m_times[j]);
                new_functions.push_back(m_functions[i] + other.m_functions[j]);
                j++;
            } else {    // this and other have the same time, only push back the time once
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] + other.m_functions[j]);
                i++;
                j++;
            }
        }
        if (i==m_times.size()){
            for (size_t k=j ; k<other.m_times.size() ; ++k){
                new_times.push_back(other.m_times[k]);
                new_functions.push_back(m_functions[i-1] + other.m_functions[k]);
            }
        } else if (j==other.m_times.size()){
            for (size_t k=i ; k<m_times.size() ; ++k){
                new_times.push_back(m_times[k]);
                new_functions.push_back(m_functions[k] + other.m_functions[j-1]);
            }
        }
        return PiecewiseFunction(new_times, new_functions);
    }

    PiecewiseFunction& operator+=(const PiecewiseFunction& other) {
        std::vector<double> new_times = {};
        std::vector<Subfunction> new_functions = {};
        
        int i=0, j=0;

        while (i<m_times.size() && j<other.m_times.size())
        {
            if (m_times[i] < other.m_times[j]) {
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] + other.m_functions[j]);
                i++;
            } else if (m_times[i] > other.m_times[j]) {
                new_times.push_back(other.m_times[j]);
                new_functions.push_back(m_functions[i] + other.m_functions[j]);
                j++;
            } else {    // this and other have the same time, only push back the time once
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] + other.m_functions[j]);
                i++;
                j++;
            }
        }
        if (i==m_times.size()){
            for (size_t k=j ; k<other.m_times.size() ; ++k){
                new_times.push_back(other.m_times[k]);
                new_functions.push_back(m_functions[i-1] + other.m_functions[k]);
            }
        } else if (j==other.m_times.size()){
            for (size_t k=i ; k<m_times.size() ; ++k){
                new_times.push_back(m_times[k]);
                new_functions.push_back(m_functions[k] + other.m_functions[j-1]);
            }
        }
        m_times = new_times;
        m_functions = new_functions;
        return *this;
    }

    PiecewiseFunction operator-(const PiecewiseFunction& other) const {
        std::vector<double> new_times = {};
        std::vector<Subfunction> new_functions = {};
        
        int i=0, j=0;

        while (i<m_times.size() && j<other.m_times.size())
        {
            if (m_times[i] < other.m_times[j]) {
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] - other.m_functions[j]);
                i++;
            } else if (m_times[i] > other.m_times[j]) {
                new_times.push_back(other.m_times[j]);
                new_functions.push_back(m_functions[i] - other.m_functions[j]);
                j++;
            } else {    // this and other have the same time, only push back the time once
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] - other.m_functions[j]);
                i++;
                j++;
            }
        }
        if (i==m_times.size()){
            for (size_t k=j ; k<other.m_times.size() ; ++k){
                new_times.push_back(other.m_times[k]);
                new_functions.push_back(m_functions[i-1] - other.m_functions[k]);
            }
        } else if (j==other.m_times.size()){
            for (size_t k=i ; k<m_times.size() ; ++k){
                new_times.push_back(m_times[k]);
                new_functions.push_back(m_functions[k] - other.m_functions[j-1]);
            }
        }
        return PiecewiseFunction(new_times, new_functions);
    }

    PiecewiseFunction operator*(const PiecewiseFunction& other) const {
        std::vector<double> new_times = {};
        std::vector<Subfunction> new_functions = {};
        
        int i=0, j=0;

        while (i<m_times.size() && j<other.m_times.size())
        {
            if (m_times[i] < other.m_times[j]) {
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] * other.m_functions[j]);
                i++;
            } else if (m_times[i] > other.m_times[j]) {
                new_times.push_back(other.m_times[j]);
                new_functions.push_back(m_functions[i] * other.m_functions[j]);
                j++;
            } else {    // this and other have the same time, only push back the time once
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] * other.m_functions[j]);
                i++;
                j++;
            }
        }
        if (i==m_times.size()){
            for (size_t k=j ; k<other.m_times.size() ; ++k){
                new_times.push_back(other.m_times[k]);
                new_functions.push_back(m_functions[i-1] * other.m_functions[k]);
            }
        } else if (j==other.m_times.size()){
            for (size_t k=i ; k<m_times.size() ; ++k){
                new_times.push_back(m_times[k]);
                new_functions.push_back(m_functions[k] * other.m_functions[j-1]);
            }
        }
        return PiecewiseFunction(new_times, new_functions);
    }

    PiecewiseFunction operator/(const PiecewiseFunction& other) const {
        std::vector<double> new_times = {};
        std::vector<Subfunction> new_functions = {};
        
        int i=0, j=0;

        while (i<m_times.size() && j<other.m_times.size())
        {
            if (m_times[i] < other.m_times[j]) {
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] / other.m_functions[j]);
                i++;
            } else if (m_times[i] > other.m_times[j]) {
                new_times.push_back(other.m_times[j]);
                new_functions.push_back(m_functions[i] / other.m_functions[j]);
                j++;
            } else {    // this and other have the same time, only push back the time once
                new_times.push_back(m_times[i]);
                new_functions.push_back(m_functions[i] / other.m_functions[j]);
                i++;
                j++;
            }
        }
        if (i==m_times.size()){
            for (size_t k=j ; k<other.m_times.size() ; ++k){
                new_times.push_back(other.m_times[k]);
                new_functions.push_back(m_functions[i-1] / other.m_functions[k]);
            }
        } else if (j==other.m_times.size()){
            for (size_t k=i ; k<m_times.size() ; ++k){
                new_times.push_back(m_times[k]);
                new_functions.push_back(m_functions[k] / other.m_functions[j-1]);
            }
        }
        return PiecewiseFunction(new_times, new_functions);
    }

    PiecewiseFunction operator-() const {
        std::vector<Subfunction> new_functions;
        for (const auto& func : m_functions) {
            new_functions.push_back(-func);
        }
        return PiecewiseFunction(m_times, new_functions);
    }

    PiecewiseFunction& operator=(const PiecewiseFunction& other) {
        if (this != &other) {
            m_times = other.m_times;
            m_functions = other.m_functions;
        }
        return *this;
    }

    // PiecewiseFunction& operator=(PiecewiseFunction&& other) {
    //     if (this != &other) {
    //         m_times = std::move(other.m_times);
    //         m_functions = std::move(other.m_functions);
    //     }
    //     return *this;
    // }

    friend PiecewiseFunction operator*(double scalar, const PiecewiseFunction& pf) {
        std::vector<Subfunction> new_functions;
        for (const auto& func : pf.m_functions) {
            new_functions.push_back(scalar * func);
        }
        return PiecewiseFunction(pf.m_times, new_functions);
    }

    friend PiecewiseFunction operator*(const PiecewiseFunction& pf, double scalar) {
        std::vector<Subfunction> new_functions;
        for (const auto& func : pf.m_functions) {
            new_functions.push_back(func * scalar);
        }
        return PiecewiseFunction(pf.m_times, new_functions);
    }

    friend PiecewiseFunction exp(const PiecewiseFunction& pf);

    PiecewiseFunction integral() const {
        std::vector<Subfunction> integral_functions;
        double previous_time = 0.0;
        double accumulated_integral = 0.0;

        for (size_t i = 0; i < m_functions.size(); ++i) {
            Subfunction integral_func = m_functions[i].integral(previous_time);
            integral_functions.push_back(Subfunction([=](double t) {
                return accumulated_integral + integral_func(t);
            }));
            accumulated_integral += integral_func(m_times[i]);
            previous_time = m_times[i];
        }

        return PiecewiseFunction(m_times, integral_functions);
    }
    
    PiecewiseFunction shift(double shift) const {   // f.shift(shift) returns the shifted function f(x + shift)
        std::vector<double> new_times = {};
        std::vector<Subfunction> new_functions = {};
        
        for (size_t i=0 ; i<m_times.size() ; ++i){
            if (m_times[i] > shift){
                new_times.push_back(m_times[i] - shift);
                new_functions.push_back(m_functions[i].shift(shift));
            }
        }
        return PiecewiseFunction(new_times, new_functions);
    }

private:
    std::vector<double> m_times;
    std::vector<Subfunction> m_functions;
};

Subfunction exp(const Subfunction& f) {
    auto shared_f = std::make_shared<Subfunction>(f);
    return Subfunction([=](double t) {
        return std::exp((*shared_f)(t));
    });
}


PiecewiseFunction exp(const PiecewiseFunction& pf) {
    std::vector<Subfunction> new_functions;
    for (const auto& func : pf.m_functions) {
        new_functions.push_back(exp(func));
    }
    return PiecewiseFunction(pf.m_times, new_functions);
}


class MatrixPiecewiseFunction {
public:
    MatrixPiecewiseFunction() 
        : m_matrix(3, std::vector<PiecewiseFunction>(3)) {}     // matrix of zero functions by default

    MatrixPiecewiseFunction(const MatrixPiecewiseFunction& other) 
        : m_matrix(other.m_matrix) {} 

    MatrixPiecewiseFunction(MatrixPiecewiseFunction&& other) 
        : m_matrix(std::move(other.m_matrix)) {} 

    // ctor to fill the matrix with the same piecewise function
    MatrixPiecewiseFunction(const PiecewiseFunction& pf) 
        : m_matrix(3, std::vector<PiecewiseFunction>(3, pf)) {}

    MatrixPiecewiseFunction cofactorMatrix() const {
        MatrixPiecewiseFunction cofactor;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::vector<std::vector<PiecewiseFunction>> minorMatrix(2, std::vector<PiecewiseFunction>(2));
                
                // Fill the minor matrix
                int minorRow = 0;
                for (int row = 0; row < 3; ++row) {
                    if (row == i) continue;
                    int minorCol = 0;
                    for (int col = 0; col < 3; ++col) {
                        if (col == j) continue;
                        minorMatrix[minorRow][minorCol] = m_matrix[row][col];
                        minorCol++;
                    }
                    minorRow++;
                }
                // Compute the determinant of the minor matrix
                PiecewiseFunction minorDet = minorMatrix[0][0] * minorMatrix[1][1] - minorMatrix[0][1] * minorMatrix[1][0];

                // Compute the cofactor
                cofactor.m_matrix[i][j] = ((i + j) % 2 == 0 ? 1 : -1) * minorDet;
            }
        }

        return cofactor;
    }

    PiecewiseFunction determinant() const {
        MatrixPiecewiseFunction cofactor = this->cofactorMatrix();

        // Compute the determinant using the first row
        PiecewiseFunction det = m_matrix[0][0] * cofactor.m_matrix[0][0] +
                                m_matrix[0][1] * cofactor.m_matrix[0][1] +
                                m_matrix[0][2] * cofactor.m_matrix[0][2];
        
        return det;
    }

    MatrixPiecewiseFunction inverse() const {
        PiecewiseFunction det = this->determinant();
        // assuming determinant is not zero)
        MatrixPiecewiseFunction cofactor = this->cofactorMatrix();
        
        // Create a MatrixPiecewiseFunction from the transposed cofactor matrix
        MatrixPiecewiseFunction adjugate = this->transpose();
        MatrixPiecewiseFunction inverse;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                inverse.m_matrix[i][j] = adjugate.m_matrix[i][j] / det;
            }
        }
        return inverse;
    }

    MatrixPiecewiseFunction transpose() const {
        MatrixPiecewiseFunction transposed;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                transposed.m_matrix[i][j] = m_matrix[j][i];
            }
        }
        return transposed;
    }

    MatrixPiecewiseFunction integral() const {
        MatrixPiecewiseFunction result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m_matrix[i][j] = m_matrix[i][j].integral();
            }
        }
        return result;
    }

    MatrixPiecewiseFunction operator+(const MatrixPiecewiseFunction& other) const {
        MatrixPiecewiseFunction result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m_matrix[i][j] = m_matrix[i][j] + other.m_matrix[i][j];
            }
        }
        return result;
    }

    MatrixPiecewiseFunction operator-(const MatrixPiecewiseFunction& other) const {
        MatrixPiecewiseFunction result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m_matrix[i][j] = m_matrix[i][j] - other.m_matrix[i][j];
            }
        }
        return result;
    }

    MatrixPiecewiseFunction operator*(const MatrixPiecewiseFunction& other) const {
        MatrixPiecewiseFunction result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    result.m_matrix[i][j] = result.m_matrix[i][j] + m_matrix[i][k] * other.m_matrix[k][j];
                }
            }
        }
        return result;
    }

    MatrixPiecewiseFunction& operator=(const MatrixPiecewiseFunction& other) {
        if (this != &other) {
            m_matrix = other.m_matrix;
        }
        return *this;
    }

    MatrixPiecewiseFunction& operator=(const MatrixPiecewiseFunction&& other) {
        if (this != &other) {
            m_matrix = std::move(other.m_matrix);
        }
        return *this;
    }

    friend MatrixPiecewiseFunction operator*(const std::vector<std::vector<double>>& c, const MatrixPiecewiseFunction& mpf) {
        MatrixPiecewiseFunction result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    result.m_matrix[i][j] = result.m_matrix[i][j] + c[i][k] * mpf.m_matrix[k][j];
                }
            }
        }
        return result;
    }

    friend MatrixPiecewiseFunction operator*(const MatrixPiecewiseFunction& mpf, const std::vector<std::vector<double>>& c) {
        MatrixPiecewiseFunction result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    result.m_matrix[i][j] = result.m_matrix[i][j] + c[i][k] * mpf.m_matrix[k][j];
                }
            }
        }
        return result;
    }

    std::vector<std::vector<double>> evaluate(const double& t) const {
        std::vector<std::vector<double>> evaluated(3, std::vector<double>(3));
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                evaluated[i][j] = m_matrix[i][j](t);
            }
        }
        return evaluated;
    }

    void printEvaluated(const double& t) const {
        std::vector<std::vector<double>> res = this->evaluate(t);
        for (size_t i = 0; i < 3; ++i) {
            for (size_t j = 0; j < 3; ++j) {
                std::cout << res[i][j] << ", ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    // Overload operator[] to access rows and elements
    std::vector<PiecewiseFunction>& operator[](size_t index) {
        return m_matrix[index];
    }

    const std::vector<PiecewiseFunction>& operator[](size_t index) const {
        return m_matrix[index];
    }

private:
    std::vector<std::vector<PiecewiseFunction>> m_matrix;
};