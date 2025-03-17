#include <gtest/gtest.h>
#include "piecewisefunction.hpp"

TEST(PiecewiseFunctionTest, MultiplicationAndAddition) {
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> zeros = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PiecewiseFunction f1(times, zeros), f2(times, zeros), g1(times, zeros), g2(times, zeros);

    PiecewiseFunction h = f1 * f2 + g1 * g2;
    EXPECT_NO_THROW(h(3.5));
}

TEST(PiecewiseFunctionTest, MultiplicationEvaluation) {
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> zeros = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PiecewiseFunction f1(times, zeros), f2(times, zeros);

    EXPECT_NO_THROW((f1 * f2)(3.5));
}

TEST(PiecewiseFunctionTest, AdditionEvaluation) {
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> zeros = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PiecewiseFunction f1(times, zeros), f2(times, zeros), g1(times, zeros), g2(times, zeros);

    EXPECT_NO_THROW((f1 * f2 + g1 * g2)(3.5));
}

TEST(PiecewiseFunctionTest, AssignmentAndEvaluation) {
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> zeros = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PiecewiseFunction f1(times, zeros), f2(times, zeros), g1(times, zeros), g2(times, zeros);

    PiecewiseFunction q(times, zeros);
    q = f1 * f2 + g1 * g2;
    EXPECT_NO_THROW(q(3.5));
}

TEST(MatrixPiecewiseFunctionTest, MultiplicationAndEvaluation) {
    std::vector<double> times = {0.09, 0.25, 0.5, 1, 2, 3, 5, 10, 20, 30};
    std::vector<double> zeros = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    PiecewiseFunction zeroFunction(times, zeros);

    MatrixPiecewiseFunction A(zeroFunction), B(zeroFunction);
    EXPECT_NO_THROW((A * B).printEvaluated(3.5));
}

TEST(SubfunctionTest, IntegralCalculation) {
    Subfunction f([](double t) { return exp(t); });
    Subfunction integral_f = f.integral(1.0);
    double result = integral_f(3.0); // Computes the integral of exp(t) from 1 to 3

    // Expected value of the integral of exp(t) from 1 to 5 is exp(3) - exp(1)
    double expected_result = exp(3.0) - exp(1.0);
    EXPECT_NEAR(result, expected_result, 1e-6);
}

TEST(PiecewiseFunctionTest, IntegralCalculation) {
    std::vector<double> times_ = {1.0, 2.0};

    // Define some non-constant subfunctions for each interval
    std::vector<Subfunction> subfunctions = {
        Subfunction([](double t) { return t; }),          // f(t) = t for 0 <= t < 1
        Subfunction([](double t) { return t * t; }),      // f(t) = t^2 for 1 <= t < 2
    };

    // Create a PiecewiseFunction with the given times and subfunctions
    PiecewiseFunction f(times_, subfunctions);

    // Compute the integral of the PiecewiseFunction from 0 to t
    PiecewiseFunction integral_f = f.integral();

    // Evaluate the integral at different points
    double t1 = 1.0;
    double t2 = 2.0;

    double result1 = integral_f(t1);
    double result2 = integral_f(t2);

    // Expected values
    double expected_result1 = 0.5; // Integral of t from 0 to 1
    double expected_result2 = 0.5 + (1.0/3.0) * (2.0*2.0*2.0 - 1.0*1.0*1.0); // Integral of t from 0 to 1 plus integral of t^2 from 1 to 2

    EXPECT_NEAR(result1, expected_result1, 1e-6);
    EXPECT_NEAR(result2, expected_result2, 1e-6);
}


int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}