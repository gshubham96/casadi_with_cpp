#include <casadi/casadi.hpp>

int main() {
    // Define variables
    casadi::SX x = casadi::SX::sym("x", 2);
    casadi::SX y = casadi::SX::sym("y", 1);

    // Define function f
    casadi::SX f = x(0) * x(1);

    // Define function y
    y = f;

    // Create function
    casadi::Function fun("fun", {x}, {y});

    // Print function
    std::cout << fun << std::endl;

    return 0;
}
