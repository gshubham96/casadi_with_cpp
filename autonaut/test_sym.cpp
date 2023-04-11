#include <casadi/casadi.hpp>

int main() {
    // Define symbolic variables
    casadi::SX x = casadi::SX::sym("x", 2);
    casadi::SX y = casadi::SX::sym("y", 1);

    // Define function f as f = x(0)*x(1)
    casadi::Function f = casadi::Function("f", {x}, {x(0)*x(1)});

    // Define y symbolically as y = f(x)
    y = f({x});

    // Print the symbolic expression of y
    std::cout << "Symbolic expression of y: " << y << std::endl;

    return 0;
}
