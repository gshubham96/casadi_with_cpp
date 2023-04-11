#include <casadi/casadi.hpp>

int main() {
    // Define symbolic variables
    casadi::SX x = casadi::SX::sym("x", 2);
    casadi::SX y = casadi::SX::sym("y", 1);

    // Define function f
    casadi::Function f("f", {x}, {x(0) * x(1)});

    // Define function g
    casadi::Function g("g", {x}, {y});
    g.init();

    // Evaluate g using f
    casadi::SXDict args;
    args["x"] = casadi::DM({1, 2}); // Set input value for x
    g.map(args); // Evaluate g using the input value for x
    casadi::DM res = g.output(); // Get the output value for y

    std::cout << res << std::endl; // Output: [2]

    return 0;
}
