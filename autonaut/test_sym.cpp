#include <casadi/casadi.hpp>

int main() {
  // Define symbolic variables
  casadi::SX x = casadi::SX::sym("x", 2);
  casadi::SX y = casadi::SX::sym("y", 1);

  // Define function f = x[0]*x[1]
  casadi::SX f = x[0]*x[1];
  casadi::Function f_func("f", {x}, {f});

  // Define y = f(x)
  casadi::SXDict args;
  args["x"] = x;
  casadi::SX f_eval = f_func(args);
  y = f_eval;
  std::cout << "y = " << y << std::endl;

  return 0;
}
