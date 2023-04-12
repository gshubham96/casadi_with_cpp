#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <casadi/casadi.hpp>

#include <random>
#include <algorithm>

namespace fs = std::filesystem;

class MyCallback : public casadi::Callback {
 public:
   // Constructor
   MyCallback(void) {
        std::string file_name = "gen.c -o ";
        std::string prefix_code = fs::current_path().parent_path().string() + "/autonaut/matlab_gen/";
        std::string prefix_lib = fs::current_path().parent_path().string() + "/build/";
        std::string lib_full_name = prefix_lib + "lib_autonaut.so";

        // Load Casadi-dynamics function
        x_dot = casadi::external("obj_ms", lib_full_name);
        construct("x_dot");
   }

   // Destructor
   ~MyCallback() override { std::cout << "MyCallback is destroyed here." << std::endl; };

   // Initialize the object
   void init() override {
     std::cout << "initializing object" << std::endl;
   }

   // Number of inputs and outputs
   casadi_int get_n_in() override { return 465;}
   casadi_int get_n_out() override { return 1;}

   // Evaluate numerically
   std::vector<casadi::DM> eval(const std::vector<casadi::DM>& arg) const override {
     casadi::DM x = arg.at(0);
     std::cout << x << std::endl;
     casadi::DM f = x_dot(arg);
     return {f};
   }

 private:
   // Data members
    casadi::Function x_dot;
 };

int main(){

    MyCallback cb ;

    casadi::SX sym_x = casadi::SX::sym("i0", 465);

    casadi::Function f = cb;
    // std::cout << f << std::endl;
    casadi::Function g = casadi::Function("g", {sym_x}, {f(sym_x)});
    // std::cout << g(sym_x) << std::endl;

    return 0;
}