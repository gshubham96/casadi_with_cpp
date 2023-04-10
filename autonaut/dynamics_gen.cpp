#include <iostream>
#include <string>
#include <filesystem>
#include <casadi/casadi.hpp>

namespace fs = std::filesystem;

int main(){

    // location for c code generated from matlab and output lib
    std::string file_name = "gen.c -o ";
    std::string prefix_code = fs::current_path().parent_path().string() + "/autonaut/matlab_gen/";
    std::string prefix_lib = fs::current_path().parent_path().string() + "/build/";
    std::string lib_full_name = prefix_lib + "lib_autonaut.so";

    // compile c code and add it to a shared library
    std::string compile_command = "gcc -fPIC -shared -O3 " + 
        prefix_code + file_name + lib_full_name;

    std::cout << compile_command << std::endl;

    int compile_flag = std::system(compile_command.c_str());
    casadi_assert(compile_flag==0, "Compilation failed!");
    std::cout << "Compilation successed!" << std::endl;

    // Load Casadi-dynamics function
    casadi::Function x_dot = casadi::external("x_dot", lib_full_name);

    casadi::SX x = casadi::SX::sym("x");
    casadi::SX y = casadi::SX::sym("y");
    casadi::Function test = casadi::Function("test",{x,y},{x+y});

    // set initial state
    std::vector<double> x0(4, 0);
    x0[1] = 0.9;

    // set input
    std::vector<double> u(1, 0);
    u[0] = 0.0;

    // set params
    std::vector<double> p(11, 0);
    p[9] = 0.9;

    std::cout << "x0 = " << x0 << std::endl;
    std::cout << "u  = " << x0 << std::endl;
    std::cout << "p  = " << x0 << std::endl;

    std::vector<casadi::DM> arg_1 = {casadi::DM(x0), casadi::DM(u), casadi::DM(p)};

    // std::cout << "f_test  = " << test({0.2, 1}) << std::endl;
    std::cout << "x_dot  = " << x_dot(arg_1) << std::endl;


    return 0;
}