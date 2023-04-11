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

    casadi::SX sym_x = casadi::SX::sym("i0", 4);
    casadi::SX sym_u = casadi::SX::sym("i1", 1);
    casadi::SX sym_p = casadi::SX::sym("i11", 11);

    // std::vector<casadi::SX> args;
    // std::vector<casadi::SX> args = {sym_x, sym_u, sym_p};
    casadi::SXDict args, res;
    args["i0"] = sym_x;
    args["i1"] = sym_u;
    args["i2"] = sym_p;

    // casadi::SXDict f_eval = xdot(args);

    std::cout << "x_dot = " << x_dot << std::endl;
    std::cout << "x_dot(args) = " << x_dot(args) << std::endl;



    return 0;
}