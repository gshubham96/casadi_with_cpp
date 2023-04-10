#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <casadi/casadi.hpp>

namespace fs = std::filesystem;

int main(){

    // Load optimization function
    std::string prefix_lib = fs::current_path().parent_path().string() + "/build/";
    std::string lib_full_name = prefix_lib + "lib_autonaut.so";

    // use this function
    casadi::Function fun_obj = casadi::external("obj_ms", lib_full_name);
    casadi::Function fun_eql = casadi::external("eql_ms", lib_full_name);

    // Optimization variables - contains both state and input variables -- [X00 X01 X02 X03 X10 X11 X12 X13 ... U0 U1 U2 ....]
    casadi::SX X = casadi::SX::sym("X", 454);

    // Parameters
    casadi::SX sym_p = casadi::SX::sym("sym_p", 11);

    // input arguements
    std::vector<casadi::SX> arg = {X, sym_p};

    std::cout << "function arg : " << fname_n_in(fun_obj) << std::endl;


    // // Objective
    // casadi::SX f = fun_obj(arg);

    // // Constraints
    // casadi::SX g = fun_eql(arg);

    // std::cout << "function obj : " << f << std::endl;
    // std::cout << "function g   : " << g << std::endl;


    // // Create an NLP solver instance
    // casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", X}, {"f", f}, {"g", g}});

    // // file name
    // std::string file_name = "autonaut_nlp";
    // // code predix
    // std::string prefix_code = fs::current_path().string() + "/code_gen/";

    // // Generate C code for the NLP functions
    // solver.generate_dependencies(file_name + ".c");

    // // compile c code to a shared library
    // std::string compile_command = "gcc -fPIC -shared -O3 " + 
    //     prefix_code + file_name + ".c -o " +
    //     prefix_lib + file_name + ".so";

    // std::cout << compile_command << std::endl;
    // int compile_flag = std::system(compile_command.c_str());
    // casadi_assert(compile_flag==0, "Compilation failed");
    // std::cout << "Compilation successed!" << std::endl;

    return 0;
}