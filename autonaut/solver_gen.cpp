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

    // Optimization variables
    casadi::SX X = casadi::SX::sym("X", 4, 91);
    casadi::SX U = casadi::SX::sym("U", 1, 90);

    std::cout << X << "  :  " <<vertsplit(X) << std::endl;
    return 0;

    // Parameters
    // casadi::SX sym_p = casadi::SX::sym("sym_p", 11);

    // // Objective
    // std::vector<casadi::DM> arg_1 = {X, U, sym_p};
    // casadi::MX f = fun_obj(arg_1);

    // // Constraints
    // casadi::MX g = fun_eql(arg_1);

    // // Create an NLP solver instance
    // casadi::Function solver = casadi::nlpsol("solver", "ipopt", {{"x", X}, {"f", f}, {"g", g}});

    // // file name
    // std::string file_name = "nlp_code";
    // // code predix
    // std::string prefix_code = fs::current_path().string() + "/";

    // // Generate C code for the NLP functions
    // solver.generate_dependencies(file_name + ".c");

    // // shared library prefix
    // std::string prefix_lib = fs::current_path().string() + "/";

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