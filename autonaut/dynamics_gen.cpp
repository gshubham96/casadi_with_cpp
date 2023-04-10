#include <iostream>
#include <string>
#include <filesystem>
#include <casadi/casadi.hpp>

namespace fs = std::filesystem;

int main(){

    // symbolic variables
    casadi::SX X = casadi::SX::sym("X", 4, 91);
    casadi::SX U = casadi::SX::sym("U", 1, 90);
    casadi::SX sym_p = casadi::SX::sym("sym_p", 11);

    //! functions to be generated 
    casadi::Function f_fun_obj("fun_obj", {posi, posi_des}, {fun_obj});
    casadi::Function f_fun_obj_grad("fun_obj_grad", {posi, posi_des}, {fun_obj_grad});

    // location for c code
    std::string prefix_code = fs::current_path().parent_path().string() + "/autonaut/matlab_gen/";

    // compile c code to a shared library
    std::string prefix_lib = fs::current_path().parent_path().string() + "/build/";
    std::string compile_command = "gcc -fPIC -shared -O3 " + 
        prefix_code + "gen.c -o " +
        prefix_lib + "lib_autonaut.so";

    std::cout << compile_command << std::endl;

    int compile_flag = std::system(compile_command.c_str());
    casadi_assert(compile_flag==0, "Compilation failed!");
    std::cout << "Compilation successed!" << std::endl;

    return 0;
}