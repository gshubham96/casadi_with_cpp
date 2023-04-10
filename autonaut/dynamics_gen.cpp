#include <iostream>
#include <string>
#include <filesystem>
#include <casadi/casadi.hpp>

namespace fs = std::filesystem;

int main(){

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