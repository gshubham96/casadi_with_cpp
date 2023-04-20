#include <CourseController.h>

int main(){

    // instantiate a controller with default values
    NMPC::CourseController nmpc;

    // set default parameters
    std::map<std::string, double> params_d;
    params_d["Vc"] = 0.35;      params_d["beta_c"] = 1.57;
    params_d["Vw"] = 5;         params_d["beta_w"] = 1.57;
    params_d["k_1"] = 0.9551;   params_d["k_2"] = -0.031775;
    params_d["Q"] = 4.5;        params_d["R"] = 1.5;
    nmpc.updateMpcParams(params_d);

    // update MPC state
    std::map<std::string, double> state_d;
    state_d["psi"] = 0.091855;
    state_d["u"] = 0.9821;
    state_d["v"] = 0.19964;
    state_d["r"] = 0.031876;
    nmpc.updateMpcState(state_d);

    // update MPC reference
    double chi = 0.53;
    nmpc.updateMpcReference(chi);

    // solve the optimization problem
    if(!nmpc.optimizeMpcProblem())
        std::cerr << "optimization FAILED :(" << std::endl;
    else
        std::cout << "Optimal input is: " << nmpc.getOptimalInput() << std::endl;
        
    nmpc.saveTrajectoryToFile();
    nmpc.print_details();

    return 0;
}