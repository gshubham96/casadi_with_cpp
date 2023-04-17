#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <casadi/casadi.hpp>

#include <random>
#include <algorithm>
// floor
#include <cmath>    

#define PI 3.14159
#define EPS 1e-9
#define DEG2RAD(angle) ((angle) * M_PI / 180.0)

namespace fs = std::filesystem;

namespace NMPC{

    class CourseController {

        private:
            // initialization variable
            bool initialized, state_update, config_update;
            // lengths of state, input and paramter vectors
            int nx, nu, np, N;
            // time of last update
            double t_update, Tp, Ts;
            // course reference for the controller
            double reference_;
            // read system parameters from "mat" file and load here
            std::map<std::string, double> system_;
            // config parameters, runtime paramters and state for MPC
            std::map<std::string, double> config_, params_, state_;
            // initial guess for warm start
            std::map<std::string, std::vector<double>> args_;
            // NLP Solver
            casadi::Function solver;
            // optimized input trajectory
            std::vector<double> input_traj_;

            double ssa(double diff) {
                while (diff < -PI) diff += 2 * PI;
                while (diff > PI) diff -= 2 * PI;
                return diff;
            }

            casadi::SX ssa(casadi::SX diff) {
                return diff;
                // return fmod(diff, 2*PI) - PI;
            }

            // Function to define and compile the NLP Optimization Problem
            bool defineMpcProblem(void){

                // assign configuration parameters
                nx = config_["nx"]; nu = config_["nu"]; np = config_["np"];
                Tp = config_["Tp"]; Ts = config_["Ts"];
                double model_dim = config_["model_dim"], model_type = config_["model_type"], cost_type = config_["cost_type"];

                std::cout << model_dim << std::endl;

                N = floor(Tp / Ts);

                // assign system parameters
                double D11 = system_["D11"], R11 = system_["R11"], INV_M11 = system_["INV_M11"];
                double D22 = system_["D22"], R22 = system_["R22"], INV_M22 = system_["INV_M22"], INV_M23 = system_["INV_M23"];
                double D33 = system_["D33"], R33 = system_["R33"], INV_M32 = system_["INV_M32"], INV_M33 = system_["INV_M33"];

                // named symbolica vars
                casadi::SX 
                    psi = casadi::SX::sym("psi", 1),
                    u = casadi::SX::sym("u", 1), 
                    v = casadi::SX::sym("v", 1), 
                    r = casadi::SX::sym("r", 1),
                    delta = casadi::SX::sym("delta", 1);

                // optim vars for a single shooting period 
                casadi::SX 
                    sym_x = vertcat(psi, u, v, r),
                    sym_u = delta,
                    sym_p = casadi::SX::sym("p_x0", np);

                // environmental parameters that are constant over a given horizon
                casadi::SX
                    chi_d = sym_p(nx),
                    Vc = sym_p(nx+1),
                    beta_c = sym_p(nx+2),
                    Vw = sym_p(nx+3),
                    beta_w = sym_p(nx+4),
                    k_1 = sym_p(nx+5),
                    k_2 = sym_p(nx+6),
                    Q = sym_p(nx+7),
                    R = sym_p(nx+8);

                // set initial position
                casadi::SX sym_x0;
                for(int i = 0; i < nx; i++){
                    sym_x0(i) = sym_p(i);
                }

                // detived states
                casadi::SX
                    u_e = u + EPS,
                    u_c = Vc * cos(beta_c - psi),
                    v_c = Vc * sin(beta_c - psi),
                    u_r = u_e - u_c,
                    v_r = v - v_c,
                    U_r2 = pow(u_r, 2) + pow(v_r, 2),
                    beta = atan(v / u_e);

                // ################################################
                // ###----------------DYNAMICS------------------###
                // ################################################

                // CURRENTS
                casadi::SX nu_c_dot_u = v_c * r;
                casadi::SX nu_c_dot_v = -u_c * r;

                // DAMPING
                casadi::SX damping_u  = D11;
                casadi::SX damping_v  = D22;
                casadi::SX damping_r  = D33;

                // YAW
                casadi::SX yaw_dot = r;

                // WAVE FOILS
                casadi::SX tau_foil_u = (k_1 + k_2*cos(psi - beta_w - PI)) * D11;

                // RUDDER
                casadi::SX tau_rudr_u, tau_rudr_v, tau_rudr_r;
                // If the model is nonlinear, consider nonlinear dynamics of the rudder
                if(model_type == 0){
                    casadi::SX alpha_r = delta - atan(v_r/u_r);
                    tau_rudr_u = R11 * U_r2 * sin(alpha_r) * sin(delta) ;
                    tau_rudr_v = R22 * U_r2 * sin(alpha_r) * cos(delta) ;
                    tau_rudr_r = R33 * U_r2 * sin(alpha_r) * cos(delta) ;
                }
                // else consider the approximated linear equations
                else if(model_type == 1){
                    tau_rudr_u = R11 * U_r2 * delta * delta ;
                    tau_rudr_v = R22 * U_r2 * delta * 0.5 ;
                    tau_rudr_r = R33 * U_r2 * delta * 0.5 ;
                }

                // #TODO WIND    

                // #TODO CORIOLIS 

                // dynamics of surge
                casadi::SX u_dot = nu_c_dot_u + INV_M11*(tau_foil_u + tau_rudr_u - damping_u*u_r);

                // dynamics of sway
                casadi::SX v_dot = nu_c_dot_v 
                                        + INV_M22*(tau_rudr_v - damping_v*v_r)
                                        + INV_M23*(tau_rudr_r - damping_r*r);

                // dynamics of yaw rate
                casadi::SX r_dot = 0 
                                    + INV_M32*(tau_rudr_v - damping_v*v_r)
                                    + INV_M33*(tau_rudr_r - damping_r*r);

                // full state dynamics
                casadi::SX nu_dot = vertcat(yaw_dot, u_dot, v_dot, r_dot);
                // expressed as a function for loop evaluation
                casadi::Function x_dot("x_dot", {sym_x, sym_u}, {nu_dot});

                // ################################################
                // ###----------------LOOP SETUP----------------###
                // ################################################

                // optimization variables
                casadi::SX
                    X = casadi::SX::sym("X", nx, N+1),
                    U = casadi::SX::sym("U", N),
                    optims = casadi::SX::sym("optims", nx*(N+1) + nu*N);

                // objective function, equlity constraints
                casadi::SX 
                    obj = 0,
                    g = casadi::SX::sym("g", nx*(N+1));

                // casadi constraints vector
                casadi::SX chi_t_dot, chi_t = chi_d;

                // casadi loop helper vars
                casadi::SX sym_du, sym_dx = casadi::SX::sym("sym_dx", nx);

                // Constraint MPC to start the trajectory from current position
                for(int j = 0; j < nx; j++)
                    sym_dx(j) = X(j,0) - sym_x0(j);
                sym_dx(0) = ssa(sym_dx(0));

                // fill in the constraint vector
                for(int j = 0; j < nx; j++)
                    g(j) = sym_dx(j);

                // optimization loop
                for(int i = 0; i < N; i++){

                    // assign current state
                    for(int j = 0; j < nx; j++)
                        sym_x(j) = X(j,i);

                    // assign current input or difference in input
                    sym_u = U(i);
                    if(i > 0)
                        sym_du = U(i) - U(i-1);
                    else
                        sym_du = U(i);

                    // trajectory
                    chi_t_dot = ssa(chi_d - chi_t);
                    chi_t = ssa(chi_t + Ts * chi_t_dot);

                    // assign states for readibility
                    casadi::SX
                        psi_p = sym_x(0),
                        u_p = sym_x(1) + EPS,
                        v_p = sym_x(2),
                        r_p = sym_x(3);

                    // 0 minimizes the different in course angle'
                    casadi::SX delta_x;
                    if(cost_type == 0){
                        casadi::SX beta = atan(sym_x(2) / sym_x(1) + EPS);
                        delta_x = ssa(chi_d - sym_x(0) - beta);
                    }
                    else if(cost_type == 1){
                        casadi::SX
                            x_dot = u_p * cos(psi_p) - v * sin(psi_p),
                            y_dot = u_p * sin(psi_p) + v * cos(psi_p),
                            U = sqrt(pow(x_dot,2) + pow(y_dot, 2)),
                            vec_chi_p = casadi::SX::sym("vec_chi_p", 2);
                        vec_chi_p(0) = 1/U * x_dot;
                        vec_chi_p(1) = 1/U * y_dot;

                        casadi::SX vec_chi_d = casadi::SX::sym("vec_chi_d", 2);
                        vec_chi_d(0) = cos(chi_t);
                        vec_chi_d(0) = sin(chi_t);

                        delta_x = 1 - mtimes(vec_chi_d.T(), vec_chi_p);
                    }

                    casadi::SX cost_x  = delta_x * Q * delta_x;
                    casadi::SX cost_u  = sym_du * R * sym_du;
                    obj = obj + cost_u + cost_x;

                    // multiple shooting using Runge-Kutta4
                    casadi::SXDict args, f_eval;
                    // Stage 1
                    args["i0"] = sym_x;
                    args["i1"] = sym_u;
                    f_eval = x_dot(args);
                    casadi::SX rk1 = f_eval["o0"];

                    // Stage 2
                    args["i0"] = sym_x + 0.5*Ts*rk1;
                    args["i1"] = sym_u;
                    f_eval = x_dot(args);
                    casadi::SX rk2 = f_eval["o0"];

                    // Stage 3
                    args["i0"] = sym_x + 0.5*Ts*rk2;
                    args["i1"] = sym_u;
                    f_eval = x_dot(args);
                    casadi::SX rk3 = f_eval["o0"];

                    // Stage 4
                    args["i0"] = sym_x + Ts*rk3;
                    args["i1"] = sym_u;
                    f_eval = x_dot(args);
                    casadi::SX rk4 = f_eval["o0"];

                    // next state
                    casadi::SX sym_x_rk4 = sym_x + (Ts/6) * (rk1 + 2*rk2 + 2*rk3 + rk4);

                    // introduce dynamics to constraints
                    for(int j = 0; j < nx; j++)
                        sym_dx(j) = X(j,i+1) - sym_x_rk4(j);
                    sym_dx(0) = ssa(sym_dx(0));

                    for(int j = 0; j < nx; j++)
                        g(nx*(i+1) + j) = sym_dx(j);

                    // push into the main vector being optimized
                    for(int j = 0; j < nx; j++)
                        optims(nx*i + j) = X(j,i);
                    optims(nx*(N+1) + i) = U(i);
                }

                for(int j = 0; j < nx; j++)
                    optims(nx*N + j) = X(j,N);

                // nlp problem
                casadi::SXDict nlp = {{"x", optims}, {"f", obj}, {"g", g}, {"p", sym_p}};

                // nlp options
                casadi::Dict opts;
                opts["ipopt.max_iter"] = 300;
                opts["ipopt.print_level"] = 0;
                opts["ipopt.acceptable_tol"] = 1e-8;
                opts["ipopt.acceptable_obj_change_tol"] = 1e-6;
                opts["ipopt.warm_start_init_point"] = "yes";

                solver = casadi::nlpsol("solver", "ipopt", nlp, opts);
                solver.generate_dependencies("nlp.c");

                // Just-in-time compilation?
                bool jit = false;
                if (jit) {
                    // Create a new NLP solver instance using just-in-time compilation
                    // casadi::Dict optsi = {"compiler": "shell", "jit": True, "jit_options": {"compiler": "gcc","flags": ["-O3"]}};
                    solver = casadi::nlpsol("solver", "ipopt", "nlp.c");
                } else {
                    // Compile the c-code
                    int flag = system("gcc -fPIC -shared -O3 nlp.c -o nlp.so");
                    casadi_assert(flag==0, "Compilation failed");

                    // Create a new NLP solver instance from the compiled code
                    solver = casadi::nlpsol("solver", "ipopt", "nlp.so");
                }

                // define state bounds
                std::vector<double> ubx, lbx, ubg, lbg;
                for(int i = 0; i < nx*(N+1); i++){
                    lbx.push_back(-casadi::inf);
                    ubx.push_back(casadi::inf);
                    lbg.push_back(0);
                    ubg.push_back(0); 
                }
                for(int i = nx*(N+1); i < nx*(N+1)+nu*N; i++){
                    lbx.push_back(-DEG2RAD(40));
                    ubx.push_back(DEG2RAD(40));
                }

                args_["lbx"] = lbx;
                args_["ubx"] = ubx;
                args_["lbg"] = lbg;
                args_["ubg"] = ubg;

                args_["x0"] = generate_random_vector(nx*(N+1)+nu*N);
                args_["lam_x0"] = generate_random_vector(nx*(N+1)+nu*N);
                args_["lam_g0"] = generate_random_vector(nx*(N+1));
            }

            // Function to load defaults for config, params and system dynamics
            bool loadDefaults(){

                // set state;
                state_["psi"] = 0;  // [rad]
                state_["u"] = 0.9;  // [m/s]
                state_["v"] = 0;    // [m/s]
                state_["r"] = 0;    // [rad/s]

                // get system dynamics
                std::string file = "system.csv";
                std::cout << loadDefaultsFromFile(file, system_) << std::endl;
                
                // get MpcConfigurationParameters
                file = "config.csv";
                std::cout << loadDefaultsFromFile(file, config_) << std::endl;

                return true;

            }

            bool loadDefaultsFromFile(const std::string &file_name, std::map<std::string, double> &data_from_file){

                std::string file = fs::current_path().parent_path().string() + "/autonaut/matlab_gen/" + file_name;

                std::ifstream myFile(file);
                std::string line;

                if (myFile.fail()){
                    std::cout << "ERROR: FILE OPEN FAILED. " << myFile.is_open() << std::endl;
                    std::cout << "ERROR: LOOKING AT: " << file << std::endl;
                    return false;               
                }

                while (std::getline(myFile, line)) {

                    // create a stringstream to read the data
                    std::istringstream iss(line);

                    // create key and value variables to store the data
                    std::string key;
                    double value;

                    // skip this line if unable to read both key and value
                    if (!(iss >> key >> value))
                        continue;                   

                    // store the key and value to map 
                    data_from_file[key] = value;
                }
                myFile.close();        

                return true;    
            }

            bool reWriteParams(std::vector<double> &param){

                // set initial state
                param[0] = state_["psi"];
                param[1] = state_["u"];
                param[2] = state_["v"];
                param[3] = state_["r"];
                // set other params                
                param[nx] = reference_;
                param[nx+1] = params_["Vc"];
                param[nx+2] = params_["beta_c"];
                param[nx+3] = params_["Vw"];
                param[nx+4] = params_["beta_w"];
                param[nx+5] = params_["k_1"];
                param[nx+6] = params_["k_2"];
                param[nx+7] = config_["Q"];
                param[nx+8] = config_["R"];

                return true;
            }

            // generates random vector for warm start
            std::vector<double> generate_random_vector(int n) {
                std::vector<double> result(n);
                std::random_device rd; // obtain a random seed from the OS
                std::mt19937 gen(rd()); // seed the generator
                std::uniform_real_distribution<> distr(EPS, 1.0); // define the range
                for (int i = 0; i < n; ++i) {
                    result[i] = distr(gen); // generate the random number and assign it to the vector
                }
                return result;
            }

        public:
            bool updateMpcParams(const std::map<std::string, double> &param){
                // Controller is ready to run when parameters are set
                if(initialized == false)
                    initialized = true;

                // Update defaul parameter list
                for (auto i : param) 
                    params_[i.first] = i.second;
            }

            bool updateMpcState(const std::map<std::string, double> &state){
                // flag to check if state was updated
                state_update = true;
                // Update vehicle state
                for (auto i : state) 
                    state_[i.first] = i.second;
                
            }

            bool updateMpcConfig(const std::map<std::string, double> &config){
                config_update = true;
                // Update Mpc Configuration parameters
                for (auto i : config) 
                    config_[i.first] = i.second;
            }

            bool updateMpcReference(const double &reference){
                // Update vehicle course reference
                reference_ = reference;
            }

            bool optimizeMpcProblem(){

                if(state_update == false){
                    std::cerr << "Update state before solving NLP\n";
                    return false;
                }

                std::map<std::string, casadi::DM> arg, res;
                // set state and input constraints
                arg["lbx"] = args_["lbx"];
                arg["ubx"] = args_["ubx"];
                arg["lbg"] = args_["lbg"];
                arg["ubg"] = args_["ubg"];
                
                // set Mpc parameters
                std::vector<double> param(np, 0);
                reWriteParams(param);
                arg["p"] = param;
                
                // set initial trajectory for warm start
                arg["x0"]  = args_["x0"];
                arg["lam_x"]  = args_["lam_x"];
                arg["lam_g"]  = args_["lam_g"];

                res = solver(arg);
                t_update = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

                // update variables for warm start
                // args_["x0"]  = res.at("x0");
                // args_["lam_x"]  = res.at("lam_x");
                // args_["lam_g"]  = res.at("lam_g");

                // get optimal input trajectory
                std::vector<double> optimized_vars(res.at("x"));
                // input_traj_.clear();
                // for(int i = 0; i < nu*N; i++){
                //     int index = floor(nx*(N+1) + i);
                //     double test = optimized_vars.at(0);
                //     // input_traj_.push_back();
                //     // input_traj_.push_back(optimized_vars(nx*(N+1) + i));
                // }

                return true;
            }

            bool getOptimalInput(double &u_star){
                double t_now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
                double t_elapsed = (t_update - t_now);
                // fail if NLP has not been run for a long time
                if(t_elapsed > 0.25*Tp){
                    std::cerr << "time since last NLP run exceeds threshold\n";
                    return false;
                }
                // otherwise, find the closest time index and send that input
                int t_ind = floor(t_elapsed/Ts);
                u_star = input_traj_[t_ind];
            }

    // Destructor
    ~CourseController() { 
        std::cout << "MyCallback is destroyed here." << std::endl; };
    };
}

int main(){

    NMPC::CourseController NMPC ;

    // std::cout<< NMPC.solveProblem() << std::endl;

    // std::cout << g(sym_x) << std::endl;

    return 0;
}