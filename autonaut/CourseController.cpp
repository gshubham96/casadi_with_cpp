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

#define PI M_PI
#define EPS 1e-9
#define DEG2RAD(angle) ((angle) * M_PI / 180.0)
#define RAD2DEG(angle) ((angle) * 180.0 / M_PI)

namespace fs = std::filesystem;

namespace NMPC{

    class CourseController {

        private:
            // ##################################
            // ##-------MEMBER VARIABLES-------##
            // ##################################

            // initialization variable
            int initialized;
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
            // dynamics
            casadi::Function x_dot;
            // NLP Solver
            casadi::Function solver;
            // optimized input trajectory
            std::vector<double> input_traj_;
            
            // ##################################
            // ##-------MEMBER FUNCTIONS-------##
            // ##################################

            // wrap the angle between [-pi, pi]
            double ssa(double diff) {
                diff -= (2*PI) * floor((diff + PI) * (1 / 2*PI));
                return diff;
            }

            // wrap the angle between [-pi, pi] for SX Symbolics
            casadi::SX ssa(casadi::SX diff) {
                return diff;
                diff -= (2*PI) * floor((diff + PI) * (1 / 2*PI));
                return diff;
            }

            // Function to define and compile the NLP Optimization Problem
            bool defineMpcProblem(void){

                // mpc params
                nx = config_["nx"]; nu = config_["nu"]; np = config_["np"];
                Tp = config_["Tp"]; Ts = config_["Ts"];
                N = floor(Tp / Ts);
                double model_dim = config_["model_dim"], model_type = config_["model_type"], cost_type = config_["cost_type"];

                // named symbolica vars
                casadi::SX psi = casadi::SX::sym("psi", 1),
                        u = casadi::SX::sym("u", 1),
                        v = casadi::SX::sym("v", 1),
                        r = casadi::SX::sym("r", 1),
                        delta = casadi::SX::sym("delta", 1);

                // optim vars for each shooting period 
                casadi::SX sym_x = vertcat(psi, u, v, r);
                casadi::SX sym_u = delta;
                casadi::SX sym_p = casadi::SX::sym("p", np);

                // environmental parameters that are constant over a given horizon
                casadi::SX
                    chi_d   = sym_p(nx),
                    Vc      = sym_p(nx+1),
                    beta_c  = sym_p(nx+2),
                    Vw      = sym_p(nx+3),
                    beta_w  = sym_p(nx+4),
                    k_1     = sym_p(nx+5),
                    k_2     = sym_p(nx+6),
                    Q       = sym_p(nx+7),
                    R       = sym_p(nx+8);

                // system params
                // assign system parameters
                const double D11 = system_["D11"], R11 = system_["R11"], INV_M11 = system_["INV_M11"];
                const double D22 = system_["D22"], R22 = system_["R21"], INV_M22 = system_["INV_M22"], INV_M23 = system_["INV_M23"];
                const double D33 = system_["D33"], R33 = system_["R31"], INV_M32 = system_["INV_M32"], INV_M33 = system_["INV_M33"];
                const double CR12 = system_["CR12"], CR21 = system_["CR21"]; 

                // detived states
                casadi::SX 
                    u_e = u + EPS,
                    u_c = Vc * cos(beta_c - psi),
                    v_c = Vc * sin(beta_c - psi),
                    u_r = u_e - u_c,
                    v_r = v - v_c,
                    nu_r = vertcat(u_r, v_r, r),
                    U_r2 = pow(u_r, 2) + pow(v_r, 2),
                    beta = atan(v / u_e);

                // ################################################
                // ###----------------DYNAMICS------------------###
                // ################################################

                // CURRENTS
                casadi::SX 
                    nu_c_dot_u = v_c * r,
                    nu_c_dot_v = -u_c * r;

                // DAMPING
                casadi::SX 
                    damping_u  = D11,
                    damping_v  = D22,
                    damping_r  = D33;

                // #TODO CORIOLIS 
                casadi::SX coriolis_u, coriolis_v, coriolis_r;

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
                casadi::SX tau_wind_u, tau_wind_v, tau_wind_r;

                // dynamics of yaw
                casadi::SX yaw_dot = r;

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

                casadi::SX nu_dot = vertcat(yaw_dot, u_dot, v_dot, r_dot);
                x_dot = casadi::Function("x_dot", {sym_x, sym_u, sym_p}, {nu_dot});

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

                // casadi loop helper vars
                casadi::SX sym_du, sym_dx = casadi::SX::sym("sym_dx", nx);

                // set initial state
                for(int j = 0; j < nx; j++)
                    sym_dx(j) = X(j,0) - sym_p(j);
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

                    // assign states for readibility
                    casadi::SX
                        psi_p = sym_x(0),
                        u_p = sym_x(1) + EPS,
                        v_p = sym_x(2),
                        r_p = sym_x(3);
                    casadi::SX U = sqrt( pow(u_p,2) + pow(v_p,2) );

                    casadi::SX delta_x;
                    if(cost_type == 0){
                        casadi::SX beta = asin(v_p / U);
                        casadi::SX delta_x = ssa(chi_d - psi_p - beta);
                    }
                    else if(cost_type == 1){
                        casadi::SX
                            x_dot = u_p * cos(psi_p) - v_p * sin(psi_p),
                            y_dot = u_p * sin(psi_p) + v_p * cos(psi_p),
                            vec_chi_p = casadi::SX::sym("vec_chi_p", 2);
                        vec_chi_p(0) = 1/U * x_dot;
                        vec_chi_p(1) = 1/U * y_dot;

                        casadi::SX vec_chi_d = casadi::SX::sym("vec_chi_d", 2);
                        vec_chi_d(0) = cos(chi_d);
                        vec_chi_d(1) = sin(chi_d);

                        delta_x = 1 - mtimes(vec_chi_d.T(), vec_chi_p);
                    }
                    else if(cost_type == 2){
                        delta_x = chi_d - psi_p;
                    }

                    casadi::SX cost_x  = delta_x * Q * delta_x;
                    casadi::SX cost_u  = sym_du * R * sym_du;
                    obj = obj + cost_u + cost_x;

                    // multiple shooting using Runge-Kutta4
                    casadi::SXDict args, f_eval;
                    // Stage 1
                    args["i0"] = sym_x;
                    args["i1"] = sym_u;
                    args["i2"] = sym_p;
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

                    // push into main vector being optimized
                    for(int j = 0; j < nx; j++)
                        optims(nx*i + j) = sym_x(j);
                    optims(nx*(N+1) + i) = sym_u;

                }
                for(int j = 0; j < nx; j++)
                    optims(nx*N + j) = X(j,N);

                // nlp problem
                casadi::SXDict nlp = {{"x", optims}, {"f", obj}, {"g", g}, {"p", sym_p}};

                // nlp options
                casadi::Dict opts;
                opts["ipopt.max_iter"] = 1000;
                opts["ipopt.print_level"] = 3;
                opts["ipopt.acceptable_tol"] = 1e-8;
                opts["ipopt.acceptable_obj_change_tol"] = 1e-6;
                opts["ipopt.warm_start_init_point"] = "yes";

                solver = casadi::nlpsol("solver", "ipopt", nlp, opts);

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

                // setup lower and upper bounds for constraints as well as warm start
                args_["lbx"] = lbx;
                args_["ubx"] = ubx;
                args_["lbg"] = lbg;
                args_["ubg"] = ubg;

                args_["x0"] = generate_random_vector(nx*(N+1)+nu*N);
                args_["lam_x0"] = generate_random_vector(nx*(N+1)+nu*N);
                args_["lam_g0"] = generate_random_vector(nx*(N+1));

                return true;
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
                if(loadDefaultsFromFile(file, system_))
                    std::cout << "Data loading from file " << file << " succeeded!" << std::endl;             
                else{
                    std::cout << "Data loading from file " << file << " FAILED!" << std::endl;             
                    return false;
                }

                file = "config.csv";
                if(loadDefaultsFromFile(file, config_))
                    std::cout << "Data loading from file " << file << " succeeded!" << std::endl;             
                else{
                    std::cout << "Data loading from file " << file << " FAILED!" << std::endl;             
                    return false;
                }
                
                return true;

            }

            // reads data from file and stores in passed arg
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

            // cosntructs a mpc-friendly format parameter vector
            std::vector<double> reWriteParams(){
                
                std::vector<double> param_vector(np, 0);

                // set initial state
                param_vector[0] = state_["psi"];
                param_vector[1] = state_["u"];
                param_vector[2] = state_["v"];
                param_vector[3] = state_["r"];
                // set desired state
                param_vector[4] = reference_;
                // set env params                
                param_vector[5] = params_["Vc"];
                param_vector[6] = params_["beta_c"];
                param_vector[7] = params_["Vw"];
                param_vector[8] = params_["beta_w"];
                param_vector[9] = params_["k_1"];
                param_vector[10] = params_["k_2"];
                // set costs
                param_vector[11] = config_["Q"];
                param_vector[12] = config_["R"];

                return param_vector;
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
            // updates parameters such as wind, currents, etc
            // need to do it atlease once
            bool updateMpcParams(const std::map<std::string, double> &param){
                // Controller is ready to run when parameters are set
                switch(initialized){
                    case -1: std::cout << "configure problem first!\n"; return false;
                    case  0: initialized++; break;
                    default: break;
                }

                // Update defaul parameter list
                for (auto i : param) 
                    params_[i.first] = i.second;

                return true;
            }

            // updates mpc state
            bool updateMpcState(const std::map<std::string, double> &state){

                // flag to check if state was updated
                switch(initialized){
                    case -1: std::cout << "configure problem first!\n"; return false;
                    case  0: std::cout << "update parameters first!\n"; return false;
                    case  1: initialized++; break;
                    default: break;
                }

                // Update vehicle state
                for (auto i : state) 
                    state_[i.first] = i.second;

                return true;
            }

            // updates config parameters if user wants to change the NLP
            bool updateMpcConfig(const std::map<std::string, double> &config){

                // recompile NLP if configuration parameters are changed
                for(auto i : config){
                    if(i.first.compare("Q") == 0 || i.first.compare("R") == 0)
                        continue;
                    initialized = -1;
                }

                // Update Mpc Configuration parameters
                for (auto i : config) 
                    config_[i.first] = i.second;

                // relaunch the configuration function
                if(defineMpcProblem())
                    std::cout << "Problem re-configured succesfully" << std::endl;
                else{
                    std::cout << "configuration failed!\n";
                    return false;
                }
                return true;
            }

            // updates course reference angle [rad]
            bool updateMpcReference(const double &reference){
                // Update vehicle course reference
                reference_ = reference;
                return true;
            }

            bool optimizeMpcProblem(){

                // flag to check if nlp was set up and parameters were updated
                switch(initialized){
                    case -1: std::cout << "configure problem first!\n"; return false;
                    case  0: std::cout << "update parameters first!\n"; return false;
                    case  1: std::cout << "update state first!\n"; return false;
                    default: break;
                }

                std::map<std::string, casadi::DM> arg, res;
                // set state and input constraints
                arg["lbx"] = args_["lbx"];
                arg["ubx"] = args_["ubx"];
                arg["lbg"] = args_["lbg"];
                arg["ubg"] = args_["ubg"];

                // set Mpc parameters
                std::vector p = reWriteParams();
                arg["p"] = p;
                std::cout << "params: " << p << std::endl;
                
                // set initial trajectory for warm start
                arg["x0"] = args_["x0"];
                arg["lam_x0"] = args_["lam_x0"];
                arg["lam_g0"] = args_["lam_g0"];

                res = solver(arg);
                t_update = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();

                // TODO CAN BE MADE MORE EFFICIENT 
                // get optimal input trajectory
                std::vector<double> optimized_vars(res.at("x"));

                // ################# DEBUG
                // for(int i = 0; i < N; i=i+10){
                //     std::cout << "N: " << i << ", st: ";
                //     for(int j = 0; j < nx; j++)
                //         std::cout << optimized_vars[nx * i + j] << ", ";                    
                //     std::cout << "cn: " << optimized_vars[nx*(N+1)+i] << std::endl;                    
                // }
                // std::cout << "N: " << N << ", st: ";
                // for(int j = 0; j < nx; j++)
                //     std::cout << optimized_vars[nx * N + j] << ", ";                    

                std::cout.precision(3);
                double psi = optimized_vars[nx*N];
                double u = optimized_vars[nx*N+1];
                double v = optimized_vars[nx*N+2];
                double beta = atan(v/u);
                double chi = psi + beta;

                std::cout << "\n\ndesired angle      : " << p[nx] << std::endl;
                std::cout << "final heading angle: " << psi << std::endl;
                std::cout << "final course  angle: " << chi  << std::endl;

                std::ofstream file;
                std::string filename = fs::current_path().parent_path().string();
                filename = filename + "/results/solver_gen.m";
                file.open(filename.c_str());
                file << "% Results file from " __FILE__ << std::endl;
                file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
                file << std::endl;
                file << "optims = " << optimized_vars << ";" << std::endl;

                // ################# DEBUG

                input_traj_.clear();
                for(int i = 0; i < nu*N; i++)
                    input_traj_.push_back(optimized_vars[nx*(N+1) + i]);
                
                // TODO CAN BE MADE MORE EFFICIENT 
                // update variables for warm start
                std::vector<double> lam_x(res.at("lam_x"));
                std::vector<double> lam_g(res.at("lam_g"));
                args_["x0"]  = optimized_vars;
                args_["lam_x0"]  = lam_x;
                args_["lam_g0"]  = lam_g;

                initialized--;
                return true;
            }

            bool getOptimalInput(double &u_star){
                // get current time
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

                return true;
            }

            void test_dynamics(){

                // initial position
                std::vector<double> pn(4,0), p0 = {0.091855, 0.9821, 0.19964, 0.031876};

                // get rk configs
                double Ts = config_["Ts"];

                // set default parameters
                // std::map<std::string, double> params_d;
                params_["Vc"] = 0.35; params_["beta_c"] = 1.57;
                params_["Vw"] = 5; params_["beta_w"] = 1.57;
                params_["k_1"] = 0.9551; params_["k_2"] = -0.031775;
                params_["Q"] = 4.5; params_["R"] = 3;

                // this->updateMpcParams(params_d);

                // multiple shooting using Runge-Kutta4
                casadi::DMDict args, f_eval;

                // Stage 1
                args["i0"] = p0;
                args["i1"] = 0;
                args["i2"] = reWriteParams();
                f_eval = x_dot(args);
                std::vector<double> rk1(f_eval["o0"]);

                std::cout << "params used: " << reWriteParams() << std::endl;
                std::cout << "rk1: " << rk1 << std::endl;

                // Stage 2
                for(int i=0; i<nx; i++)
                    pn[i] = p0[i] + 0.5*Ts*rk1[i];
                args["i0"] = pn;
                f_eval = x_dot(args);
                std::vector<double> rk2(f_eval["o0"]);

                // Stage 3
                for(int i=0; i<nx; i++)
                    pn[i] = p0[i] + 0.5*Ts*rk2[i];
                args["i0"] = pn;
                f_eval = x_dot(args);
                std::vector<double> rk3(f_eval["o0"]);

                // Stage 4
                for(int i=0; i<nx; i++)
                    pn[i] = p0[i] + Ts*rk3[i];
                args["i0"] = pn;
                f_eval = x_dot(args);
                std::vector<double> rk4(f_eval["o0"]);

                // next state
                std::vector<double> x_rk4(4,0);
                for(int i=0; i<nx; i++){
                    rk2[i] *= 2;
                    rk3[i] *= 2;

                    x_rk4[i] = p0[i] + (Ts/6) * (rk1[i] + rk2[i] + rk3[i] + rk4[i]);
                }
                std::cout << "state: " << x_rk4 << std::endl;

            }

    // Constructor
    CourseController(){
        initialized = -1;
        // loading defaults from a file
        if(loadDefaults())
            std::cout << "default options loaded!\n";
        else{
            std::cout << "could not load default options, exiting since user did not specify initilization variables\n";
            return;
        }
        // configuring the problem with default vars
        if(defineMpcProblem()){
            std::cout << "Problem configured succesfully" << std::endl;
            initialized++;
        }
        else{
            std::cout << "configuration failed!\n";
            return;
        }
    }

    // allow user to skip configuration
    CourseController(bool flag){
        initialized = -1;
        // loading defaults from a file
        if(loadDefaults())
            std::cout << "default options loaded!\n";
        else{
            std::cout << "could not load default options, make sure to update ALL MPC configs parameters before configuration";
        }
        if(flag)
            std::cout << "skipping configuration for now!" << std::endl;
        else{
            // relaunch the configuration function
            if(defineMpcProblem()){
                std::cout << "Problem configured succesfully" << std::endl;
                initialized++;
            }
            else
                std::cout << "configuration failed!\n";
        }

    }

    // Destructor
    ~CourseController() { 
        // std::cout << "My class is destroyed here. :(" << std::endl; 
        };

    };
}

int main(){

    // instantiate a controller with default values
    NMPC::CourseController nmpc;

    // set default parameters
    std::map<std::string, double> params_d;
    params_d["Vc"] = 0.35; params_d["beta_c"] = 1.57;
    params_d["Vw"] = 0; params_d["beta_w"] = 1.57;
    params_d["k_1"] = 0.9551; params_d["k_2"] = -0.031775;
    params_d["Q"] = 4.5; params_d["R"] = 1.5;

    // nmpc.updateMpcParams(params_d);

    // update MPC state
    std::map<std::string, double> state_d;
    state_d["psi"] = 0.091855;
    state_d["u"] = 0.9821;
    state_d["v"] = 0.19964;
    state_d["r"] = 0.031876;

    nmpc.updateMpcState(state_d);

    // update MPC reference
    double chi_ref = 0.537;
    nmpc.updateMpcReference(chi_ref);

    // solve the optimization problem
    if(nmpc.optimizeMpcProblem()){
        std::cout << "optimization succesful" << std::endl;
        double u_opt;
        if (nmpc.getOptimalInput(u_opt))
            std::cout << "optimal output is: " << u_opt << std::endl;
    }
    else
        std::cout << "optimization failed :(" << std::endl;

    // nmpc.updateMpcState(state_d);
    // if(nmpc.optimizeMpcProblem())
    //     std::cout << "optimization succesful" << std::endl;
    // else
    //     std::cout << "optimization failed :(" << std::endl;


    return 0;
}