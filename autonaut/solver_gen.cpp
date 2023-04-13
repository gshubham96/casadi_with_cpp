#include <iostream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <filesystem>
#include <casadi/casadi.hpp>

#include <random>
#include <algorithm>
#include <cmath>

#define PI 3.14159
#define EPS 1e-9
#define DEG2RAD(angle) ((angle) * M_PI / 180.0)

namespace fs = std::filesystem;

// linear, 4, chi_d
class MpcProblem {

    private:
        // costs
        int nx, nu, np;
        double Tp, Ts, N;
        double Q, R;

        // params
        double chi_d;
        double Vc, beta_c, Vw, beta_w, k_1, k_2;
        // std::vector<double> x0;

        // bounds
        std::vector<double> ubx, lbx, ubu, lbu, lbg, ubg;

        // state vectors 
        casadi::SX psi, u, v, r;
        casadi::SX delta;
        casadi::SX sym_x, sym_u, p_x0;

        casadi::SX u_e, u_c, v_c, u_r, v_r, nu_r, U_r2, beta;

        // Optimization variables
        casadi::SX X, U;

        // Objective Function
        casadi::SX obj;

        // constraints (multiple-shooting)
        casadi::SX optims, g;

        // helper vars
        casadi::SX sym_du, sym_dx;

        double ssa(double diff) {
            while (diff < -PI) diff += 2 * PI;
            while (diff > PI) diff -= 2 * PI;
            return diff;
        }

        casadi::SX ssa(casadi::SX diff) {
            return fmod(diff, 2*PI) - PI;
        }

        casadi::Function solver;

    public:

    MpcProblem(void){

        // mpc params
        Tp = 45; Ts = 0.5; N = Tp / Ts;
        Q = 4.5; R = 3;

        // named symbolica vars
        psi = casadi::SX::sym("psi", 1);
        u = casadi::SX::sym("u", 1);
        v = casadi::SX::sym("v", 1);
        r = casadi::SX::sym("r", 1);
        delta = casadi::SX::sym("delta", 1);

        // optim vars for each shooting period 
        sym_x = vertcat(psi, u, v, r);
        sym_u = delta;

        // model dims
        nx = 4; nu = 1; np = 11;

        // constant parameters for test - Vc, beta_c, Vw, beta_w,
        // Vc = 0.35; beta_c = 1.57;
        // Vw = 5;    beta_w = 1.57;
        // k_1 = 0.1; k_2 = 0.1;
        chi_d = 0;
        Vc = 0.35; beta_c = 1.57; Vw = 5; beta_w = 1.57; k_1 = 0.9551; k_2 = -0.031775;

        // system params
        double D11, R11, INV_M11;
        D11 = 286.7200;
        R11 = -53.2158;
        INV_M11 = 0.0035;

        double D22, R22, INV_M22, INV_M23;
        D22 = 194.56;
        R22 = -100.72;
        INV_M22 = 0.0026042; 
        INV_M23 = -0.00017773;

        double D33, R33, INV_M32, INV_M33;
        D33 = 1098.6;
        R33 = 207.74;
        INV_M32 = -0.00017773; 
        INV_M33 = 0.000922343;

        // detived states
        u_e = u + EPS;
        u_e = u;    // #! TODO
        u_c = Vc * cos(beta_c - psi);
        v_c = Vc * sin(beta_c - psi);
        u_r = u_e - u_c;
        v_r = v - v_c;
        nu_r = vertcat(u_r, v_r, r);
        U_r2 = pow(u_r, 2) + pow(v_r, 2);
        beta = atan(v / u_e);

        std::cout << "1 u_c & v_r = " << u_c << " --  " << v_r << std::endl;

        // dynamics of yaw
        casadi::SX yaw_dot = r;

        // dynamics of surge
        casadi::SX nu_c_dot_u = v_c * r;
        casadi::SX tau_foil_u = (k_1 + k_2*cos(psi - beta_w - PI)) * D11;
        casadi::SX tau_rudr_u = R11 * U_r2 * delta * delta ;
        casadi::SX damping_u  = D11;

        casadi::SX u_dot = nu_c_dot_u + INV_M11*(tau_foil_u + tau_rudr_u - damping_u*u_r);

        // dynamics of sway
        casadi::SX nu_c_dot_v = -u_c * r;
        casadi::SX tau_rudr_v = R22 * U_r2 * delta * 0.5 ;
        casadi::SX damping_v  = D22;

        // dynamics of yaw rate
        casadi::SX tau_rudr_r = R33 * U_r2 * delta * 0.5 ;
        casadi::SX damping_r  = D33;

        casadi::SX v_dot = nu_c_dot_v 
                                + INV_M22*(tau_rudr_v - damping_v*v_r)
                                + INV_M23*(tau_rudr_r - damping_r*r);

        casadi::SX r_dot = 0 
                            + INV_M32*(tau_rudr_v - damping_v*v_r)
                            + INV_M33*(tau_rudr_r - damping_r*r);

        casadi::SX nu_dot = vertcat(yaw_dot, u_dot, v_dot, r_dot);
        casadi::Function x_dot("x_dot", {sym_x, sym_u}, {nu_dot});

        std::cout << "2 x_dot = " << x_dot << std::endl;

        // std::vector<double> x(4, 0), u(1,0);
        // x[1] = 0.49;  x[2] = 0.13; x[3] = -0.34;
        // std::vector<casadi::DM> args = {casadi::DM(x), casadi::DM(u)};
        // std::vector<casadi::DM> result = x_dot(args);

        // std::cout << std::fixed;
        // std::cout << std::setprecision(5);
        // std::cout << "3 x_dot(args) = " << result << std::endl;

        // optimization variables
        X = casadi::SX::sym("X", 4, N+1);
        U = casadi::SX::sym("U", 1, N);
        p_x0 = casadi::SX::sym("p_x0", 1, nx);

        obj = 0;
        g = casadi::SX::sym("g", nx*(N+1));
        optims = casadi::SX::sym("optims", nx*(N+1) + nu*N);

        // set initial state
        sym_dx = casadi::SX::sym("sym_dx", 4);
        for(int j = 0; j < nx; j++)
            sym_dx(j) = X(j,0) - p_x0(j);
        sym_dx(0) = ssa(sym_dx(0));

        for(int j = 0; j < nx; j++)
            g(j) = sym_dx(j);

        // optimization loop
        // for(int i = 0; i < N; i++){
        for(int i = 0; i < 2; i++){
            
            for(int j = 0; j < nx; j++)
                sym_x(j) = X(j,i);

            sym_u = U(i);
            if(i > 0)
                sym_du = U(i) - U(i-1);
            else
                sym_du = U(i);

            casadi::SX psi_p = sym_x(0);
            casadi::SX u_p = sym_x(1) + EPS;
            casadi::SX v_p = sym_x(2);
            casadi::SX r_p = sym_x(3);
            
            casadi::SX cost_x = (chi_d - psi_p - atan(v_p/u_p)) * Q * (chi_d - psi_p - atan(v_p/u_p));
            casadi::SX cost_u = sym_du * R * sym_du;

            obj = obj + cost_u + cost_x;

            // multiple shooting using Runge-Kutta4
            casadi::SXDict args, f_eval;
            // // Stage 1
            // args["i0"] = sym_x;
            // args["i1"] = sym_u;
            // f_eval = x_dot(args);
            // casadi::SX rk1 = f_eval["o0"];

            // // Stage 2
            // args["i0"] = sym_x + 0.5*Ts*rk1;
            // args["i1"] = sym_u;
            // f_eval = x_dot(args);
            // casadi::SX rk2 = f_eval["o0"];

            // // Stage 3
            // args["i0"] = sym_x + 0.5*Ts*rk2;
            // args["i1"] = sym_u;
            // f_eval = x_dot(args);
            // casadi::SX rk3 = f_eval["o0"];

            // // Stage 4
            // args["i0"] = sym_x + Ts*rk3;
            // args["i1"] = sym_u;
            // f_eval = x_dot(args);
            // casadi::SX rk4 = f_eval["o0"];

            // Euler
            args["i0"] = sym_x;
            args["i1"] = sym_u;
            f_eval = x_dot(args);
            casadi::SX eul = f_eval["o0"];

           // next state
            // casadi::SX sym_x_rk4 = sym_x + (Ts/6) * (rk1 + 2*rk2 + 2*rk3 + rk4);
            casadi::SX sym_x_rk4 = sym_x + Ts * eul;

            casadi::SX x_n = casadi::SX::sym("x_n", 4);
            for(int j = 0; j < nx; j++)
                x_n(j) = X(j,i+1);

            for(int j = 0; j < nx; j++)
                sym_dx(j) = X(j,i) - sym_x_rk4(j);
            sym_dx(0) = ssa(sym_dx(0));

            for(int j = 0; j < nx; j++){
                g(nx*(i+1) + j) = sym_dx(j);
                std::cout << "g." << nx*(i+1) + j << " = " << g(nx*(i+1) + j) << std::endl;
            }

            // push into main vector being optimized
            for(int j = 0; j < nx; j++)
                optims((nx+nu)*i + j) = X(j,i);
            optims((nx+nu)*i + nx) = U(i);
        }

        // for(int j = 0; j < nx*(N+1); j++)
        // std::cout << "5 obj = " << obj << std::endl;

        for(int j = 0; j < nx; j++)
            optims((nx+nu)*N + j) = X(j,N);

        // nlp problem
        casadi::SXDict nlp = {{"x", optims}, {"f", obj}, {"g", g}, {"p", p_x0}};

        // nlp options
        // casadi::Dict opts;
        // opts["max_iter"] = 300;
        // opts["print_level"] = 3;
        // opts["acceptable_tol"] = 1e-8;
        // opts["acceptable_obj_change_tol"] = 1e-6;
        // TODO first try withut warm start
        // opts["warm_start_init_point"] = "yes";

        // solver = nlpsol("solver", "ipopt", nlp, opts);
        // solver = nlpsol("solver", "ipopt", nlp);
        // solver.print_options();

        // define state bounds
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
    }

    bool solveProblem(){

        // TODO 
        std::vector<double> p0 = {0.091855, 0.9821, 0.19964, 0.031876};

        std::map<std::string, casadi::DM> arg, res;
        arg["lbx"] = lbx;
        arg["ubx"] = ubx;
        arg["lbg"] = lbg;
        arg["ubg"] = ubg;
        // arg["x0"] = x0;
        arg["p"] = p0;

        // res = solver(arg);
        // std::vector<double> optimized_vars = res["x"];
        // std::cout << "optimal input found that is: " << optimized_vars(0) << std::endl;
        
        // std::ofstream file;
        // std::string filename = "casadi_ipopt_results.m";
        // file.open(filename.c_str());
        // file << "% Results file from " __FILE__ << std::endl;
        // file << "% Generated " __DATE__ " at " __TIME__ << std::endl;
        // file << std::endl;
        // file << "optimized trajectory = " << optimized_vars << ";" << std::endl;

        return false;

    }

   // Destructor
   ~MpcProblem() { std::cout << "MyCallback is destroyed here." << std::endl; };


 };

int main(){

    MpcProblem NMPC ;

    std::cout<< NMPC.solveProblem() << std::endl;

    // std::cout << g(sym_x) << std::endl;

    return 0;
}