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
            std::vector<double> optimized_vars_, input_traj_;
            // file handling
            int filecount;
            std::ofstream file;
            std::string filename;
            
            // ##################################
            // ##-------MEMBER FUNCTIONS-------##
            // ##################################

            // wrap the angle between [-pi, pi]
            double ssa(double diff);
            casadi::SX ssa(casadi::SX diff);

            // Function to define and compile the NLP Optimization Problem
            bool defineMpcProblem(void);

            // Function to load defaults for config, params and system dynamics
            bool loadDefaults();

            // reads data from file and stores in passed arg
            bool loadDefaultsFromFile(const std::string &file_name, std::map<std::string, double> &data_from_file);

            // cosntructs a mpc-friendly format parameter vector
            std::vector<double> reWriteParams();

            // generates random vector for warm start
            std::vector<double> generate_random_vector(int n);

        public:
            // updates parameters such as wind, currents, etc
            // need to do it atlease once
            bool updateMpcParams(const std::map<std::string, double> &param);

            // updates mpc state
            bool updateMpcState(const std::map<std::string, double> &state);

            // updates config parameters if user wants to change the NLP
            bool updateMpcConfig(const std::map<std::string, double> &config);

            // updates course reference angle [rad]
            bool updateMpcReference(const double &reference);

            // Optimize the NLP with updated state and parameters
            bool optimizeMpcProblem();

            // Get the result from the optimized input
            double getOptimalInput();

            // debug function to test system dynamics
            void test_dynamics();

            // debug function to save trajectory to file
            void saveTrajectoryToFile();

            // debug function to print trajectory and other info on screen
            void print_details();

    // Constructor
    CourseController();

    // allow user to skip configuration
    CourseController(bool flag);

    // Destructor
    ~CourseController();

    };
}