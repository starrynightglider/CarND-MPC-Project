#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"

using CppAD::AD;

size_t N = 9;
double dt = 0.12;
int latencyState = 2;
double refVel = 60;
size_t xS = 0;
size_t yS = xS + N;
size_t psiS = yS + N;
size_t velS = psiS + N;
size_t cteS = velS + N;
size_t epsiS = cteS + N;
size_t deltaS = epsiS + N;
size_t accS = deltaS + N - 1;
// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;
// Setup array offset

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {

    fg[0] = 0;

    // Reference state cost
    for (int i = 0; i < N; ++i){
      fg[0] += 100*CppAD::pow(vars[cteS+i], 2);
      fg[0] += 20*CppAD::pow(vars[epsiS+i], 2);
      fg[0] += 1*CppAD::pow(vars[velS+i] - refVel, 2);
    }

    // Prefer smaller delta and acc
    for (int i = 0; i < N-1; ++i){
      fg[0] += 400*CppAD::pow(vars[deltaS+i], 2);
      fg[0] += 100*CppAD::pow(vars[accS+i], 2);
      //fg[0] += 200*CppAD::pow(vars[deltaS+i]*vars[velS+i],2);
    }
    // Reduce change rate
    for (int i = 0; i < N - 2; ++i){
      fg[0] += 600*CppAD::pow(vars[deltaS+i+1]-vars[deltaS+i], 2);
      fg[0] += 1*CppAD::pow(vars[accS+i+1] - vars[accS+i], 2);
    }

    // Setup Constraints

    // Initial constraints
    // Add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    fg[1 + xS] = vars[xS];
    fg[1 + yS] = vars[yS];
    fg[1 + psiS] = vars[psiS];
    fg[1 + velS] = vars[velS];
    fg[1 + cteS] = vars[cteS];
    fg[1 + epsiS] = vars[epsiS];


    for (int i = 1; i < N; ++i){
      AD<double> x1 = vars[xS+i];
      AD<double> y1 = vars[yS+i];
      AD<double> psi1 = vars[psiS+i];
      AD<double> vel1 = vars[velS+i];
      AD<double> cte1 = vars[cteS+i];
      AD<double> epsi1 = vars[epsiS+i];

      AD<double> x0 = vars[xS+i-1];
      AD<double> y0 = vars[yS+i-1];
      AD<double> psi0 = vars[psiS+i-1];
      AD<double> vel0 = vars[velS+i-1];
      AD<double> cte0 = vars[cteS+i-1];
      AD<double> epsi0 = vars[epsiS+i-1];

      AD<double> acc = vars[accS];
      AD<double> delta = vars[deltaS];

      // use previous actuations
      if (i > 1) {
        acc = vars[accS + i - 2];
        delta = vars[deltaS + i - 2];
      }
      AD<double> f0 = coeffs[0] + coeffs[1] * x0 + coeffs[2] * CppAD::pow(x0, 2)
                                + coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psiDes0 = CppAD::atan(coeffs[1] + 2*coeffs[2]*x0 +
                                                   3*coeffs[3]*CppAD::pow(x0, 2));
      // Updates
      // x_[t+1] = x[t] + v[t] * cos(psi[t]) * dt
      // y_[t+1] = y[t] + v[t] * sin(psi[t]) * dt
      // psi_[t+1] = psi[t] + v[t] / Lf * delta[t] * dt
      // v_[t+1] = v[t] + a[t] * dt
      // cte[t+1] = f(x[t]) - y[t] + v[t] * sin(epsi[t]) * dt
      // epsi[t+1] = psi[t] - psides[t] + v[t] * delta[t] / Lf * dt

      fg[1 + xS + i] = x1 - (x0 + vel0 * CppAD::cos(psi0) * dt);
      fg[1 + yS + i] = y1 - (y0 + vel0 * CppAD::sin(psi0) * dt);
      fg[1 + psiS + i] = psi1 - (psi0 - vel0/Lf * delta * dt);
      fg[1 + velS + i] = vel1 - (vel0 + acc * dt);
      fg[1 + cteS + i] = cte1 - ((f0 - y0) + (vel0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsiS + i] = epsi1 - ((psi0 - psiDes0) - vel0/Lf * delta * dt);
    }
  }
};

// MPC class definition implementation.

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  // Set the number of model variables (includes both states and inputs).
  double x = state[0];
  double y = state[1];
  double psi = state[2];
  double vel = state[3];
  double cte = state[4];
  double epsi = state[5];
  // For example: If the state is a 4 element vector, the actuators is a 2
  // element vector and there are 10 timesteps. The number of variables is:
  // 4 * 10 + 2 * 9
  size_t n_vars = N*6 + (N-1)*2;
  // Set the number of constraints
  size_t n_constraints = N*6;

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  vars[xS] = x;
  vars[yS] = y;
  vars[psiS] = psi;
  vars[velS] = vel;
  vars[cteS] = cte;
  vars[epsiS] = epsi;

  for (int i = 0; i < deltaS; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of delta are set to -25 and 25
  // degrees (values in radians).

  for (int i = deltaS; i < accS; ++i) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
}

  // Acceleration/decceleration upper and lower limits.

  for (int i = accS; i < n_vars; ++i) {
      vars_lowerbound[i] = -1.0;
      vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; ++i) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[xS] = x;
  constraints_lowerbound[yS] = y;
  constraints_lowerbound[psiS] = psi;
  constraints_lowerbound[velS] = vel;
  constraints_lowerbound[cteS] = cte;
  constraints_lowerbound[epsiS] = epsi;

  constraints_upperbound[xS] = x;
  constraints_upperbound[yS] = y;
  constraints_upperbound[psiS] = psi;
  constraints_upperbound[velS] = vel;
  constraints_upperbound[cteS] = cte;
  constraints_upperbound[epsiS] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);

  //
  // NOTE: You don't have to worry about these options
  //
  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.

  vector<double> result;

  result.push_back(solution.x[deltaS]);
  result.push_back(solution.x[accS]);

  for (int i = 1; i < N; i++) {
    result.push_back(solution.x[xS + i]);
    result.push_back(solution.x[yS + i]);
  }

  return result;
}
