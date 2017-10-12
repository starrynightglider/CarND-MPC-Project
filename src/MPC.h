#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  MPC():prevDelta_(0), prevAcc_(0){}
  virtual ~MPC(){}

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
  void setPrevActuator(double delta, double acc){
    prevDelta_ = delta;
    prevAcc_ = acc;
  }

private:
  double prevDelta_;
  double prevAcc_;
};

#endif /* MPC_H */
