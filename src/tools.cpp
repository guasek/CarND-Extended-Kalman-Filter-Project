#include <iostream>
#include "tools.h"
#include "float.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    if (estimations.size() == 0 || estimations.size() != ground_truth.size()){
        return rmse;
    }

    for(int i=0; i < estimations.size(); ++i){
        VectorXd residual = estimations[i] - ground_truth[i];
        rmse = rmse + (residual.array() * residual.array()).matrix();
    }

    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj(3,4);
  Hj << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;

  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  if (px == 0) {
    px = DBL_MIN;
  }

  if (py == 0) {
    py = DBL_MIN;
  }

  double px2_py2_sum = pow(px, 2) + pow(py, 2);

  if (px2_py2_sum == 0) {
    px2_py2_sum = DBL_MIN;
  }

  double powered = pow(px2_py2_sum, 1.5);
  if (powered == 0) {
    powered = DBL_MIN;
  }

  cout << "sqrt: " << pow(px2_py2_sum, 1.5) << endl;
  Hj(0, 0) = px / (sqrt(px2_py2_sum));
  Hj(0, 1) = py / (sqrt(px2_py2_sum));

  Hj(1, 0) = -py / px2_py2_sum;
  Hj(1, 1) = px / px2_py2_sum;

  Hj(2, 0) = (py * (vx*py - vy*px)) / (powered);
  Hj(2, 1) = (px * (vy*px - vx*py+)) / (powered);
  Hj(2, 2) = px / (sqrt(px2_py2_sum));
  Hj(2, 3) = py / (sqrt(px2_py2_sum));

  return Hj;
}
