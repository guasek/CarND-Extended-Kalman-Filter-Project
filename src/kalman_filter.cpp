#include <iostream>
#include "kalman_filter.h"
#include "float.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict(double dt, int noise_ax, int noise_ay) {
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  Q_ << pow(dt, 4) / 4 * noise_ax, 0                        , pow(dt, 3) / 2 * noise_ax, 0                        ,
        0                        , pow(dt, 4) / 4 * noise_ay, 0                        , pow(dt, 3) / 2 * noise_ay,
        pow(dt, 3) / 2 * noise_ax, 0                        , pow(dt, 2) * noise_ax    , 0                        ,
        0                        , pow(dt, 3) / 2 * noise_ay, 0                        , pow(dt, 2) * noise_ay    ;

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z, const MatrixXd &R) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const MatrixXd &R, const MatrixXd &Hj) {

  VectorXd z_pred = VectorXd(3);

  double rho = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
  if (rho == 0) {
    rho = DBL_MIN;
  }
  z_pred(0) = rho;
  z_pred(1) = atan2(x_(1), x_(0));
  z_pred(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / (rho);


  VectorXd y = z - z_pred;
  cout << "H: " << Hj << endl;
  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * Hj) * P_;

}
