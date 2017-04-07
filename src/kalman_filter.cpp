#include <iostream>
#include "kalman_filter.h"
#include "float.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Initialize(Eigen::VectorXd &initial_x, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R_laser,
                              const Eigen::MatrixXd &R_radar, int noise_ax, int noise_ay) {
  x_ = initial_x;

  P_ = MatrixXd(4, 4);
  P_ << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1000, 0,
        0, 0, 0, 1000;

  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  long x_size = initial_x.size();
  I_ = MatrixXd::Identity(x_size, x_size);

  H_ = H;
  Q_ = MatrixXd(4, 4);

  R_laser_ = R_laser;
  R_radar_ = R_radar;
  noise_ax_ = noise_ax;
  noise_ay_ = noise_ay;
}

void KalmanFilter::Predict(double dt) {
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  double dt4 = pow(dt, 4) / 4; 
  double dt3 = pow(dt, 3) / 2;
  double dt2 = pow(dt, 2);
  double dt3_noise_ay = dt3 * noise_ay_;
  double dt3_noise_ax = dt3 * noise_ax_;

  Q_ << dt4 * noise_ax_, 0              , dt3_noise_ax, 0,
        0              , dt4 * noise_ay_, 0              , dt3_noise_ay,
        dt3_noise_ax   , 0              , dt2 * noise_ax_, 0,
        0              , dt3_noise_ay   , 0              , dt2 * noise_ay_;

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd z_pred = VectorXd(3);
  MatrixXd Hj = tools.CalculateJacobian(x_);

  double rho = sqrt(pow(x_(0), 2) + pow(x_(1), 2));
  if (rho == 0) {
    rho = DBL_MIN;
  }
  z_pred(0) = rho;
  z_pred(1) = atan2(x_(1), x_(0));
  z_pred(2) = (x_(0) * x_(2) + x_(1) * x_(3)) / (rho);


  VectorXd y = z - z_pred;
  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  P_ = (I_ - K * Hj) * P_;
}