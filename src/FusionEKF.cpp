#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;
  previous_timestamp_ = 0;

  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225,      0,
              0     , 0.0225;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0     ,    0,
              0,    0.0009,    0,
              0,    0     , 0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if (!is_initialized_) {
    VectorXd initial_x(4);
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double px, py;
      tie(px, py) = tools.PolarToCartesian(
        measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1]
      );
      initial_x << px, py, 0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      initial_x << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    kalman_filter_.Initialize(initial_x, H_laser_, R_laser_, R_radar_, 9, 9);
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;

    return;
  }

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  kalman_filter_.Predict(dt);
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    kalman_filter_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    kalman_filter_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << kalman_filter_.x_ << endl;
  cout << "P_ = " << kalman_filter_.P_ << endl;
}
