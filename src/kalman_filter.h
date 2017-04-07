#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_
#include "Eigen/Dense"
#include "tools.h"


class KalmanFilter {

public:

  Eigen::VectorXd x_;
  Eigen::MatrixXd P_;
  Eigen::MatrixXd F_;
  Eigen::MatrixXd Q_;
  Eigen::MatrixXd H_;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd I_;

  int noise_ax_;
  int noise_ay_;

  /**
   * Constructor
   */
  KalmanFilter();

  /**
   * Destructor
   */
  virtual ~KalmanFilter();

  /**
   * Initializes a Kalman Filter with a given values.
   */
  void Initialize(Eigen::VectorXd &initial_x, const Eigen::MatrixXd &H, const Eigen::MatrixXd &R_laser,
                  const Eigen::MatrixXd &R_radar, int noise_ax, int noise_ay);

  /**
   * Prediction Predicts the state and the state covariance
   * using the process model
   * @param dt Time between k and k+1 in s
   */
  void Predict(double dt);

  /**
   * Updates the state by using standard Kalman Filter equations
   * @param z The measurement at k+1
   */
  void Update(const Eigen::VectorXd &z);

  /**
   * Updates the state by using Extended Kalman Filter equations
   * @param z The measurement at k+1
   */
  void UpdateEKF(const Eigen::VectorXd &z);

private:
    Tools tools;
};

#endif /* KALMAN_FILTER_H_ */
