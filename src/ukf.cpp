#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.6;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  n_z_lidar_ = 2;
  n_z_radar_ = 3;
  
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  
  double var_px = pow(std_laspx_,2);
  double var_py = pow(std_laspy_,2);

  double var_rho_ = pow(std_radr_,2);
  double var_radphi_ = pow(std_radphi_,2);
  double var_rho_dot_ = pow(std_radrd_,2);
  
  R_lidar = MatrixXd(n_z_lidar_,n_z_lidar_);
  R_lidar << var_px, 0,
            0, var_py;

  R_radar = MatrixXd(n_z_radar_,n_z_radar_);
  R_radar << var_rho_, 0, 0,
            0, var_radphi_, 0,
            0, 0, var_rho_dot_;
  
  x_ = VectorXd(n_x_);
  P_ = MatrixXd(n_x_,n_x_);

  x_.fill(0.0);

  P_.fill(0.0);
  P_ = MatrixXd::Identity(5,5);
  // Weight
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = (double) lambda_/(lambda_ + n_aug_);  
  for(int i = 1; i < 2*n_aug_+1; ++i)
  {
    weights_(i) = (double) 1.0/2.0/(lambda_ + n_aug_);
  }
  
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
    double dt = (meas_package.timestamp_ - time_us_)/1e6;
    time_us_ = meas_package.timestamp_;
    Prediction(dt);

  if(meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if(meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
}

void UKF::angleNomalizer(double &phi)
{
  phi = atan2(sin(phi),cos(phi));
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Augmentation statement
  VectorXd x_aug = VectorXd(n_aug_);
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  // Put Augmented state
  x_aug.fill(0.0);
  x_aug.head(5) = x_;

  // Put Augmented covariacne
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = pow(std_a_,2);
  P_aug(6,6) = pow(std_yawdd_,2);

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  
  const double const_aug = sqrt(lambda_ + n_aug_);
  
  for(int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i+1) = x_aug + const_aug*L.col(i);
    Xsig_aug.col(i+n_aug_+1) = x_aug - const_aug*L.col(i);
  }
  // Augmented State sigma points calculation completed

  // Predict State sigma point 

  for(int i = 0; i < 2*n_aug_ + 1; ++i)
  {
    // Get augmetned state variables
    double px = Xsig_aug(0,i);
    double py = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    double px_pred, py_pred;

    // Avoid division by zero
    if(fabs(yawd) > 0.0001)
    {
      px_pred = px + v/yawd*(sin(yaw+yawd*delta_t)-sin(yaw));
      py_pred = py + v/yawd*(-cos(yaw+yawd*delta_t)+cos(yaw));
      
    }
    else
    {
      px_pred = px + v * delta_t * cos(yaw);
      py_pred = py + v * delta_t * sin(yaw);
    }
    
    double v_pred = v;
    double yaw_pred = yaw + yawd *delta_t;
    double yawd_pred = yawd;

    px_pred += 0.5*nu_a*pow(delta_t,2)*cos(yaw);
    py_pred += 0.5*nu_a*pow(delta_t,2)*sin(yaw);
    v_pred += nu_a*delta_t;

    yaw_pred += 0.5*nu_yawdd*pow(delta_t,2);
    yawd_pred += nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_pred;
    Xsig_pred_(1,i) = py_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = yaw_pred;    
    Xsig_pred_(4,i) = yawd_pred;  
  }

  // Predict the mean state
  x_.fill(0.0);
  /**
  for(int i = 0; i < 2*n_aug_ + 1; ++i)
  {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  **/
  x_ = Xsig_pred_*weights_;
  // Predict the covariance matrix
  P_.fill(0.0);
  for(int i = 0; i < 2*n_aug_ + 1; ++i)
  {
    VectorXd x_tild = Xsig_pred_.col(i) - x_;
    while (x_tild(3) > M_PI) x_tild(3) -=2.*M_PI;
    while (x_tild(3) < -M_PI) x_tild(3) +=2.*M_PI;
    P_ = P_ + weights_(i)*x_tild*x_tild.transpose();
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  MatrixXd Zsig_lidar = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

  VectorXd z_pred_lidar = VectorXd(n_z_lidar_);

  MatrixXd S_lidar = MatrixXd(n_z_lidar_, n_z_lidar_);

  MatrixXd Tc_lidar = MatrixXd(n_x_,n_z_lidar_);

  MatrixXd K_lidar = MatrixXd(n_x_,n_z_lidar_);

  // Calculate sigma points of predicted state w.r.t lidar meas state variables
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);

    Zsig_lidar(0,i) = px;
    Zsig_lidar(1,i) = py;
  }

  // Calculate the mean value of the sigma points
  z_pred_lidar.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred_lidar = z_pred_lidar + weights_(i) * Zsig_lidar.col(i);
  }

  // Calculate Cross correlation covariance and innovation covariance matrices
  Tc_lidar.fill(0.0);
  S_lidar.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_tild_lidar = Zsig_lidar.col(i) - z_pred_lidar;
    S_lidar = S_lidar + weights_(i)*z_tild_lidar*z_tild_lidar.transpose();

    VectorXd x_tild = Xsig_pred_.col(i) - x_;
    angleNomalizer(x_tild(3));

    Tc_lidar = Tc_lidar + weights_(i) * x_tild * z_tild_lidar.transpose();
  }

  S_lidar = S_lidar + R_lidar;
  K_lidar = Tc_lidar * S_lidar.inverse();

  VectorXd z_innov_lidar_ = meas_package.raw_measurements_ - z_pred_lidar;

  // Update state and covariance
  x_ = x_ + K_lidar*z_innov_lidar_;
  P_ = P_ - K_lidar * S_lidar * K_lidar.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  MatrixXd Zsig_radar = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

  VectorXd z_pred_radar = VectorXd(n_z_radar_);

  MatrixXd S_radar = MatrixXd(n_z_radar_, n_z_radar_);

  MatrixXd Tc_radar = MatrixXd(n_x_,n_z_radar_);

  MatrixXd K_radar = MatrixXd(n_x_,n_z_radar_);

  // Calculate sigma points of predicted state w.r.t radar meas state variables
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double rho = sqrt(px*px + py*py);

    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double vx = v * cos(yaw);
    double vy = v * sin(yaw);

    double rhod = 0.0;
    
    if(fabs(rho)>0.001)
    rhod = (px*vx + py*vy) / rho;

    Zsig_radar(0,i) = rho;
    Zsig_radar(1,i) = atan2(py,px);
    Zsig_radar(2,i) = rhod;
  }

  // Calculate the mean value of the sigma points
  z_pred_radar.fill(0.0);
  z_pred_radar = Zsig_radar*weights_;

  // Calculate Cross correlation covariance and innovation covariance matrices
  Tc_radar.fill(0.0);
  S_radar.fill(0.0);
  for(int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    VectorXd z_tild_radar = Zsig_radar.col(i) - z_pred_radar;
    angleNomalizer(z_tild_radar(1));
    S_radar = S_radar + weights_(i)*z_tild_radar*z_tild_radar.transpose();

    VectorXd x_tild = Xsig_pred_.col(i) - x_;
    angleNomalizer(x_tild(3));
    
    Tc_radar = Tc_radar + weights_(i) * x_tild * z_tild_radar.transpose();
  }

  S_radar = S_radar + R_radar;
  K_radar = Tc_radar * S_radar.inverse();

  VectorXd z_innov_radar_ = meas_package.raw_measurements_ - z_pred_radar;
  angleNomalizer(z_innov_radar_(1));
  // Update state and covariance
  x_ = x_ + K_radar*z_innov_radar_;
  P_ = P_ - K_radar * S_radar * K_radar.transpose();
}