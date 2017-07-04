#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const float PI = 3.14;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initially stto false, set to tru in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // store preciously recorded time in u
  time_us_ = 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.0;

  // Process noise st1ndard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;

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

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;
  
  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);

  // Create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Z vector
  z_laser_ = VectorXd(2);

  // Initialize H_laser
  H_laser_ = MatrixXd(2, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  // Initialize measurement noise for covariance matrix
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
             0, std_laspy_ * std_laspy_;
}

UKF::~UKF() {}

// double UKF::NormalizeAngle(double x) {
//   // if (x > M_PI)
//   //   x = fmod(x - M_PI, 2*M_PI) - M_PI;
//   // if (x < -M_PI)
//   //   x = fmod(x + M_PI, 2*M_PI) + M_PI;
//   while (x> M_PI) x-=2.*M_PI;
//   while (x<-M_PI) x+=2.*M_PI;
//   return x;
// }


double UKF::NormalizeAngle(double& pValue) {
  // if (x > M_PI)
  //   x = fmod(x - M_PI, 2*M_PI) - M_PI;
  // if (x < -M_PI)
  //   x = fmod(x + M_PI, 2*M_PI) + M_PI;
  // while (x> M_PI) x-=2.*M_PI;
  // while (x<-M_PI) x+=2.*M_PI;
  // return x;
  return pValue;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  float rho;
  float phi;
  float px;
  float py;

  if (!is_initialized_) {
    cout << "UKF: " << endl;

    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 0.5, 0,
          0, 0, 0, 0, 0.5;
    
    // weights_(0) = lambda_/(lambda_+(n_aug_));
    // for (int i = 1; i < 2*n_aug_+1; i++) {
    //   weights_(i) = 0.5*(lambda_+n_aug_);
    // }

    double weight_0 = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0;
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
      double weight = 0.5/(n_aug_+lambda_);
      weights_(i) = weight;
    }


    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Read the rho and phi input from RADER censor
      rho = meas_package.raw_measurements_[0];
      phi = meas_package.raw_measurements_[1];
      // Compute px and py from rho and phi
      px = rho * cos(phi);
      py = rho * sin(phi);

      if (fabs(px) < 0.001 && fabs(py)<0.001) {
        px = 0.001;
        py = 0.001;
      }
      if (px < 0.001) {
        px = 0.001;
      }

      if (py < 0.001) {
        py = 0.001;
      }
      // float rho_dot = meas_package.raw_measurements_[2];
      // if (rho_dot < 0.001) {
      //   rho_dot = 0.001;
      // }
      // float vx = rho_dot * cos(phi);
      // float vy = rho_dot * sin(phi);
      // float v = sqrt(vx*vx + vy*vy);

      float v = 0;
      float psi = 0;
      float psid = 0;
      x_ << px, py, v, psi, psid;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[0];
      if (px < 0.001) {
        px = 0.001;
      }
      if (py < 0.001) {
        py = 0.001;
      }
      if (fabs(px) < 0.001 && fabs(py)<0.001) {
        px = 0.001;
        py = 0.001;
      }
      x_ << px, py, 0, 0, 0;
    }
    time_us_ = meas_package.timestamp_;
    
    is_initialized_ = true;
    
    return;
  }

  double secondsDivisor = 1000000.0; // 1,000,000
  double delta_t = (meas_package.timestamp_ - time_us_) / secondsDivisor;

  time_us_ = meas_package.timestamp_; 

  while (delta_t > 0.1) {
    const double dt = 0.05;
    Prediction(dt);
    delta_t -= dt;
  } 
  Prediction(delta_t);

  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  } else if(use_laser_) {
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // Create augmented mean vector
  VectorXd x_aug_ = VectorXd(7);
  // Create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(7, 7);
  // Create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Create augmented mean state
  x_aug_.fill(0.0);
  x_aug_.head(n_x_) = x_;
  x_aug_(n_x_) = 0; // nu_a
  x_aug_(n_x_+1) = 0; // nu_yawdd

  // Create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_, n_x_) = P_;
  P_aug_(n_x_, n_x_) = std_a_ * std_a_;
  P_aug_(n_x_+1, n_x_+1) = std_yawdd_ * std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug_.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug_.fill(0.0);
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L.col(i);
  }

  // Write predicted sigma points into right column
  for (int i = 0; i < 2 * n_aug_+1; i++) {
    double px = Xsig_aug_(0,i);
    double py = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);

    double yaw = Xsig_aug_(3,i);
    yaw = NormalizeAngle(yaw);
    
    double yawd = Xsig_aug_(4,i);
    
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    double px_p, py_p;

    if (fabs(yawd) > 0.001) {
      px_p = px + (v/yawd)*(sin(yaw+yawd*delta_t) - sin(yaw));
      py_p = py + (v/yawd)*(cos(yaw) - cos(yaw+yawd*delta_t));
    } else {
      px_p = px + v*cos(yaw)*delta_t;
      py_p = py + v*sin(yaw)*delta_t;
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    yaw_p = NormalizeAngle(yaw_p);
    double yawd_p = yawd;

    px_p = px_p + 0.5*delta_t*delta_t*cos(yaw)*nu_a;
    py_p = py_p + 0.5*delta_t*delta_t*sin(yaw)*nu_a;
    v_p = v_p + delta_t*nu_a;

    yaw_p = yaw_p + 0.5*delta_t*delta_t*nu_yawdd;
    yaw_p = NormalizeAngle(yaw_p);
    yawd_p = yawd_p + delta_t*nu_yawdd;
    
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  x_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++ ){
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  // VectorXd x = VectorXd(n_x_);
  // x.fill(0.0);
  // for (int i = 0; i < 2*n_aug_+1; i++ ){
  //   x = x + (weights_(i) * Xsig_pred_.col(i));
  // }

  P_.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++ ) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    // while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    // while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
  // MatrixXd P = MatrixXd(n_x_, n_x_);
  // P.fill(0.0);
  // for (int i = 0; i < 2*n_aug_+1; i++ ) {
  //   VectorXd x_diff = Xsig_pred_.col(i) - x;
  //   //angle normalization
  //   // while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
  //   // while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
  //   x_diff(3) = NormalizeAngle(x_diff(3));
  //   P = P + weights_(i) * x_diff * x_diff.transpose() ;
  // }
  // x_ = x;
  // P_ = P;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // float px = meas_package.raw_measurements_[0];
  // float py = meas_package.raw_measurements_[1];
  // z << px, py;

  z_laser_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);

  VectorXd z_pred = H_laser_ * x_;

  VectorXd y = z_laser_ - z_pred;

  MatrixXd Ht = H_laser_.transpose();
  MatrixXd PHt = P_ * Ht;

  MatrixXd S = H_laser_ * PHt + R_laser_;
  MatrixXd Si = S.inverse();

  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
  MatrixXd const Q = P_.topLeftCorner(2,2);
  MatrixXd const Qi = Q.inverse();
  NIS_laser_ = y.transpose() * Qi * y;

  // Create matrix for sigma points in masurement space
  // MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  // // double px, py, v, yaw, yawd;
  // for (int i = 0; i < 2 * n_aug_ + 1; i++) {
  //   double px = Xsig_pred_(0,i);
  //   double py = Xsig_pred_(1,i);

  //   if(fabs(px) < 0.001 && fabs(py)<0.001) {
  //     Zsig(0,i) = 0.0;
  //     Zsig(1,i) = 0.0;
  //   } else {
  //     Zsig(0,i) = px;
  //     Zsig(1,i) = py;
  //   }
  // }

  // VectorXd z_pred = VectorXd(n_z);
  // z_pred.fill(0.0);
  // for (int i = 0; i < 2*n_aug_+1; i++) {
  //   z_pred += weights_(i) * Zsig.col(i);
  // }

  // // VectorXd z_pred = H_ * x_;
  // VectorXd y = z - z_pred;
  // // while (y(1) < -PI)
  // //   y(1) += 2 * PI;
  // // while (y(1) > PI)
  // //   y(1) -= 2 * PI;
  // // y(1) = NormalizeAngle(y(1));
  // // MatrixXd Ht = H_.transpose();
  // MatrixXd S = MatrixXd(n_z, n_z);
  // S.fill(0.0);
  // for (int i = 0; i < 2*n_aug_+1; i++) {
  //   VectorXd z_diff = Zsig.col(i) - z_pred;
  //   z_diff(1) = NormalizeAngle(z_diff(1));
  //   S += weights_(i) * z_diff * z_diff.transpose();
  // }

  // S += R_;
  // // MatrixXd S = H_ * P_ * Ht + R_;
  // // MatrixXd Si = S.inverse();
  // MatrixXd Tc = MatrixXd(n_x_, n_z);
  // Tc.fill(0.0);
  // for (int i = 0; i < 2*n_aug_+1; i++) {
  //   VectorXd z_diff = Zsig.col(i) - z_pred;
  //   z_diff(1) = NormalizeAngle(z_diff(1));
  //   VectorXd x_diff = Xsig_pred_.col(i) - x_;
  //   x_diff(3) = NormalizeAngle(x_diff(3));
  //   Tc += weights_(i) * x_diff * z_diff.transpose();
  // }

  // MatrixXd Si = S.inverse();
  // // MatrixXd K = Tc*Si;
  // // MatrixXd PHt = P_ * Ht;
  // // MatrixXd K = PHt * Si;

  // MatrixXd K = Tc*Si;

  // x_ = x_ + (K * y);
  // long x_size = x_.size();
  // MatrixXd I = MatrixXd::Identity(x_size, x_size);
  // P_ = (I - K * H_) * P_;

  // // VectorXd yt = y.transpose();
  // NIS_laser_ = y.transpose()*Si*y;
  // // MatrixXd const Q = P_.topLeftCorner(2,2);
  // // MatrixXd const Q_inverse  = Q.inverse();
  // NIS_laser_ = y_.transpose() * Q_inverse * y_;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3;

  // VectorXd z = VectorXd(n_z);
  // float ro = meas_package.raw_measurements_[0];
  // float phi = meas_package.raw_measurements_[1];
  // phi = NormalizeAngle(phi);
  // float ro_dot = meas_package.raw_measurements_[2];
  // z << ro, phi, ro_dot;
  // // Create matrix for sigma points in masurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  
  // double px, py, v, yaw, yawd;
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    yaw = NormalizeAngle(yaw);
    double yawd = Xsig_pred_(4,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    if(fabs(px) < 0.001 && fabs(py)<0.001) {
      Zsig(0,i) = 0.0;
      Zsig(1,i) = 0.0;
      Zsig(2,i) = 0.0;
    } else {
      // ro = 
      // phi = 
      // phi = NormalizeAngle(phi);
      // ro_dot = (px*v1 + py*v2) / ro;
      
      Zsig(0,i) = sqrt(px*px + py*py);
      double temp_phi = atan2(py,px);
      Zsig(1,i) = NormalizeAngle(temp_phi);
      // Zsig(1,i) = NormalizeAngle(atan2(py,px));
      Zsig(2,i) = (px*v1 + py*v2 ) / sqrt(px*px + py*py); 
    }
  }

  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  z_pred(1) = NormalizeAngle(z_pred(1));
  
  // S.Zero(n_z,n_z);
  // Measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    S += weights_(i) * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0.0);
  // R(0,0) = std_radr_ * std_radr_;
  // R(1,1) = std_radphi_* std_radphi_;
  // R(2,2) = std_radrd_ * std_radrd_;
  R << std_radr_*std_radr_, 0, 0,
       0, std_radphi_*std_radphi_, 0,
       0, 0,std_radrd_*std_radrd_;
  
  S += R;

  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2*n_aug_+1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = NormalizeAngle(z_diff(1));
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  MatrixXd Si = S.inverse();
  MatrixXd K = Tc*Si;

  double const rho_in = meas_package.raw_measurements_(0);
  double const phi_in = NormalizeAngle(meas_package.raw_measurements_(1));
  // phi_in = NormalizeAngle(phi_in);
  double const rhod_in = meas_package.raw_measurements_(2);


  VectorXd z(n_z);
  z << rho_in, phi_in, rhod_in;

  VectorXd z_diff = z - z_pred;

  z_diff(1) = NormalizeAngle(z_diff(1));

  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;

  // VectorXd y = z - z_pred;
  // y(1) = NormalizeAngle(y(1));

  // while (y(1) < -PI)
  //   y(1) += 2 * PI;
  // while (y(1) > PI)
  //   y(1) -= 2 * PI;

  // x_ = x_ + K*y;
  // P_ = P_ - K*S*K.transpose();

  // VectorXd yt = y.transpose();
  // MatrixXd test = yt*Si;
  // MatrixXd test2 = test*y;
  // NIS_radar_ = y.transpose()*Si*y;
}
