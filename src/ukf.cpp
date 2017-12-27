#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

//#define DEBUG_PRED_3  // Used for testing prediction
//#define DEBUG_RAD_1
//#define DEBUG_RAD_2
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true ;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // TODO

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.1 ;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...

  */
  n_x_   = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_; 

  /* Weights of sigma points */
  weights_ = VectorXd(2*n_aug_ +1);
  Xsig_pred_ = MatrixXd(n_x_,2*n_aug_+1);
  is_initialized_ = false;
  #ifdef DEBUG_PRED_1
    std_a_ = 0.2;
    std_yawdd_ = 0.2 ;
  #endif
  #ifdef DEBUG_RAD_1
    std_radr_ = 0.3;
    std_radphi_ = 0.0175;
    std_radrd_ = 0.1;
  #endif
}

UKF::~UKF() {}

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
  if (!is_initialized_) { 
  // TODO : Complete initialization
  // set x and P

  P_  << .51,0,0,0,0, 
         0,.51,0,0,0,
	 0,0,1,0,0,
	 0,0,0,1,0,
	 0,0,0,0,1;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR ) {
    /**
    Convert radar from polar to cartesian coordinates and initialize state.
    x = rho.cos(phi), y = rho.sin(phi)
    */
    cout  << "Initializing for Radar data"<<endl;
    float theta = meas_package.raw_measurements_(1);
    float ro    = meas_package.raw_measurements_(0); 
    x_  <<  ro*cos(theta),ro*sin(theta),0,0,0;
    previous_timestamp_  = meas_package.timestamp_;
    is_initialized_ = true;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER ) {
    /**
    Initialize state with position and 0 velocity
    */
    cout  << "Initializing for Laser data"<<endl;
    x_ << meas_package.raw_measurements_(0),meas_package.raw_measurements_(1),0,0.0,0.0 ;
    previous_timestamp_  = meas_package.timestamp_;
    is_initialized_ = true;
  }


  // done initializing, no need to predict or update
  cout  << "Initializing Completed"<<endl;
  return;
  }

  float dt = 0.0f;
  dt = float((meas_package.timestamp_ - previous_timestamp_)/1000000.0); 
  previous_timestamp_ = meas_package.timestamp_;
  #ifdef DEBUG_PRED_2
    dt = 0.1;
  #endif
  cout  << "Calling prediction with dt = "<<dt<<endl;
  Prediction(dt);
  cout  << "Prediction completed "<<endl;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    cout  << "Caling Radar Update "<<endl;
    UpdateRadar (meas_package);
    cout  << "Radar Update completed "<<endl;
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    cout  << "Caling Lidar Update "<<endl;
    UpdateLidar (meas_package);
    cout  << "Lidar Update completed "<<endl;
  }

  
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
  // Step-1 : Use the existing state x k|k to create sigma points X(cal) k|k
  //        : We will have to create augmented sigma points
  // Step-2 : Pass these points through the process model to get X(cal) K+1|k 
  //        : These are reduced dimension since they are not augmented
  // Step-3 : Predict mean and covariance of points from STep2
  //        : This gives x K+1|k and P k+1|k

  
  #ifdef DEBUG_PRED_1
    x_ <<   5.7441,
         1.3800,
         2.2049,
         0.5015,
         0.3528;
    P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;

  #endif
  VectorXd x_aug;
  x_aug = VectorXd(n_aug_); 
  
  x_aug.head(5) = x_; // Assign first 5 elmeents as x
  x_aug(5)      = 0;
  x_aug(6)      = 0; //Assign last 2 elements as 0 (noise mean =0)

  MatrixXd P_aug; //Augment process covariance matrix
  P_aug = MatrixXd (n_aug_, n_aug_);

  P_aug.fill (0.0); // Fill with zeros
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) =  std_a_ * std_a_;
  P_aug(6,6) =  std_yawdd_ * std_yawdd_;

  // Step-1 : Create sigma points
  /* Create X (cal) k |k using x k|k,  x k|k +/- scale.sqrt(P)
   * scale = sqrt(lambda + nx) , nx = x.size()
   * X (cal) k|k is stored as X_sig_aug
  */
  
  MatrixXd L = P_aug.llt().matrixL(); // Sqrt matrix (7x7)


  MatrixXd X_sig_aug = MatrixXd(n_aug_,n_aug_*2 +1); //7x15
  // Each column represents 1 sigma point and we have 15

  X_sig_aug.col(0) = x_aug;
  for (int i=0;i<n_aug_;i++) {
    X_sig_aug.col(1+i) = x_aug + sqrt(lambda_ + n_aug_)* L.col(i);
    X_sig_aug.col(n_aug_+1+i) = x_aug - sqrt(lambda_ + n_aug_)* L.col(i);
  }

  // Step-2: Pass sigma points through process function
  // Create X(cal) k+1 |k
  // Stored as Xsig_pred_ 5x15
 
  #ifdef DEBUG_PRED_2
  X_sig_aug <<    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
      1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
         0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
         0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;
  #endif
  //cout << "x Prior Prediction = " << x_ << endl;
  /* Loop through each sigma point */ 
  for (int i=0;i<2*n_aug_+1;i++) {
    double p_x      = X_sig_aug(0,i);
    double p_y      = X_sig_aug(1,i);
    double v        = X_sig_aug(2,i);
    double yaw      = X_sig_aug(3,i);
    double yawd     = X_sig_aug(4,i);
    double nu_a     = X_sig_aug(5,i);
    double nu_yawdd = X_sig_aug(6,i);
#if 0
    double cos_yaw   = cos(yaw);
    double sin_yaw   = sin(yaw);
    double dt_cos    = delta_t * cos_yaw;
    double dt_sin    = delta_t * sin_yaw;
    // Check for divide by zero for yawdd and normalize yaw
    if (fabs(yawd) <= 0.001) {
      Xsig_pred_(0,i) = p_x + v*dt_cos + 0.5*delta_t*dt_cos*nu_a;
      Xsig_pred_(1,i) = p_y + v*dt_sin + 0.5*delta_t*dt_sin*nu_a;
      Xsig_pred_(2,i) = v + 0 + delta_t*nu_a;
      Xsig_pred_(3,i) = yaw + 0 + 0.5 * delta_t*delta_t*nu_yawdd;
      Xsig_pred_(4,i) = yawd +0 + delta_t * nu_yawdd;
    }
    else {
      Xsig_pred_(0,i) = p_x + ((v/yawd)*(sin(yaw+ yawd*delta_t) - sin_yaw)) + (0.5 * delta_t*dt_cos*nu_a);  
      Xsig_pred_(1,i) = p_y + ((v/yawd)*(-1*cos(yaw+ yawd*delta_t) +cos_yaw)) + (0.5 * delta_t*dt_sin*nu_a);
      Xsig_pred_(2,i) = v + 0 + delta_t*nu_a;
      Xsig_pred_(3,i) = yaw + yawd*delta_t + 0.5 * delta_t*delta_t*nu_yawdd;
      Xsig_pred_(4,i) = yawd +0 + delta_t * nu_yawdd;
    }
#else
    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
#endif


    while (Xsig_pred_(3,i) > 3.14) {
      Xsig_pred_(3,i) = Xsig_pred_(3,i) - 2*3.14;
    }
    while (Xsig_pred_(3,i) < -3.14) {
      Xsig_pred_(3,i) = Xsig_pred_(3,i) + 2*3.14;
    }
  } //for i


  // Step-3: Calculate mean and covariance to get x k+1|k and P k+1|k
  #ifdef DEBUG_PRED_3  //Use precomputed Xsig_pred_ values 
  Xsig_pred_ << 5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  #endif

  #ifdef DEBUG_RAD_1  //Use precomputed Xsig_pred_ values 
  Xsig_pred_ << 5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  #endif


  /* Calculate weights */
  weights_(0) = lambda_ / (lambda_ +  n_aug_);
  for (int i=1; i< (2* n_aug_) +1; i++){
    weights_(i)  = 1/(2*(lambda_ + n_aug_));
  }

  /* Calculate new mean x_ k+1|k */
  x_  << 0.0,0.0,0.0,0.0,0.0;
  for (int i=0; i< (2* n_aug_) +1; i++){
    x_ = x_ + (weights_(i) * Xsig_pred_.col(i)); 
  }
  while (x_(3) > 3.14) {
    x_(3) = x_(3) - 2*3.14;
  }
  while (x_(3) < -3.14) {
    x_(3) = x_(3) + 2*3.14;
  }
  //cout << "Predicted x = " << x_ << endl;
  /* Calculate new Process Covariance (P) k+1 |k */



  P_.fill(0.0);
  VectorXd temp_x = VectorXd(n_x_);
  for (int i=0; i< (2* n_aug_) +1; i++){
    temp_x  = Xsig_pred_.col(i) - x_;
    while (temp_x(3) > 3.14) {
      temp_x(3) = temp_x(3) - 2*3.14;
    }
    while (temp_x(3) < -3.14) {
      temp_x(3) = temp_x(3) + 2*3.14;
    }
    P_ = P_ + weights_(i) * temp_x * temp_x.transpose();
  }

  #ifdef DEBUG_PRED_1  
    cout <<"Debug X_sig_aug "<< X_sig_aug <<endl;
    exit(-1);
  #endif
  #ifdef DEBUG_PRED_2  
    cout <<"Debug Xsig_pred_ "<< Xsig_pred_ <<endl;
    exit(-1);
  #endif
  #ifdef DEBUG_PRED_3  
    cout <<"x_ "<< x_ <<endl;
    cout <<"P_ "<< P_ <<endl;
    exit(-1);
  #endif
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Step-1 : Take sigma points and process them through State transformation
  // funciton (Converting from state space to measurement space). This gives
  // Z (cal) k+1 |k -> Denoted as Zsig 3x15
  // Step-2 : Calculate mean and variance of the transformed sigma points
  // This provides z k+1|k and S k+1|k - covariance matrix of predicted meas
  // z k+1|k is denoted as z_pred
  
  // Step-1
  MatrixXd Zsig = MatrixXd (2, (2*n_aug_) +1); // Lidar meas is 2 x15
  for (int i=0; i< (2* n_aug_) +1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);

    // measurement model
    Zsig(0,i) = p_x;                        //x
    Zsig(1,i) = p_y;                        //y  

  }

  int n_z = 2;

  // Step2.1 : Mean z k+1|k
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // Step2.1 : Mean S k+1|k

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_laspx_*std_laspx_, 0, 
          0, std_laspy_*std_laspy_;
  S = S + R;

  // Step3: Kalman filter eqs
  VectorXd z = VectorXd(n_z);
  z(0)       =  meas_package.raw_measurements_(0);
  z(1)       =  meas_package.raw_measurements_(1);
  MatrixXd Tc = MatrixXd(n_x_, n_z); //Cross Correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    while (x_diff(3)> 3.14) x_diff(3)-=2.*3.14;
    while (x_diff(3)<-3.14) x_diff(3)+=2.*3.14;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = Tc * S.inverse();
  //residual
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  while (x_(3)> 3.14) x_(3)-=2.*3.14;
  while (x_(3)<-3.14) x_(3)+=2.*3.14;
  P_ = P_ - K*S*K.transpose();
//  double NIS_E = z_diff.transpose * S.inverse() * z_diff; 

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
  
  // Step-1 : Take sigma points and process them through State transformation
  // funciton (Converting from state space to measurement space). This gives
  // Z (cal) k+1 |k -> Denoted as Zsig 3x15
  //
  // Step-2 : Calculate mean and variance of the transformed sigma points
  // This provides z k+1|k and S k+1|k - covariance matrix of predicted meas
  // z k+1|k is denoted as z_pred
  //
  // Step-3: Apply standard kalman filter eqs though Kalman gain is modified
  

  // Step-1
  MatrixXd Zsig = MatrixXd (3, (2*n_aug_) +1); // Radar meas is 3 x15
  for (int i=0; i< (2* n_aug_) +1; i++){
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    while (Zsig(1,i)> 3.14) Zsig(1,i)-= 2.*3.14;
    while (Zsig(1,i)< -3.14) Zsig(1,i)+= 2.*3.14;
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot 

  }
  int n_z = 3;

  //cout << "Radar Step -1 Completed " << endl;

  // Step2.1 : Mean z k+1|k
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  while (z_pred(1)> 3.14)  z_pred(1)-= 2.*3.14;
  while (z_pred(1)< -3.14) z_pred(1)+= 2.*3.14;

  // Step2.1 : Mean S k+1|k

  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> 3.14) z_diff(1)-=2.*3.14;
    while (z_diff(1)<-3.14) z_diff(1)+=2.*3.14;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  //cout << "Radar Step -2 Completed " << endl;

  #ifdef DEBUG_RAD_1
    cout <<"z_pred = " <<z_pred <<endl;
    cout <<" S = " <<S <<endl;
    exit (-1);
  #endif


  // Step3: Kalman filter eqs
  VectorXd z = VectorXd(n_z);
  z(0)       =  meas_package.raw_measurements_(0);
  z(1)       =  meas_package.raw_measurements_(1);
  z(2)       =  meas_package.raw_measurements_(2);

  #ifdef DEBUG_RAD_2
    Xsig_pred_ <<5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
     x_ <<      5.93637,
               1.49035,
               2.20528,
              0.536853,
              0.353577;
     P_ << 0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
          -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
          0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
         -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
         -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;
     Zsig <<      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
                  0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
                  2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

     z_pred <<       6.12155,
                     0.245993,
                     2.10313;

     S <<  0.0946171, -0.000139448,   0.00407016,
          -0.000139448,  0.000617548, -0.000770652,
           0.00407016, -0.000770652,    0.0180917;
      z <<       5.9214,
                 0.2187,
                 2.0062; 
  #endif
  //cout << "z measured = " << z<< endl;
  MatrixXd Tc = MatrixXd(n_x_, n_z); //Cross Correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    /*
    if (i==0) {
      cout << "z_pred =" << z_pred<<endl;
      cout << "Zsig.col(0) ="  << Zsig.col(0)<<endl;
      cout << "Z_diff ="  << z_diff<<endl;
      cout << "xsigpred =" <<Xsig_pred_.col(0) <<endl;
    }
    */
    
    //angle normalization
    while (z_diff(1)> 3.14) {
      z_diff(1)-=2.*3.14;
    }
    while (z_diff(1)<-3.14) {
      z_diff(1)+=2.*3.14;
    }

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> 3.14) x_diff(3)-=2.*3.14;
    while (x_diff(3)<-3.14) x_diff(3)+=2.*3.14;
    
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  MatrixXd K = Tc * S.inverse();
  //residual
  VectorXd z_diff = z - z_pred;
  //angle normalization
  while (z_diff(1)> 3.14) z_diff(1)-=2.*3.14;
  while (z_diff(1)<-3.14) z_diff(1)+=2.*3.14;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  while (x_(3)> 3.14) x_(3)-=2.*3.14;
  while (x_(3)<-3.14) x_(3)+=2.*3.14;
  P_ = P_ - K*S*K.transpose();

  #ifdef DEBUG_RAD_2
    cout << "x_ = " << x_ <<endl;
    cout << "P_ = " << P_ <<endl;
    exit (-1);
  #endif
//  double NIS_E = z_diff.transpose * S.inverse() * z_diff; 
}
