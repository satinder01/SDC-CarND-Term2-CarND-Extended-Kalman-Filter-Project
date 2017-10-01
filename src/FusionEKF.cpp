#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#define EPS 0.0001 // A small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  float previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
   H_laser_ << 1,0,0,0,
               0,1,0,0;

   Hj_ << 0,0,0,0,
          0,0,0,0,
          0,0,0,0;


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  float rho;
  float phi;
  float rhodot;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF initialise: " << endl;

    //state covariance matrix P
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<  1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;

    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0.0,0.0,0.0,0.0;
 
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */

      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);

      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      ekf_.x_ << px,py,vx,vy;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
        float px = measurement_pack.raw_measurements_(0);
        float py = measurement_pack.raw_measurements_(1);
        ekf_.x_ << px,py,0.0,0.0;
    }
    if (fabs(ekf_.x_(0)) < EPS and fabs(ekf_.x_(1)) < EPS){
      ekf_.x_(0) = EPS;
      ekf_.x_(1) = EPS;
    }
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1, 0,
               0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt =  (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  float c1 = dt_4 / 4;
  float c2 = dt_3 / 2;
  float noise_ax = 9;
  float noise_ay = 9;


  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ << 1, 0, dt, 0,
	     0, 1, 0, dt,
	     0, 0, 1, 0,
	     0, 0, 0, 1;

  ekf_.Q_ << (c1 * noise_ax), 0, (c2 * noise_ax), 0,
             0, (c1 * noise_ay), 0, (c2 * noise_ay),
             (c2 * noise_ax), 0, (dt_2 * noise_ax), 0,
             0, (c2 * noise_ay), 0, (dt_2 * noise_ay);

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    VectorXd zR = VectorXd(3);
    zR << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1],
				measurement_pack.raw_measurements_[2];
    
    //ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, Hj_, R_radar_, ekf_.Q_);
    //ekf_.H_= Hj_;
    ekf_.R_= R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
   //if (!Hj_.isZero(0))
      ekf_.UpdateEKF(zR);

  } else {
    // Laser updates

    VectorXd zL(2);
    zL << measurement_pack.raw_measurements_(0),measurement_pack.raw_measurements_(1);
    //ekf_.Init(ekf_.x_, ekf_.P_, ekf_.F_, H_laser_, R_laser_, ekf_.Q_);
    ekf_.H_= H_laser_;
    ekf_.R_= R_laser_;
    ekf_.Update(zL);

  }

  // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
