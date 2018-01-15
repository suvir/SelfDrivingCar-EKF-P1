#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define PI 3.14159265

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

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

    step = 0;
}

FusionEKF::~FusionEKF() {}

float normalize_angle(float theta) {
    while (theta < -PI) {
        theta += 2 * PI;
    }
    while (theta > PI) {
        theta -= 2 * PI;
    }
    return theta;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    step++;

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 0, 0, 5, 0;

        ekf_.F_ = MatrixXd(4, 4);
        ekf_.F_ << 1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;

        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1000, 0,
                0, 0, 0, 1000;

        ekf_.H_ = MatrixXd(2, 4);
        ekf_.H_ << 1, 0, 0, 0,
                0, 1, 0, 0;

        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            float rho = measurement_pack.raw_measurements_(0);
            float theta = normalize_angle(measurement_pack.raw_measurements_(1));
            float rho_p = measurement_pack.raw_measurements_(2);
            //px and py
            ekf_.x_(0) = rho * cos(theta);
            ekf_.x_(1) = rho * sin(theta);
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
             Initialize state.
             */
            //px and py
            ekf_.x_(0) = measurement_pack.raw_measurements_(0);
            ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/
    float noise_ay = 9.0;
    float noise_ax = 9.0;
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; // in sec
    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    previous_timestamp_ = measurement_pack.timestamp_;
    // update F_ (state transition matrix)
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;
    // update Q_ (process noise covariance matrix)
    ekf_.Q_ = MatrixXd(4, 4);

    ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
            0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
            dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
            0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar updates
        float theta = normalize_angle(measurement_pack.raw_measurements_(1));
        VectorXd meas_data = VectorXd(3);
        meas_data << measurement_pack.raw_measurements_(0), theta, measurement_pack.raw_measurements_(2);
        ekf_.R_ = R_radar_;
        ekf_.Hj_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.Update(meas_data, 1);
    } else {
        // Laser updates
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_, 0);
    }
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}