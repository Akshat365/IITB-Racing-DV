/*
    Formula Student Driverless Project (FSD-Project).
    Copyright (c) 2018:
     - Sonja Brits <britss@ethz.ch>
    FSD-Project is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    FSD-Project is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with FSD-Project.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <ros/ros.h>
#include "slam.hpp"
#include <sstream>
#include <cmath>
#include <eigen>

namespace ns_slam {
// Constructor
Slam::Slam() {
  initializeState();
};

// Getters
std::vector<fsd_common_msgs::Cone> Slam::getMap() const { return cone_map_; }
geometry_msgs::Pose2D Slam::getState() const { return slam_state_; }

// Setters
void Slam::setMaxMapSize(int size) { max_map_size_ = size; }

void Slam::initializeState() {
  slam_state_.x = 0;
  slam_state_.y = 0;
  slam_state_.theta = 0;
  slam_state_cov_ << 100, 0, 0, 0, 100, 0, 0, 0, 100;
  slam_full_state_ << 0, 0, 0;
  slam_full_state_cov_ << 100, 0, 0, 0, 100, 0, 0, 0, 100;
}

void Slam::updateMap(const fsd_common_msgs::ConeDetections &cones) {
  // ROS_INFO("Map size is %d", (int) cone_map_.size());
  // int excess_cone_length = cone_map_.size() + cones.cone_detections.size() - max_map_size_;
  // if (cone_map_.size() + cones.cone_detections.size() > max_map_size_) {
  //   ROS_INFO("Removing %d elements to accomodate new elements", (int) excess_cone_length);
  //   cone_map_.erase(cone_map_.begin(), cone_map_.begin() + excess_cone_length);
  // }
  // for (int i = 0; i < cones.cone_detections.size(); i++) {
  //   ROS_INFO("cone add to map");
  //   cone_map_.push_back(cones.cone_detections[i]);
  // }

  int prev_N = cone_map_.size(); // vector

  Eigen::Matrix Fx = Matrix<double, 3, 2*prev_N+3>::Identity();

  Eigen::Matrix<double, 2*prev_N+3, 1> prev_state;
  prev_full_state = slam_full_state_;

  Eigen::Matrix<double, 3, 1> motion_update;
  motion_update << vt*Ts*std::cos(slam_state_.theta + velocity.car_state_dt.theta*Ts/2), vt*Ts*std::sin(slam_state_.theta + velocity.car_state_dt.theta*Ts/2), velocity.car_state_dt.theta*Ts;

  Eigen::Matrix<double, 2*prev_N+3, 1> pred_state;
  pred_full_state = prev_full_state + Fx.transpose()*motion_update;

  Eigen::Matrix<double, 3, 3> gk;
  gk << 0, 0, , -vt*Ts*std::sin(slam_state_.theta + velocity.car_state_dt.theta*Ts/2),
        0, 0, , vt*Ts*std::cos(slam_state_.theta + velocity.car_state_dt.theta*Ts/2),
        0, 0, 0;

  Eigen::Matrix<double, 2*prev_N+3, 2*prev_N+3> Gk;
  Gk = Eigen::Matrix<double, 2*prev_N+3, 2*prev_N+3>::Identity() + Fx.transpose()*gk*Fx; 

  Eigen::Matrix<double, 2, 2> R; // process noise
  R << 1, 0, 0,
       0, 1, 0,
       0, 0, 1;

  Eigen::Matrix<double, 2*prev_N+3, 2*prev_N+3> pred_full_state_cov;
  pred_full_state_cov = slam_full_state_cov_ + Fx.transpose()*R*Fx; 
  
  slam_full_state_ = pred_full_state;
  slam_full_state_cov_ = pred_full_state_cov;
  
  // prediction of new state using controls --- completed

  Eigen::Matrix<double, 2, 2> Q; // measurement noise
  Q << 1, 0,
       0, 1;

  double alpha = 0.5;

  int num_measurements = sizeof(cones.cone_detections)/sizeof(cones.cone_detections[0]); // array
  

  for (i = 0; i<num_measurements; i++){

    Matrix<double, 2, 1> new_landmark;
    r_i = cones.cone_detections[i].position.x;
    phi_i = cones.cone_detections[i].position.y;
    theta = pred_full_state(3,1);
    new_landmark << pred_full_state(1,1) + r_i*(std::cos(phi_i + theta)), pred_full_state(2,1) + r_i*(std::sin(phi_i + theta));

    Eigen::Matrix<double, (prev_N+1), 1> pi;
    pi(prev_N+1, 1) = alpha;
    double pi_min = alpha;
    int pi_min_index = 0;

    Eigen::Matrix<double, 2, 1> zti;
    zti << cones.cone_detections[i].position.x, cones.cone_detections[i].position.y; // r, phi  
    
    
    for (int k = 0; k<prev_N; k++){
      
      double del_kx = (cone_map_[k].position.x - pred_full_state(1,1));
      double del_ky = (cone_map_[k].position.y - pred_full_state(2,1)); 

      double q = std::pow(del_kx,2) + std::pow(del_ky,2);
      double r = std::sqrt(q);

      Eigen::Matrix<double, 2, 1> est_ztk;
      est_ztk << r, std::atan2(del_kx, del_ky)-theta;

      Eigen::Matrix Fxk = Eigen::Matrix::Zero(5, 2*prev_N+3);
      Fxk(0,0) = 1;
      Fxk(1,1) = 1;
      Fxk(2,2) = 1;
      Fxk(3,2*k+3) = 1;
      Fxk(4,2*k+4) = 1;
      
      Eigen::Matrix<double, 2, 5> htk;
      htk << -del_kx/r, -del_ky/r, 0, del_kx/r, del_ky/r,
              del_ky/q, -del_kx/q, -1, -del_ky/q, del_kx/q;
      
      Eigen::Matrix<double, 2, 2*prev_N+3> Htk;
      Htk = htk*Fxk;

      Eigen::Matrix<double, 2, 2*prev_N+3> psi_tk;
      psi_tk = Htk*slam_state_cov_*(Htk*transpose()) + Q;

      pi(k, 1) = ((zti-est_ztk).transpose())*(psi_tk.inverse())*(zti-est_ztk);
      if (pi(k, 1) <= pi_min){
        pi_min = pi(k, 1);
        pi_min_index = k;
      }
    }

    if (pi_min == alpha){
      // add new landmark and update slam_full_state_, slam_full_state_cov_ in the end
      fsd_common_msgs::Cone new_cone;
      new_cone.position.x = new_landmark(1,1);
      new_cone.position.y = new_landmark(2,1); 
      new_cone.color = cones.cone_detections[i].color;
      cone_map_.cone_detections.push_back(new_cone);
    }
    else{
      // update state values using the correspondence

      double del_kx = (cone_map_[pi_min_index].position.x - pred_full_state(1,1));
      double del_ky = (cone_map_[pi_min_index].position.y - pred_full_state(2,1)); 

      double q = std::pow(del_kx,2) + std::pow(del_ky,2);
      double r = std::sqrt(q);

      Eigen::Matrix<double, 2, 1> est_ztk;
      est_ztk << r, std::atan2(del_kx, del_ky)-theta;

      Eigen::Matrix Fxk = Eigen::Matrix::Zero(5, 2*prev_N+3);
      Fxk(0,0) = 1;
      Fxk(1,1) = 1;
      Fxk(2,2) = 1;
      Fxk(3,2*pi_min_index+3) = 1;
      Fxk(4,2*pi_min_index+4) = 1;
      
      Eigen::Matrix<double, 2, 5> htk;
      htk << -del_kx/r, -del_ky/r, 0, del_kx/r, del_ky/r,
              del_ky/q, -del_kx/q, -1, -del_ky/q, del_kx/q;
      
      Eigen::Matrix<double, 2, 2*prev_N+3> Htk;
      Htk = htk*Fxk;

      Eigen::Matrix<double, 2, 2*prev_N+3> psi_tk;
      psi_tk = Htk*slam_state_cov_*(Htk*transpose()) + Q;

      Eigen::Matrix<double, 2*prev_N+3, 2> K_tk;
      K_tk = slam_full_state_cov_*(Htk.transpose())*(psi_tk.inverse());

      slam_full_state_ = slam_full_state_cov_ + K_tk*(zti-est_ztk);
      slam_full_state_cov_ = ( Eigen::Matrix::Zero(2*prev_N+3, 2*prev_N+3) - K_tk*Htk)*slam_full_state_cov_;
    }
  }

  int new_cones_detected = cone_map_.size()-N_prev;
  for (int n=N_prev; n<cone_map_.size(); n++){
    slam_full_state_(2*n+3, 1) = cone_map_[n].position.x;
    slam_full_state_(2*n+4, 1) = cone_map_[n].position.y;
    slam_full_state_cov_(2*n+3, 2*n+3) = 10;
    slam_full_state_cov_(2*n+4, 2*n+4) = 10;
  }

}



void Slam::calculateState(const fsd_common_msgs::CarStateDt &velocity, const fsd_common_msgs::ConeDetections &cones) {
  // slam_state_.x += velocity.car_state_dt.x;
  // slam_state_.y += velocity.car_state_dt.y;
  // slam_state_.theta += velocity.car_state_dt.theta;  =====> velocity.car_state_dt.theta is the angular velocity
  
  Eigen::Matrix<double, 3, 1> prev_state;
  Eigen::Matrix<double, 3, 3> Gt;
  Eigen::Matrix<double, 3, 2> Vt;
  Eigen::Matrix<double, 2, 2> Mt;
  Eigen::Matrix<double, 3, 1> pred_state;
  Eigen::Matrix<double, 3, 3> pred_state_cov;
  Eigen::Matrix<double, 2, 2> Qt;

  double Ts = 1;
  double alpha1 = 0.5;
  double alpha2 = 0.5;
  double alpha3 = 0.5;
  double alpha4 = 0.5;
  double vt = std::sqrt(std::pow(velocity.car_state_dt.x, 2) + std::pow(velocity.car_state_dt.y, 2));

  prev_state << slam_state_.x, slam_state_.y, slam_state_.theta;

  Gt << 1, 0, , -vt*Ts*std::sin(slam_state_.theta + velocity.car_state_dt.theta*Ts/2),
        0, 1, , vt*Ts*std::cos(slam_state_.theta + velocity.car_state_dt.theta*Ts/2),
        0, 0, 1;
  Vt << std::cos(slam_state_.theta + velocity.car_state_dt.theta*Ts/2), -0.5*std::sin(slam_state_.theta + velocity.car_state_dt.theta*Ts/2),
        std::sin(slam_state_.theta + velocity.car_state_dt.theta*Ts/2), 0.5*std::cos(slam_state_.theta + velocity.car_state_dt.theta*Ts/2),
        0, 1;
  Mt << alpha1*vt*vt+alpha2*velocity.car_state_dt.theta*velocity.car_state_dt.theta, 0,
        0, alpha3*vt*vt+alpha4*velocity.car_state_dt.theta*velocity.car_state_dt.theta;

  pred_state << slam_state_.x + vt*Ts*std::cos(slam_state_.theta + velocity.car_state_dt.theta*Ts/2), slam_state_.y + vt*Ts*std::sin(slam_state_.theta + velocity.car_state_dt.theta*Ts/2), slam_state_.theta+velocity.car_state_dt.theta*Ts;
  pred_state_cov = Gt*slam_state_cov_*(Gt.transpose()) + Vt*Mt*(Vt.transpose());

  Qt << 1, 0,
        0, 1;
  
  int num_measurements = sizeof(cones.cone_detections)/sizeof(cones.cone_detections[0]); // array
  int map_landmarks = cone_map_.size(); // vector
 
  for (int i = 0; i<num_measurements; i++){

    double prob[map_landmarks];
    int map_index = 0;
    
    Eigen::Matrix<double, 2, 1> zti;
    zti << cones.cone_detections[i].position.x, cones.cone_detections[i].position.y; // r, phi  
    
    for (int k = 0; k<map_landmarks; k++){

      q = std::pow((cone_map_[k].position.x - pred_state(1,1)),2) + std::pow((cone_map_[k].position.y - pred_state(1,2)),2);
      
      Eigen::Matrix<double, 2, 1> est_ztk;
      Eigen::Matrix<double, 2, 3> Htk;
      Eigen::Matrix<double, 3, 3> Stk;
      
      est_ztk << std::sqrt(q), std::atan2((cone_map_[k].position.x - pred_state(1,1)), (cone_map_[k].position.y - pred_state(1,2)))-pred_state(1,3);
      Htk << -(cone_map_[k].position.x - pred_state(1,1))/std::sqrt(q), -(cone_map_[k].position.y - pred_state(1,2))/std::sqrt(q), 0,
              (cone_map_[k].position.y - pred_state(1,2))/q, -(cone_map_[k].position.x - pred_state(1,1))/q, -1;
      Stk = Htk*pred_state_cov*(Htk.transpose()) + Qt;

      prob[k] = ( std::exp(-0.5*(zti-est_ztk)*(Stk.inverse())*((zti-est_ztk).transpose())) )/std::sqrt(2*3.1416*Stk.determinant());
      
      if (prob[k] > prob[map_index]){
        map_index = k;
      }

    }
    
    q = std::pow((cone_map_[map_index].position.x - pred_state(1,1)),2) + std::pow((cone_map_[map_index].position.y - pred_state(1,2)),2);

    Eigen::Matrix<double, 2, 1> est_zti;
    Eigen::Matrix<double, 2, 3> Hti;
    Eigen::Matrix<double, 3, 3> Sti;
    Eigen::Matrix<double, 3, 2> Kti;
    
    est_zti << cones.cone_detections[i].position.x, cones.cone_detections[i].position.y; // r, phi  
    Hti << -(cone_map_[map_index].position.x - pred_state(1,1))/std::sqrt(q), -(cone_map_[map_index].position.y - pred_state(1,2))/std::sqrt(q), 0,
              (cone_map_[map_index].position.y - pred_state(1,2))/q, -(cone_map_[map_index].position.x - pred_state(1,1))/q, -1;
    Sti = Hti*pred_state_cov*(Hti.transpose()) + Qt;
    Kti = pred_state_cov*(Hti.transpose())*(Sti.inverse());
    pred_state = pred_state + Kti*(zti-est_zti);
    pred_state_cov = pred_state_cov - Kti*Hti*pred_state_cov;

  }

  slam_state_.x = pred_state(1,1);
  slam_state_.y = pred_state(2,1);
  slam_state_.theta = pred_state(3,1);

  slam_state_cov_ = pred_state_cov;

}




}