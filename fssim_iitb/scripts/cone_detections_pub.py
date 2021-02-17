#!/usr/bin/env python

import rospy
from sensor_msgs.msg import LaserScan
from fssim_iitb.msg import ConeDetections

def callback(msg):
  cone_detections = ConeDetections()
  range_arr_length = len(msg.ranges)
  range_arr = msg.ranges
  
  min_angle = msg.angle_min
  max_angle = msg.angle_max
  delta_angle = msg.angle_increment
  
  ranges = []
  angles = []
  num_vals = 0

  for i in range(range_arr_length):
    if range_arr[i] != float('inf'):
      ranges.append(range_arr[i])
      angles.append(min_angle + i*delta_angle)
      num_vals = num_vals + 1
      
  cone_detections.num_detections = num_vals
  cone_detections.ranges = ranges
  cone_detections.angles = angles
  
  print(cone_detections)


rospy.init_node('lidar_detected_cones')
sub = rospy.Subscriber('/perception/laserscan', LaserScan, callback)
rospy.spin()
