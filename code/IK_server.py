#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *
import numpy as np


#Creating rotation matrices
def rot_z(q):
            
    R_z =  Matrix([[        cos(q),         -sin(q),             0,     0],
	           [        sin(q),          cos(q),             0,     0],
	           [             0,               0,             1,     0],
                   [             0,               0,             0,     1]])
    return R_z 
            
def rot_y(q):
            
    R_y =  Matrix([[        cos(q),               0,        sin(q),     0],
	           [             0,               1,             0,     0],
	           [       -sin(q),               0,        cos(q),     0],
	           [             0,               0,             0,     1]])
    return R_y  

def rot_x(q):
            
    R_x =  Matrix([[             1,               0,             0,     0],
    		   [             0,          cos(q),       -sin(q),     0],
	     	   [             0,          sin(q),        cos(q),     0],
	           [             0,               0,             0,     1]])
    return R_x  

# Create individual transformation matrices
def DH(q, alpha, d, a):
        
    T = Matrix([[             cos(q),            -sin(q),            0,              a],
                [ sin(q)*cos(alpha), cos(q)*cos(alpha), -sin(alpha), -sin(alpha)*d],
                [ sin(q)*sin(alpha), cos(q)*sin(alpha),  cos(alpha),  cos(alpha)*d],
                [                   0,                   0,            0,               1]])
        
    return T

def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:

    ### Your FK code here
    # Create symbols
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8') #theta
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')

        # Create Modified DH parameters
        s = {alpha0 :     0,  a0 :     0, d1 :  0.75,  q1:       q1,
             alpha1 : -pi/2,  a1 :  0.35, d2 :     0,  q2 : q2-pi/2,
             alpha2 :     0,  a2 :  1.25, d3 :     0,  q3:       q3,
             alpha3 : -pi/2,  a3 : -.054, d4 :  1.50,  q4:       q4,
             alpha4 :  pi/2,  a4 :     0, d5 :     0,  q5:       q5,
             alpha5 : -pi/2,  a5 :     0, d6 :     0,  q6:       q6,
             alpha6 :     0,  a6 :     0, d7 : 0.303,  q7 :       0}

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

	    # Extract end-effector position and orientation from request
	    # px,py,pz = end-effector position
	    # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                    req.poses[x].orientation.z, req.poses[x].orientation.w])

            ### Your IK code here
	    # Compensate for rotation discrepancy between DH parameters and Gazebo
            R_z = rot_z(pi)
            R_y = rot_y(-pi/2)
            R_corr = simplify(R_z * R_y)

	    R_x = rot_x(roll)
	    R_y = rot_y(pitch)
	    R_z = rot_z(yaw)
	    Rrpy = R_z * R_y * R_x * R_corr

            #Finding the wrist center
	    wx = px - (s[d6]+s[d7])*Rrpy[0,2]
	    wy = py - (s[d6]+s[d7])*Rrpy[1,2]
	    wz = pz - (s[d6]+s[d7])*Rrpy[2,2]  
	
	    # Calculate joint angles using Geometric IK method
            #Finding the joints rotation angles
	    theta1 = atan2(wy,wx)
#            theta1 = np.clip(theta1, -pi/2, pi/2)
   
            #Solving for theta 2 using trigonometry and the cosines law

	    rtheta2   = sqrt((sqrt((wx)**2+wy**2)-s[a1])**2+(wz-s[d1])**2)
	    r = sqrt(wx**2 + wy**2) - s[a1]
	    alphaIden = acos((rtheta2**2 + s[a2]**2 - s[d4]**2) / (2 * rtheta2*s[a2])) 
	    gammaIden = acos((s[d4]**2 + s[a2]**2 - rtheta2**2) / (2 * s[d4]*s[a2]))
	    betaIden  = acos((s[d4]**2 + rtheta2**2 - s[a2]**2) / (2 * rtheta2*s[d4]))

	    theta2 = pi / 2 -alphaIden - atan2(wz-s[d1],r)


            #Solving for theta 3 using the cosines law
	    d = sqrt(s[d4]**2 - s[a3]**2)
	    pimasgamma = acos(d/s[d4])

	    theta3 = pi / 2 - (gammaIden - pimasgamma)

	    #Homogeneous  transforms
	    #base_link to link_1
	    R0_1 = DH(q1, alpha0, d1, a0).subs(s)
            #link_1 to link_2
	    R1_2 = DH(q2, alpha1, d2, a1).subs(s)
	    #link_2 to link_3
	    R2_3 = DH(q3, alpha2, d3, a2).subs(s)
	    # Extract rotation matrices from the transformation matrices (Since R0_6 = Rrpy)
	    R0_2 =  simplify(R0_1 * R1_2) #base link to link2
	    R0_3 =  simplify(R0_2 * R2_3) #base link to link3

	    R0_3 = R0_3.evalf(subs={q1:theta1, q2:theta2, q3:theta3})
	    R3_6 = R0_3.inv("LU") * Rrpy

            #Solving for theta4 to theta6 from R3_6
	    if abs(R3_6[1,2]) != 1:
		theta5 = atan2(sqrt(R3_6[2,2]**2 + R3_6[0,2]**2), R3_6[1,2])
		if (theta1 > pi/2):
		    theta4 = -pi/2
		    theta6 = atan2(-R3_6[1,1], R3_6[1,0])  

		else:
		    theta4 = atan2(R3_6[2,2], -R3_6[0,2])
		    theta6 = atan2(-R3_6[1,1], R3_6[1,0])  
	   
	    elif R3_6[1,2] == 1:
	        theta5 = pi/2
		theta4 = 0
		theta6 = atan2(-R3_6[1,1], R3_6[1,0]) 

	    elif R3_6[1,2] == -1:
                theta5 = -pi/2
		theta4 = 0
		theta6 = atan2(R3_6[1,1], -R3_6[1,0]) 

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print "Ready to receive an IK request"
    rospy.spin()

if __name__ == "__main__":
    IK_server()
