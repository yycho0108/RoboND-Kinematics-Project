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
from kuka_kin import KUKAKin
import numpy as np

# Note : see kuka_kin.py for the bulk of the implementation!
kin = KUKAKin(build=False)
errs = []

def handle_calculate_IK(req):
    global kin, errs

    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print "No valid poses received"
        return -1
    else:
        #kin.FK()

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

            jpos = kin.IK([px,py,pz], [roll, pitch, yaw])
            fk = kin.FK(jpos)

            target = [px, py, pz, roll, pitch, yaw]
            forward = list(fk[0]) + list(fk[1])
            err = np.subtract(target, forward)
            errs.append(err)

            # Populate response for the IK request
            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
	    joint_trajectory_point.positions = jpos
	    joint_trajectory_list.append(joint_trajectory_point)

        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    save = rospy.get_param('~save', False)
    rospy.loginfo('Save Enabled : {}'.format(save))
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    rospy.loginfo('Ready to receive an IK request')
    rospy.on_shutdown(lambda: np.savetxt('/tmp/err.csv', np.float32(errs)) if save else None)
    rospy.spin()

    #rate = rospy.Rate(100)
    #while not rospy.is_shutdown():
    #    rate.sleep()

    #if save:
    #    np.savetxt('err.csv', np.float32(errs))

if __name__ == "__main__":
    IK_server()
