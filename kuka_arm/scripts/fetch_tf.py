#!/usr/bin/env python2

import rospy
import tf

def main():
    rospy.init_node('fetch_tf')
    tl = tf.TransformListener()
    rate = rospy.Rate(10)
    while not rospy.is_shutdown():
        try:
            t, q = tl.lookupTransform('base_link', 'link_2', rospy.Time())
            _, _, d1 = t
            t, q = tl.lookupTransform('link_1', 'link_2', rospy.Time())
            a1, _, _= t
            t, q = tl.lookupTransform('link_2', 'link_3', rospy.Time())
            _, _, a2 = t
            t, q = tl.lookupTransform('link_3', 'link_5', rospy.Time())
            d4, _, a3 = t
            t, q = tl.lookupTransform('link_5', 'gripper_link', rospy.Time())
            d7,_,_ = t

            print '==='
            print 'd1 : {}'.format(d1)
            print 'a1 : {}'.format(a1)
            print 'a2 : {}'.format(a2)
            print 'd4 : {}'.format(d4)
            print 'a3 : {}'.format(a3)
            print 'd7 : {}'.format(d7)

        except tf.Exception as e:
            rospy.logwarn_throttle(0.5, 'TF Exception : {}'.format(e))
            pass

        rate.sleep()

if __name__ == "__main__":
    main()
