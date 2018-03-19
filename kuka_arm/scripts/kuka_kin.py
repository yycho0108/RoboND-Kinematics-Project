#!/usr/bin/env python
"""
Handles KUKA Kinematics.
Authored by Yoonyoung Cho @ 03.19.2018
"""

import tf
from mpmath import *
from sympy import *
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def rmat(axis, angle):
    """ 3x3 Rotation Matrix from Axis{x,y,z} + Angle(rad) """
    c, s = cos(angle), sin(angle)

    if axis == 'x':
        return Matrix([
            [1,0,0],
            [0,c,-s],
            [0,s,c]])
    elif axis == 'y':
        return Matrix([
            [c,0,s],
            [0,1,0],
            [-s,0,c]])
    elif axis == 'z':
        return Matrix([
            [c,-s,0],
            [s,c,0],
            [0,0,1]])

    raise NotImplementedError("Generalized rmat is not implemented yet!")
    return None

def coslaw(a,b,c):
    """ Cosine Law """
    return np.arccos((a*a+b*b-c*c) / (2*a*b))

def cxc(r1,c2,r2):
    """
    Intersection of two circles,
    First Circle located at origin with radius r1,
    Second Circle located at c2=(x2,y2) with radius r2.
    returns None if the circles do not intersect.
    """
    x,y = c2
    r11 = r1*r1
    r22 = r2*r2
    R2 = x*x+y*y 
    det = 2*(r11+r22)/R2 - np.square(r11-r22)/(R2*R2) - 1

    if det < 0:
        return None

    c = np.multiply(0.5+(r11-r22)/(2*R2), c2)
    d = np.multiply(0.5*np.sqrt(det), [y,-x])
    return c+d, c-d

class KUKAKin(object):
    """
    Forward & Inverse Kinematics
    for KUKA KR210 6DOF Robot Arm.
    """
    def __init__(self):

        # DH Parameters
        self._DH = [self.getDH(i) for i in range(7)]

        # Numerical Values for DH Parameters
        print('Constructing Numerical DH Parameters ...')
        self._s = {}
        # order : alpha, a, d, q

        # a = (z_i - z_{i-1}) length along x_{i-1}
        # d_i = (x_i - x_{i-1}) length along z_i
        # alpha = angle(z_i, z_{i-1})
        # follows coordinate definitions as shown in figures/coord.png

        self._s.update({k:v for (k,v) in zip(self._DH[0], [0, 0, 0.75])}) #0-1
        self._s.update({k:v for (k,v) in zip(self._DH[1], [-pi/2, 0.35, 0])})#1-2
        self._s.update({k:v for (k,v) in zip(self._DH[2], [0, 1.25, 0])})#2-3
        self._s.update({k:v for (k,v) in zip(self._DH[3], [-pi/2, -0.054, 1.50])})#3-4
        self._s.update({k:v for (k,v) in zip(self._DH[4], [pi/2, 0, 0])})#4-5
        self._s.update({k:v for (k,v) in zip(self._DH[5], [-pi/2, 0, 0])})#5-6
        self._s.update({k:v for (k,v) in zip(self._DH[6], [0, 0, 0.303, 0])})#6-7

        # handle exceptional case, angular offset
        q2 = symbols('q2')
        self._s[q2] = q2-pi/2

        # Transformation Matrix
        print('Constructing Raw Transformation Matrix ...')
        self._T_raw = [self.dh2T(*dh) for dh in self._DH]
        print('Substituting Constants ...')
        self._T_par = [T.subs(self._s) for T in self._T_raw]
        print('Composing Homogeneous Transforms ...')
        self._T = reduce(lambda a,b:simplify(a*b), self._T_par) # full transform
        print('Applying Final Correction ...')
        T_cor = self.dh2URDF()
        self._T = self._T*T_cor
        self._T02 = self._T_par[0] * self._T_par[1] #0->2
        self._R03 = (self._T02 * self._T_par[2])[:3,:3]

        # Variables
        self._q = symbols('q1:7')
        # == self._T.free_variables()

        # validation
        #M = self._T.subs({q:0 for q in self._q})
        #print M[:, 3]

    @staticmethod
    def dh2URDF(R=False):
        """ Correction from DH-URDF """
        Rz = rmat('z', np.pi)
        Ry = rmat('y', -np.pi/2)
        Rc = simplify(Rz*Ry)
        if R:
            return Rc
        T = Rc.row_join(Matrix([0,0,0])).col_join(Matrix([[0,0,0,1]]))
        return T

    @staticmethod
    def dh2T(alpha, a, d, q):
        """ Convert DH Parameters to Transformation Matrix """
        cq = cos(q)
        sq = sin(q)
        ca = cos(alpha)
        sa = sin(alpha)

        T = Matrix([
            [cq, -sq, 0, a],
            [sq*ca, cq*ca, -sa, -sa*d],
            [sq*sa, cq*sa, ca, ca*d],
            [0, 0, 0, 1]
            ])
        return T

    @staticmethod
    def getDH(i):
        """ Yield DH Parameters for a joint at index i."""
        return symbols('alpha{0}, a{0}, d{1}, q{1}'.format(i, i+1))

    def FK(self, q):
        """ Compute Forward Kinematics """
        # TODO : support homogeneous -> trans, rot (by option)
        T = self._T.subs({k:v for (k,v) in zip(self._q, q)})
        T = np.array(T).astype(np.float32)
        pos = tf.transformations.translation_from_matrix(T)
        rot = tf.transformations.euler_from_matrix(T)
        return pos, rot

    def IK(self, pos, rot):
        """ Compute Inverse Kinematics """

        # inverse position ...

        # compute wrist position ...
        r, p, y = rot
        Rrpy = rmat('z',y) * rmat('y',p) * rmat('x',r) * self.dh2URDF(R=True)
        n = np.asarray(Rrpy[:, 2]).astype(np.float32) # normal

        d6 = 0.0
        l = 0.303
        wpos = np.subtract(pos, (d6 + l)*n.ravel())
        #print 'wpos', wpos

        q1 = np.arctan2(wpos[1], wpos[0])
        p2 = self._T02.subs({'q1':q1})[:3,3]
        # convert to float; TODO: better way?
        p2 = (float(p2[0]), float(p2[1]), float(p2[2]))

        dx, dy, dz = np.subtract(wpos, p2)
        # don't want to care abt projections
        dr = np.sqrt(float(dx*dx)+float(dy*dy))

        # TODO : avoid hardcoding constants
        # refer to diagram in #15 for a,b,c assignments
        r_c = 1.25 #a2
        r_a = np.sqrt(1.50**2 + 0.054**2) #d4
        r_b = np.sqrt(dr*dr+dz*dz)

        a = coslaw(r_b, r_c, r_a)
        b = coslaw(r_c, r_a, r_b)

        q2 = np.pi/2 - a - np.arctan2(dz, dr)
        q3 = np.pi/2 - b - 0.03619# not exactly aligned

        # inverse rotation ...
        #R36 = self._T_par[3]*self._T_par[4]*self._T_par[5]*self._T_par[6]
        #R36 = R36[:3,:3]
        #print 'R36', simplify(R36)
        #[
        #[-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - sin(q6)*cos(q4)*cos(q5), -sin(q5)*cos(q4)],
        #[sin(q5)*cos(q6), -sin(q5)*sin(q6), cos(q5)],
        #[-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4), sin(q4)*sin(q6)*cos(q5) - cos(q4)*cos(q6), sin(q4)*sin(q5)]
        #])

        R03 = self._R03.subs({'q1':q1, 'q2':q2, 'q3':q3})
        R36 = np.array(R03.inv("LU")*Rrpy).astype(np.float32)

        q4 = np.arctan2(R36[2,2], -R36[0,2])
        q6 = np.arctan2(-R36[1,1], R36[1,0])
        q5 = np.arctan2(-R36[1,1]/np.sin(q6), R36[1,2])

        # TODO : characterize when numbers are unstable

        return q1,q2,q3,q4,q5,q6 #q4,q5,q6

def test(kin, n=256, tol=np.deg2rad(1.0)):

    good = []
    bad = []

    for i in range(n):
        m = len(good)+len(bad)
        if (m%100) == 0:
            print('{}/{}'.format(m, n))

        q = np.random.uniform(-np.pi, np.pi, size=6)
        fk = kin.FK(q)
        ik = kin.IK(*fk)
        if np.allclose(q, ik, rtol=1e-3, atol=tol):
            good.append(fk[0])
        else:
            bad.append(fk[0])

    good = np.asarray(good, dtype=np.float32)
    bad  = np.asarray(bad, dtype=np.float32)

    print np.shape(good)
    print np.shape(bad)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    if len(good) > 0:
        ax.scatter(good[:,0], good[:,1], good[:,2],  c='b')
    if len(bad) > 0:
        ax.scatter(bad[:,0], bad[:,1], bad[:,2], c='r')
    plt.show()

def main():
    kin = KUKAKin()
    test(kin)
    #ik = kin.IK([2.153, 0, 1.947], [0, 0, 0.0])
    #print kin.FK(ik)

    #res = kin.IK([2.149, 0.027, 1.906], [1.465, 0.135, 0.091])
    #print res
    #print kin.FK(res)

if __name__ == "__main__":
    main()
