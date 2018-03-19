#!/usr/bin/env python
# Handles KUKA Kinematics.
import tf
from mpmath import *
from sympy import *
import numpy as np

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
        self._T02 = self._T_par[0] * self._T_par[1]*T_cor #0->2

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
        Ry = rmat('y', -np.pi)
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

    def getDH(self, i):
        return symbols('alpha{0}, a{0}, d{1}, q{1}'.format(i, i+1))

    def FK(self, q):
        return self._T.subs({k:v for (k,v) in zip(self._q, q)})

    def IK(self, pos, rot):
        #d6 = 
        # inverse position ...
        #x,y,z = pos
        r,p,y = rot

        # compute wrist position ...
        # wx = x - (d6+l)*nx ...
        Rrpy = rmat('z',y) * rmat('y',p) * rmat('x',r) * self.dh2URDF(R=True)
        print 'rrpy', Rrpy
        print Rrpy * Matrix([0,0,1]) #z-vec rot.

        n0 = np.asarray(Rrpy[:, 0]).astype(np.float32)
        n1 = np.asarray(Rrpy[:, 1]).astype(np.float32)
        n = np.asarray(Rrpy[:, 2]).astype(np.float32)
        #n = np.asarray([1,0,0], dtype=np.float32)

        d6 = 0.0
        l = 0.303
        wpos = np.subtract(pos, (d6 + l)*n0.ravel()) # NOT [:,2]??
        #wpos1 = np.subtract(pos, (d6 + l)*n1.ravel())
        #wpos2 = np.subtract(pos, (d6 + l)*n.ravel())
        print 'wpos', wpos

        q1 = np.arctan2(wpos[1],wpos[0])
        p2 = self._T02.subs({symbols('q1'):q1})[:3,3]
        print 'p2', p2
        # convert to float; TODO: better way?
        p2 = (float(p2[0]), float(p2[1]), float(p2[2]))

        dx, dy, dz = np.subtract(wpos, p2)
        print 'dz', dz
        # don't want to care abt projections
        dr = np.sqrt(float(dx*dx)+float(dy*dy))

        # TODO : avoid hardcoding constants
        r_c = 1.25 #a2
        r_a = np.sqrt(1.50**2 + 0.054**2) #d4

        p3s = cxc(r_c, [dr,dz], r_a) #wx-py2, wz-py2?

        if p3s is None:
            # not possible
            return None

        # set preference - elbow up
        if p3s[0][1] > p3s[1][1]:
            p3 = p3s[1]
        else:
            p3 = p3s[0]

        print 'delta', dr, dz
        print 'p3', p3

        q2 = np.arctan2(p3[0], p3[1])
        q3 = np.arctan2(p3[1]-dz, dr-p3[0]) - 0.03619

        # inverse rotation ...
        #r,p,y = rot

        return q1,q2,q3,0,0,0#q4,q5,q6

def main():
    kin = KUKAKin()
    #print kin.FK([0,0,0,0,0,0])
    #print kin.IK([1.498, 0.0, 1.160], [0,0,0])
    #q1, q2, q3 = kin.IK([2.153, 0, 1.947], [0, 0, 0.0])[:3]
    q1, q2, q3 = kin.IK([-1.001, -0.728, 1.479], [-0.740, -0.367, -1.894])[:3]
    print 'q1-q2-q3', q1,q2,q3

    #t = np.array(kin.FK([q1, q2, q3, 0, 0, 0])[:3,3]).astype(np.float32)
    #t = t.ravel()
    #q1, q2, q3 = kin.IK(t, [0,0,0])[:3]
    #print 'q1-q2-q3', q1,q2,q3

if __name__ == "__main__":
    main()
