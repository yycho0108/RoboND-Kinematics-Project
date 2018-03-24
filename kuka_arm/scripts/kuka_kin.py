#!/usr/bin/env python
"""
Handles KUKA Kinematics.
Authored by Yoonyoung Cho @ 03.19.2018
"""

import tf
from mpmath import *
from sympy import *
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

import pickle
import cloudpickle as pickle
#import dill

import rospkg
import os
rospack = rospkg.RosPack()
pkg_root = rospack.get_path('kuka_arm')

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
    def __init__(self, build=True):
        # Variables
        self._q_fk = symbols('q1:7') # for FK
        self._q_ik = symbols('r,p,y') # for IK

        # DH Parameters
        self._DH = [self.getDH(i) for i in range(7)]
        self._s = self._params(self._DH)

        # building takes time - build once, then load thenceforth.
        if build:
            T, T02, T_cor, R03, R03i, Rrpy = self._build(self._s)
            self._save(self._q_fk, self._q_ik, T, T02, T_cor, R03, R03i, Rrpy)
        T, T02, T_cor, R03, R03i, Rrpy = self._load()

        self._T = T
        self._T02 = T02
        self._T_cor = T_cor
        self._R03 = R03
        self._R03i = R03i
        self._Rrpy = Rrpy

        # validation
        #M = self._T(0,0,0,0,0,0)
        #print M[:3, 3]

    def _params(self, DH):
        """
        Numerical Values for DH Parameters;
        order : alpha, a, d, q
        a = (z_i - z_{i-1}) length along x_{i-1}
        d_i = (x_i - x_{i-1}) length along z_i
        alpha = angle(z_i, z_{i-1})
        follows coordinate definitions as shown in figures/coord.png
        """
        print('Constructing Numerical DH Parameters ...')
        s = {}
        s.update({k:v for (k,v) in zip(DH[0], [0, 0, 0.75])}) #0-1
        s.update({k:v for (k,v) in zip(DH[1], [-pi/2, 0.35, 0])})#1-2
        s.update({k:v for (k,v) in zip(DH[2], [0, 1.25, 0])})#2-3
        s.update({k:v for (k,v) in zip(DH[3], [-pi/2, -0.054, 1.50])})#3-4
        s.update({k:v for (k,v) in zip(DH[4], [pi/2, 0, 0])})#4-5
        s.update({k:v for (k,v) in zip(DH[5], [-pi/2, 0, 0])})#5-6
        s.update({k:v for (k,v) in zip(DH[6], [0, 0, 0.303, 0])})#6-7

        # handle exceptional case, angular offset
        q2 = symbols('q2')
        s[q2] = q2-pi/2
        return s

    def _build(self, s):
        # Transformation Matrix
        print('Constructing Raw Transformation Matrix ...')
        T_raw = [self.dh2T(*dh) for dh in self._DH]
        print('Substituting Constants ...')
        T_par = [T.subs(s) for T in T_raw]

        def mdprint(M):
            print ' == BEGIN == '
            #print '|||'
            #print ':-:|:-:|:-:|:-:'
            print '\\begin{pmatrix}'
            n,_ = M.shape
            for i in range(n):
                s = ' & '.join(['{}'.format(e) for e in M[i,:]]) + ' \\\\'
                print s
            print '\\end{pmatrix}'

        for T in T_par:
            mdprint(T)

        print('Composing Homogeneous Transforms ...')
        T = reduce(lambda a,b:simplify(a*b), T_par) # full transform
        print('Applying Final Correction ...')
        T_cor = self.dh2URDF()
        T = T*T_cor
        mdprint(T)
        T02 = T_par[0] * T_par[1] #0->2
        R03 = (T02 * T_par[2])[:3,:3]
        #R03i = simplify(R03.inv("LU"))
        R03i = simplify(R03.T) # pretty sure this should work.
        print R03i

        r, p, y = symbols('r,p,y')
        Rrpy = rmat('z',y) * rmat('y',p) * rmat('x',r) * self.dh2URDF(R=True)

        return T, T02, T_cor, R03, R03i, Rrpy

    def _save(self, syms_fk, syms_ik, T, T02, T_cor, R03, R03i, Rrpy):
        fname = os.path.join(pkg_root, 'config', 'TF.txt') 

        # TEST: arguments validation
        #print 'Args'
        #print T.free_symbols #q[1-6]
        #print T02.free_symbols #q[1-2]
        #print T_cor.free_symbols #None
        #print R03i.free_symbols #q[1-3]
        #print Rrpy.free_symbols #r,p,y

        # create lambda functions
        f_T = lambdify(syms_fk, T)
        f_T02 = lambdify(syms_fk[:2], T02)
        f_T_cor = lambdify([], T_cor)
        f_R03 = lambdify(syms_fk[:3], R03)
        f_R03i = lambdify(syms_fk[:3], R03i)
        f_Rrpy = lambdify(syms_ik, Rrpy)

        # Dump
        pickle.dump([f_T, f_T02, f_T_cor, f_R03, f_R03i, f_Rrpy], open(fname, 'w'))

    def _load(self):
        fname = os.path.join(pkg_root, 'config', 'TF.txt') 
        f_T, f_T02, f_T_cor, f_R03, f_R03i, f_Rrpy = pickle.load(open(fname, 'r'))
        return f_T, f_T02, f_T_cor, f_R03, f_R03i, f_Rrpy

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
        T = self._T(*q)
        #T = self._T.subs({k:v for (k,v) in zip(self._q_fk, q)})
        T = np.array(T).astype(np.float32)
        pos = tf.transformations.translation_from_matrix(T)
        rot = tf.transformations.euler_from_matrix(T)
        return pos, rot

    def IK(self, pos, rot):
        """ Compute Inverse Kinematics """

        # inverse position ...

        # compute wrist position ...
        r, p, y = rot
        #Rrpy = rmat('z',y) * rmat('y',p) * rmat('x',r) * self.dh2URDF(R=True)
        #Print 'rrpy0', Rrpy
        Rrpy = self._Rrpy(r,p,y)
        #print 'rrpy1', Rrpy
        n = np.asarray(Rrpy[:, 2]).astype(np.float32) # normal
        #print n

        d6 = 0.0
        l = 0.303
        wpos = np.subtract(pos, (d6 + l)*n.reshape([3]))
        #print 'wpos', wpos

        q1 = np.arctan2(wpos[1], wpos[0])
        p2 = self._T02(q1, 0.0)[:3,3]
        #p2 = self._T02.subs({'q1':q1})[:3,3]
        # convert to float; TODO: better way?
        p2 = (float(p2[0]), float(p2[1]), float(p2[2]))

        dx, dy, dz = np.subtract(wpos, p2)

        def nrmat(a):
            ca = np.cos(a)
            sa = np.sin(a)
            return np.reshape([ca,-sa,sa,ca], (2,2))

        dr = nrmat(-q1).dot([dx, dy])[0]
        # don't want to care abt projections
        #dr = np.dot(n[:2], [dx,dy])#n.dot([dx,dy,dz])
        #dr = np.sqrt(float(dx*dx)+float(dy*dy))

        # TODO : avoid hardcoding constants
        # refer to diagram in #15 for a,b,c assignments
        r_c = 1.25 #a2
        r_a = np.sqrt(1.50**2 + 0.054**2) #d4
        r_b = np.sqrt(dr*dr+dz*dz)

        a = coslaw(r_b, r_c, r_a)
        b = coslaw(r_c, r_a, r_b)
        q2 = np.pi/2 - a - np.arctan2(dz, dr)
        #q3 = np.pi/2 - b - 0.03619# not exactly aligned
        q3 = np.pi/2 - b - np.arctan2(0.054, 1.50)
        #sol = cxc(r_c, [dr, dz], r_a)

        #if sol is None:
        #    # fallback to q2-q3??
        #    # TODO : validate
        #    print 'Fail - Unreachable'
        #else:
        #    sol_i = np.argmax([sol[0][1], sol[1][1]])
        #    #sol = sol[sol_i] # "prefer" elbow up
        #    ##sol = sol[0] # choose one
        #    sol0_q2 = np.arctan2(sol[0][0], sol[0][1])
        #    sol1_q2 = np.arctan2(sol[1][0], sol[1][1])

        #    sol0_q3 = np.arctan2(sol[0][1]-dz, dr-sol[0][0]) - q2 - 0.03619
        #    sol1_q3 = np.arctan2(sol[1][1]-dz, dr-sol[1][0]) - q2 - 0.03619

        #    i = np.argmin(np.abs([sol0_q2, sol1_q2]))

        #    q2 = [sol0_q2, sol1_q2][i]
        #    q3 = [sol0_q3, sol1_q3][i]
        #    ##print q2, [sol0_q2, sol1_q2][i]
        #    #print q3, sol0_q3, sol1_q3


        # inverse rotation ...
        #R36 = self._T_par[3]*self._T_par[4]*self._T_par[5]*self._T_par[6]
        #R36 = R36[:3,:3]
        #print 'R36', simplify(R36)
        #[
        #[-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - sin(q6)*cos(q4)*cos(q5), -sin(q5)*cos(q4)],
        #[sin(q5)*cos(q6), -sin(q5)*sin(q6), cos(q5)],
        #[-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4), sin(q4)*sin(q6)*cos(q5) - cos(q4)*cos(q6), sin(q4)*sin(q5)]
        #])

        R03i = self._R03i(q1, q2, q3)
        #R03 = self._R03(q1,q2,q3)
        #print R03
        #R03i = np.linalg.inv(R03)

        #R03i = self._R03i.subs({'q1':q1, 'q2':q2, 'q3':q3})
        R36 = np.array(np.dot(R03i, Rrpy)).astype(np.float32)

        q4 = np.arctan2(R36[2,2], -R36[0,2])
        q6 = np.arctan2(-R36[1,1], R36[1,0])
        q5 = np.arctan2(-R36[1,1]/np.sin(q6), R36[1,2])

        return q1,q2,q3,q4,q5,q6 #q4,q5,q6

def hide_axis(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

def error(kin, n=1024):
    """ Characterize FK-IK Errors """
    fig = plt.figure(figsize=(9.6, 4.8))
    print fig.get_size_inches()
    

    # dummy axis for setting super title
    ax0 = fig.add_subplot(111)
    ax0.set_title('IK Error Characterization')
    hide_axis(ax0)

    ax_p = fig.add_subplot(121, projection='3d')
    ax_p.set_title('Position (m) ')
    ax_q = fig.add_subplot(122, projection='3d')
    ax_q.set_title('Orientation (rad) ')

    ps = []
    perrs = []

    qs = []
    qerrs = []

    cm0= matplotlib.cm.get_cmap('cool')
    cm = cm0(np.arange(cm0.N))
    cm[:,-1] = np.linspace(0.1, 1, cm0.N) # minimum 0.1
    cm = ListedColormap(cm)

    for i in range(n):
        p = np.random.uniform(-2.0, 2.0, size=3)
        q = np.random.uniform(-np.pi, np.pi, size=3)
        q[1] = np.random.uniform(-np.pi/2, np.pi/2) # limit pitch to +-pi/2
        #q *= 0

        ik = kin.IK(p, q)
        fk = kin.FK(ik)
        if np.any(np.isnan(fk)):
            # impossible
            continue
        perr = np.subtract(fk[0], p) # positional error, meters
        qerr = np.subtract(fk[1], q) # rotational error, radians

        ps.append(p)
        qs.append(q)
        perrs.append(perr)
        qerrs.append(qerr)

    ps, perrs, qs, qerrs = [np.float32(e) for e in (ps,perrs,qs,qerrs)]

    ax = ax_p
    c = np.linalg.norm(perrs, axis=-1)
    s = ax.scatter(ps[:,0], ps[:,1], ps[:,2],
            c=c,
            cmap=cm#'rainbow'
            )
    fig.colorbar(s, ax=ax)
    
    ax.quiver(ps[:,0], ps[:,1], ps[:,2],
            perrs[:,0], perrs[:,1], perrs[:,2],
            color='red',
            alpha=0.5,
            )
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax = ax_q
    #c = np.linalg.norm(qerrs, axis=-1)
    #c = np.sqrt(np.sum(np.square(qerrs), axis=-1))
    #print np.sort(c[np.logical_not(np.isnan(c))])[-100:]
    #print np.sort(qerrs[np.logical_not(np.isnan(qerrs))].ravel())[-100:]
    c = np.max(np.abs(qerrs), axis=-1)
    print np.sort(qerrs.ravel())[-100:]
    print np.sort(c)[-100:]
    s = ax.scatter(qs[:,0], qs[:,1], qs[:,2],
            c=c,
            cmap=cm,
            )

    fig.colorbar(s, ax=ax)
    ax.quiver(qs[:,0], qs[:,1], qs[:,2],
            qerrs[:,0], qerrs[:,1], qerrs[:,2],
            color='red',
            alpha=0.5,
            )
    ax.set_xlabel('r')
    ax.set_ylabel('p')
    ax.set_zlabel('y')

    #ax_q.view_init(elev=90, azim=0)
    #ax_p.view_init(elev=90, azim=0)

    def init():
        ax_q.view_init(azim=0)
        ax_p.view_init(azim=0)
        return fig,

    def animate(i):
        print i
        ax_q.view_init(azim=i)
        ax_p.view_init(azim=i)
        return fig,

    #ani = FuncAnimation(fig, animate, init_func=init,
    #        frames=range(0,360,2), interval=20, blit=False)
    #ani.save('ik_errors.gif', fps=20, writer='imagemagick')

    plt.show()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = np.linalg.norm(perrs, axis=-1)
    s = ax.scatter(ps[:,0], ps[:,1],
            c=c,
            cmap=cm#'rainbow'
            )
    fig.colorbar(s, ax=ax)
    
    #ax.quiver(ps[:,0], ps[:,1],
    #        perrs[:,0], perrs[:,1],
    #        color='red',
    #        alpha=0.5,
    #        angles='xy',
    #        scale_units='xy',
    #        scale=1.0,
    #        )
    t = np.linspace(-np.pi, np.pi)
    ax.plot(
            0.303 + 0.35*np.cos(t),
            0.35*np.sin(t))
    print np.max(perrs[:,0])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.grid()
    plt.show()


def test(kin, n=1024, lim=np.pi, tol=np.deg2rad(1.0)):
    np.random.seed(0)

    good_if = [] # ik->fk
    bad_if = []

    good_fi = [] #fk->ik
    bad_fi = []

    for i in range(n):
        if (i%100) == 0:
            print('{}/{}'.format(i, n))

        # test fk-ik
        q = np.random.uniform(-lim, lim, size=6)
        fk = kin.FK(q)
        ik = kin.IK(*fk)
        if np.allclose(q, ik, rtol=1e-3, atol=tol):
            good_fi.append(q)
        else:
            bad_fi.append(q)

        # test ik-fk
        p = np.random.uniform(-2, 2, size=3) # xyz
        q = np.random.uniform(-np.pi, np.pi, size=3) #rpy
        q[1] = np.random.uniform(-np.pi/2, np.pi/2)
        ik = kin.IK(p,q)
        fk = kin.FK(ik)
        if np.allclose(p, fk[0], atol=0.01) and np.allclose(q, fk[1], atol=tol):
            good_if.append(list(p) + list(q))
        else:
            bad_if.append(list(p) + list(q))

    good_fi = np.asarray(good_fi, dtype=np.float32)
    bad_fi  = np.asarray(bad_fi, dtype=np.float32)

    good_if = np.asarray(good_if, dtype=np.float32)
    bad_if  = np.asarray(bad_if, dtype=np.float32)

    # draw ...
    fig = plt.figure()

    # fk-ik
    good = good_fi
    bad = bad_fi

    ax = fig.add_subplot(221, projection='3d')
    ax.set_title('Joint Space Viz')
    if len(good) > 0:
        ax.scatter(good[:,0], good[:,1], good[:,2],  c='b', label='good012')
    if len(bad) > 0:
        ax.scatter(bad[:,0], bad[:,1], bad[:,2], c='r', label='bad012')
    ax.set_xlabel('0')
    ax.set_ylabel('1')
    ax.set_zlabel('2')
    ax.legend()

    ax = fig.add_subplot(223, projection='3d')
    if len(good) > 0:
        ax.scatter(good[:,3], good[:,4], good[:,5],  c='b', label='good345')
    if len(bad) > 0:
        ax.scatter(bad[:,3], bad[:,4], bad[:,5], c='r', label='bad345')
    ax.set_xlabel('3')
    ax.set_ylabel('4')
    ax.set_zlabel('5')
    ax.legend()

    # ik-fk
    good = good_if
    bad = bad_if
    ax = fig.add_subplot(222, projection='3d')
    ax.set_title('Cartesian Space Viz')
    if len(good) > 0:
        ax.scatter(good[:,0], good[:,1], good[:,2],  c='b', label='goodxyz')
    if len(bad) > 0:
        ax.scatter(bad[:,0], bad[:,1], bad[:,2], c='r', label='badxyz')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.legend()

    ax = fig.add_subplot(224, projection='3d')
    if len(good) > 0:
        ax.scatter(good[:,3], good[:,4], good[:,5],  c='b', label='goodrpy')
    if len(bad) > 0:
        ax.scatter(bad[:,3], bad[:,4], bad[:,5], c='r', label='badrpy')
    ax.set_xlabel('r')
    ax.set_ylabel('p')
    ax.set_zlabel('y')
    ax.legend()

    plt.show()

def main():
    kin = KUKAKin(build=False)

    #for n in range(100):
    #    r,p,y = np.random.uniform(-np.pi, np.pi, size=3)
    #    r0 = np.linalg.inv(kin._R03(r,p,y))
    #    r1 = kin._R03i(r,p,y)
    #    print np.allclose(r0, r1)

    #r = kin._R03(0,0,1)
    #print r
    #print np.linalg.inv(r)
    #print kin._R03i(0,0,1)

    #test(kin, n=4096, lim=np.pi)
    # verified : usually ok for -pi/2 < pitch < pi/2
    error(kin, n=1024)

    #xs = []
    #ys = []
    #for r in np.linspace(-np.pi, np.pi, 200):
    #    ik = kin.IK([2.153, 0, 1.947], [0, 0, r])
    #    fk = kin.FK(ik)[1][2]
    #    print r
    #    xs.append(r)
    #    ys.append(fk)

    #plt.plot(xs,ys)
    #plt.show()

    #res = kin.IK([2.149, 0.027, 1.906], [1.465, 0.135, 0.091])
    #print res
    #print kin.FK(res)

if __name__ == "__main__":
    main()
