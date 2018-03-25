#!/usr/bin/env python

"""
Handles KUKA Kinematics.
Authored by Yoonyoung Cho @ 03.25.2018

Requirements:
    - All related ROS Libraries
    - sympy, numpy (math)
    - pickle, cloudpickle (cache)
    - matplotlib, scipy (visualization)
"""

# math ...
import tf
from mpmath import *
from sympy import *
import numpy as np

# visualization ...
import matplotlib
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull

# cache ...
import pickle
import cloudpickle as pickle

# path configurations ...
import rospkg
import os

def rmat(axis, angle):
    """ 3x3 Symbolic Rotation Matrix from Axis{x,y,z} + Angle(rad) """
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

def nrmat(a):
    """ 2x2 Numeric Rotation Matrix about Z-Axis """
    ca = np.cos(a)
    sa = np.sin(a)
    return np.reshape([ca,-sa,sa,ca], (2,2))

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

def mdprint(M):
    """ Print Matrix in LaTeX-friendly format """
    print ' == BEGIN == '
    print '\\begin{pmatrix}'
    n,_ = M.shape
    for i in range(n):
        s = ' & '.join(['{}'.format(e) for e in M[i,:]]) + ' \\\\'
        print s
    print '\\end{pmatrix}'

class KUKAKin(object):
    """
    Forward & Inverse Kinematics
    for KUKA KR210 6DOF Robot Arm.
    """
    def __init__(self, build=True, fname=None):
        """
        Initialize Kinematics Utility;
        If running for the first time, build must be set to True.
        Otherwise, there is no need to re-run the computation;
        simply set build=False to enable loading from cached routines.

        Cache is, by default, built at "$(rospack find kuka_arm)/config/TF.txt".
        If other paths are desired, simply override fname and use it consistently, including IK_server.py.
        """
        # Variables
        self._q_fk = symbols('q1:7') # for FK
        self._q_ik = symbols('r,p,y') # for IK

        # DH Parameters
        self._DH = [self.getDH(i) for i in range(7)]
        self._s = self._params(self._DH)

        # Build/Save/Load FK/IK Routines
        if fname is None:
            # figure out default path
            rospack = rospkg.RosPack()
            pkg_root = rospack.get_path('kuka_arm')
            fname = os.path.join(pkg_root, 'config', 'TF.txt') 

        self._fname = fname

        if not os.path.exists(self._fname) or build:
            # build if cache is not initialized, or forced by build=True
            T, T02, T_cor, R03, R03i, Rrpy = self._build(self._s)
            self._save(self._q_fk, self._q_ik, T, T02, T_cor, R03, R03i, Rrpy)
        T, T02, T_cor, R03, R03i, Rrpy = self._load()

        self._T = T
        self._T02 = T02
        self._T_cor = T_cor
        self._R03 = R03
        self._R03i = R03i
        self._Rrpy = Rrpy

        # remember useful numbers
        self._a2, self._a3, self._d4, self._d6, self._d7 = \
            [float(e.subs(self._s)) for e in symbols('a2, a3, d4, d6, d7')]

    def _params(self, DH):
        """
        Numerical Values for DH Parameters;
        order : alpha, a, d, q
        a = (z_i - z_{i-1}) length along x_{i-1}
        d_i = (x_i - x_{i-1}) length along z_i
        alpha = angle(z_i, z_{i-1})
        follows coordinate definitions as shown in figures/coord.png
        """
        # Define DH Parameters
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
        """ Build Kinematics Routines """
        print('Constructing Raw Transformation Matrix ...')
        T_raw = [self.dh2T(*dh) for dh in self._DH]
        print('Substituting Constants ...')
        T_par = [T.subs(s) for T in T_raw]

        print('Composing Homogeneous Transforms ...')
        T = reduce(lambda a,b:simplify(a*b), T_par) # full transform
        print('Applying Final Correction ...')
        T_cor = self.dh2URDF()
        T = T*T_cor
        T02 = T_par[0]*T_par[1] #0->2
        R03 = (T02*T_par[2])[:3,:3] # extract rotation part from p2
        #R03i = simplify(R03.inv("LU"))
        R03i = simplify(R03.T) # transpose of rotation matrix is its inverse

        # inverse kinematics end-effector orientation
        r, p, y = symbols('r,p,y')
        Rrpy = rmat('z',y) * rmat('y',p) * rmat('x',r) * self.dh2URDF(R=True)

        return T, T02, T_cor, R03, R03i, Rrpy

    def _save(self, syms_fk, syms_ik, T, T02, T_cor, R03, R03i, Rrpy):
        """ Save kinematics routines from _build() """
        # create lambda functions
        f_T = lambdify(syms_fk, T)
        f_T02 = lambdify(syms_fk[:2], T02)
        f_T_cor = lambdify([], T_cor)
        f_R03 = lambdify(syms_fk[:3], R03)
        f_R03i = lambdify(syms_fk[:3], R03i)
        f_Rrpy = lambdify(syms_ik, Rrpy)

        # Dump
        pickle.dump([f_T, f_T02, f_T_cor, f_R03, f_R03i, f_Rrpy], open(self._fname, 'w'))

    def _load(self):
        """ Load kinematics routines from _save() """
        f_T, f_T02, f_T_cor, f_R03, f_R03i, f_Rrpy = pickle.load(open(self._fname, 'r'))
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
        """
        Compute Forward Kinematics.
        Parameters:
            q : Joint Angles, (q1,q2,q3,q4,q5,q6)
        Returns:
            (x,y,z), (r,p,y) : End effector position
        """
        T = self._T(*q)
        T = np.array(T).astype(np.float32)
        pos = tf.transformations.translation_from_matrix(T)
        rot = tf.transformations.euler_from_matrix(T)
        return pos, rot

    def IK(self, pos, rot):
        """
        Compute Inverse Kinematics.
        Parameters:
            (x,y,z), (r,p,y) : End effector position
        Returns:
            q : Joint Angles, (q1,q2,q3,q4,q5,q6)
        """

        # unpack params
        a2, a3, d4, d6, d7 = self._a2, self._a3, self._d4, self._d6, self._d7

        # compute wrist position ...
        r, p, y = rot
        Rrpy = self._Rrpy(r,p,y)
        n = np.asarray(Rrpy[:, 2]).astype(np.float32) # extract normal-z

        wpos = np.subtract(pos, (d6 + d7)*n.reshape([3]))

        # obtain q1 from Trigonometry
        q1 = np.arctan2(wpos[1], wpos[0])
        p2 = self._T02(q1, 0.0)[:3,3]
        p2 = (float(p2[0]), float(p2[1]), float(p2[2]))

        dx, dy, dz = np.subtract(wpos, p2)

        # IMPORTANT : SIGNED DISTANCE,
        # NOT sqrt(dx**2 + dy**2).
        dr = nrmat(-q1).dot([dx, dy])[0]

        # obtain q2/q3 from cosine laws.
        # Refer to [diagram](figures/q2q3.png) for a,b,c assignments
        r_a = np.sqrt(a3**2 + d4**2)
        r_b = np.sqrt(dr*dr+dz*dz)
        r_c = a2

        a = coslaw(r_b, r_c, r_a)
        b = coslaw(r_c, r_a, r_b)
        q2 = np.pi/2 - a - np.arctan2(dz, dr)
        q3 = np.pi/2 - b + np.arctan2(a3, d4) # account for angle offset

        # =================================
        # Option 2 : Circular Intersections.
        # Deprecated; do not use.
        # left here only for archival purposes.

        # sol = cxc(r_c, [dr, dz], r_a)
        # if sol is None:
        #     # fallback to q2-q3??
        #     # TODO : validate
        #     print 'Fail - Unreachable'
        # else:
        #     sol_i = np.argmax([sol[0][1], sol[1][1]])
        #     #sol = sol[sol_i] # "prefer" elbow up
        #     ##sol = sol[0] # choose one
        #     sol0_q2 = np.arctan2(sol[0][0], sol[0][1])
        #     sol1_q2 = np.arctan2(sol[1][0], sol[1][1])
        #     sol0_q3 = np.arctan2(sol[0][1]-dz, dr-sol[0][0]) - q2 - 0.03619
        #     sol1_q3 = np.arctan2(sol[1][1]-dz, dr-sol[1][0]) - q2 - 0.03619
        #     i = np.argmin(np.abs([sol0_q2, sol1_q2]))
        #     q2 = [sol0_q2, sol1_q2][i]
        #     q3 = [sol0_q3, sol1_q3][i]
        #     ##print q2, [sol0_q2, sol1_q2][i]
        #     #print q3, sol0_q3, sol1_q3
        # ====================================

        # inverse rotation ...

        # Reference:
        # R36 = self._T_par[3]*self._T_par[4]*self._T_par[5]*self._T_par[6]
        # R36 = R36[:3,:3]
        # print 'R36', simplify(R36)
        # should yield:
        # [
        # [-sin(q4)*sin(q6) + cos(q4)*cos(q5)*cos(q6), -sin(q4)*cos(q6) - sin(q6)*cos(q4)*cos(q5), -sin(q5)*cos(q4)],
        # [sin(q5)*cos(q6), -sin(q5)*sin(q6), cos(q5)],
        # [-sin(q4)*cos(q5)*cos(q6) - sin(q6)*cos(q4), sin(q4)*sin(q6)*cos(q5) - cos(q4)*cos(q6), sin(q4)*sin(q5)]
        # ])

        # ... Rrpy = Rz(y)*Ry(p)*Rx(r)*Rc
        R03i = self._R03i(q1, q2, q3) # inverse of R03
        R36 = np.array(np.dot(R03i, Rrpy)).astype(np.float32)

        # obtain q4,q5,q6 as corresponding to the above reference
        q4 = np.arctan2(R36[2,2], -R36[0,2])
        q6 = np.arctan2(-R36[1,1], R36[1,0])
        q5 = np.arctan2(-R36[1,1]/np.sin(q6), R36[1,2])

        return q1,q2,q3,q4,q5,q6

def hide_axis(ax):
    """ Hide Axis (Dummy) """
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

def characterize(kin, n=1024, animate=False):
    """ Characterize FK-IK Errors """
    fig = plt.figure(figsize=(9.6, 4.8))

    # dummy axis for setting super title
    ax0 = fig.add_subplot(111)
    ax0.set_title('IK Error Characterization')
    hide_axis(ax0)

    # create axes, Layout [ 1 | 2 ]
    ax_p = fig.add_subplot(121, projection='3d')
    ax_p.set_title('Position (m) ')
    ax_q = fig.add_subplot(122, projection='3d')
    ax_q.set_title('Orientation (rad) ')

    # keep track of positions, orientations, and errors
    ps = []
    qs = []
    perrs = []
    qerrs = []

    # also keep track of unreachable positions
    bad_ps = []

    # create new color map with alpha values from 0.1-1.0
    cm0= matplotlib.cm.get_cmap('cool')
    cm = cm0(np.arange(cm0.N))
    cm[:,-1] = np.linspace(0.1, 1, cm0.N) # minimum 0.1
    cm = ListedColormap(cm)

    # run ...
    for i in range(n):
        p = np.random.uniform(-4.0, 4.0, size=3) # +-4m
        q = np.random.uniform(-np.pi, np.pi, size=3)
        q[1] = np.random.uniform(-np.pi/2, np.pi/2) # limit pitch to +-pi/2

        ik = kin.IK(p, q)
        fk = kin.FK(ik)
        if np.any(np.isnan(fk)):
            bad_ps.append(p)
            # impossible
            continue
        perr = np.subtract(fk[0], p) # positional error, meters
        qerr = np.subtract(fk[1], q) # rotational error, radians

        ps.append(p)
        qs.append(q)
        perrs.append(perr)
        qerrs.append(qerr)

    ps, perrs, qs, qerrs = [np.float32(e) for e in (ps,perrs,qs,qerrs)]
    bad_ps = np.float32(bad_ps)

    hull = ConvexHull(ps)

    print '== Workspace Boundary : =='
    print 'Max(XYZ) : {}'.format(np.max(ps, axis=0))
    print 'Min(XYZ) : {}'.format(np.min(ps, axis=0))
    print '=========================='

    # plot ...
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
    ax.plot_trisurf(ps[:,0],ps[:,1],ps[:,2],
            triangles=hull.simplices, color='green',
            alpha=0.1
            )
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax = ax_q

    c = np.max(np.abs(qerrs), axis=-1)
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

    if animate:
        def init():
            ax_q.view_init(azim=0)
            ax_p.view_init(azim=0)
            return fig,#ax_q, ax_p

        def animate(i):
            print i
            ax_q.view_init(azim=i)
            ax_p.view_init(azim=i)
            return fig,#ax_q, ax_p

        plt.draw()
        plt.pause(0.001)
        ani = FuncAnimation(fig, animate, init_func=init,
                frames=range(0,360,2), interval=20, blit=False)
        # uncomment below to save animation
        #ani.save('ik_errors.gif', fps=20, writer='imagemagick')
    plt.show()

    # also show workspace boundary
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    h_b = ax.scatter(bad_ps[:,0], bad_ps[:,1], bad_ps[:,2],
            c='red', alpha=0.05,
            label='bad'
            )
    h_g = ax.scatter(ps[:,0], ps[:,1], ps[:,2],
            c='blue', alpha=0.1,
            label='good'
            )
    h_h = ax.plot_trisurf(ps[:,0],ps[:,1],ps[:,2],
            triangles=hull.simplices, color='green',
            alpha=0.2,
            label='hull'
            )
    h_h_fake = ax.scatter([0],[0],[0], linestyle="-", c='g', marker='o', alpha=0.2) 
    # add fake scatter to make labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('KR210 Workspace Characterization')
    ax.legend([h_b, h_g, h_h_fake], ['bad','good','hull'])
    plt.show()

def main():
    kin = KUKAKin(build=False)
    #test(kin, n=4096, lim=np.pi)
    #verified : usually ok for -pi/2 < pitch < pi/2
    characterize(kin, n=8192)

if __name__ == "__main__":
    main()
