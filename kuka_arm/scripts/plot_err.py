import numpy as np
from matplotlib import pyplot as plt

def hide_axis(ax):
    ax.spines['top'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['right'].set_color('none')
    ax.tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

def main():
    M = np.loadtxt('err1.csv')
    fig = plt.figure(figsize=(12.0, 4.8))
    ax0 = fig.add_subplot(111)
    ax0.set_title('Live Session Error', y=1.05)
    hide_axis(ax0)

    ax_p = fig.add_subplot(121)
    ax_p.plot(M[:,0], label='x')
    ax_p.plot(M[:,1], label='y')
    ax_p.plot(M[:,2], label='z')
    ax_p.grid()
    ax_p.set_xlabel('Request #')
    ax_p.set_ylabel('Error (m)')
    ax_p.legend()
    ax_p.set_title('Position')

    ax_q = fig.add_subplot(122)
    ax_q.plot(M[:,3], label='r')
    ax_q.plot(M[:,4], label='p')
    ax_q.plot(M[:,5], label='y')
    ax_q.grid()
    ax_q.set_xlabel('Request #')
    ax_q.set_ylabel('Error (rad)')
    ax_q.legend()
    ax_q.set_title('Orientation')
    #plt.legend(['x','y','z','r,','p','y'])
    #plt.title('Live Session Error')
    plt.show()

if __name__ == "__main__":
    main()
