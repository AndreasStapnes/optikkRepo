import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

example_n1 = 3
example_n2 = 1

if __name__ == '__main__':
    ax: Axes
    fig: Figure
    fig, ax = plt.subplots(1,1)
    xvals = np.linspace(-0.5,1,300)
    yvals = np.array([(lambda xval: -0.5 if xval > 0 and xval < 0.6 else -0.2)(xval) for xval in xvals])
    ax.plot(xvals, yvals, linewidth=3, color='r', label=r'$-k(x)^2$')
    ax.axvline(x=0, color='k', linestyle='--', alpha=0.5)
    ax.axvline(x=0.6, color='k', linestyle='--', alpha=0.5)
    ax.axhline(y=-0.5, color='k', linestyle='--', alpha=0.5)
    ax.axhline(y=-0.2, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel("x", fontsize=16)
    #ax.set_ylabel(r"$\frac{w^2n(x)^2}{c^2}$")
    ax.set_ylim([-0.6,0])
    #ax.tick_params(axis='x', label1On=False)
    #ax.tick_params(axis='y', label1On=False)
    ax.set_xticks([0, 0.6])
    ax.set_xticklabels([0, 'b'], fontsize=14)
    ax.set_yticks([-0.2, -0.5])
    ax.set_yticklabels([r'$-\frac{w^2n_2^2}{c^2}$', r'$-\frac{w^2n_1^2}{c^2}$'], fontsize=14)
    ax.legend()
    fig.show()
    fig.savefig('eksempelBronn.svg')