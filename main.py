import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from typing import Dict

from optiskLederSol import find_effective_refractive_indexes, Mode, system
from cosmetics import customline, level

# For å modifisere output burde i grunnen kun følgende variabler endres på

n_1 = 1.9           # (1)       brytnings index for optisk leder
n_2 = 1.0           # (1)       brytnings index for området rundt den optiske leder
wavelen = 400e-9    # (meter)   bølgelengde til aktuelle bølger (slik de hadde vært i vakuum)
b = 220e-9          # (meter)   optisk leder-flate tykkelse







elmag_sys = system(n_1, n_2, wavelen, b)


# lvls : Dict[Mode, List[float]] =
# {mode: find_effective_refractive_indexes(n_1, n_2, wavelen, b, mode) for mode in Mode}
lvls: Dict[Mode, np.ndarray] = {Mode.TE: elmag_sys.TE_N_sols,
                                Mode.TM: elmag_sys.TM_N_sols}

for datavalues in lvls.values():
    for index in range(len(datavalues)):
        datavalues[index] = -(2*np.pi*datavalues[index]/wavelen)**2


def to_min_k_sq(n, wavelen):
    return -(2*np.pi*n/wavelen)**2


if __name__ == '__main__':
    ax: Axes
    fig: Figure
    fig, ax = plt.subplots(1, 1)
    xvals = np.linspace(-0.5*b, 1.5*b, 300)
    yvals = np.array([(lambda xval: to_min_k_sq(n_1, wavelen) if xval > 0 and xval < b
            else to_min_k_sq(n_2, wavelen))(xval) for xval in xvals])
    ax.plot(xvals, yvals, linewidth=4, color='r', label=r'$-k(x)^2$')
    customline(ax, x=0)
    customline(ax, x=b)
    customline(ax, y=to_min_k_sq(n_1, wavelen))
    customline(ax, y=to_min_k_sq(n_2, wavelen))
    for lvl in lvls[Mode.TE]:
        level(ax, lvl, color="xkcd:nasty green")
    ax.plot([], [], linestyle=":", color="xkcd:nasty green", label=r"$-\beta^2_{TE}$")
    for lvl in lvls[Mode.TM]:
        level(ax, lvl, color="xkcd:burple")
    ax.plot([], [], linestyle=":", color="xkcd:burple", label=r"$-\beta^2_{TM}$")

    ax.set_xlabel("x", fontsize=16)
    ax.set_ylim((to_min_k_sq(n_1, wavelen)*1.2, 0.0))
    ax.set_xticks([0, b])
    ax.set_xticklabels([0, 'b'], fontsize=14)
    ax.set_yticks([to_min_k_sq(n_2, wavelen), to_min_k_sq(n_1, wavelen)])
    ax.set_yticklabels([r'$-\frac{w^2n_2^2}{c^2}$', r'$-\frac{w^2n_1^2}{c^2}$'], fontsize=14)
    ax.legend()
    fig.show()
    fig.savefig('eksempelBronn.svg')

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    lw = 1.4
    ls = "-"
    for m in range(len(elmag_sys.TE_N_sols)):
        TE_E_y = elmag_sys[m](xvals).TE.E.y
        TE_E_y *= np.sign(TE_E_y[0])
        TM_B_y = elmag_sys[m](xvals).TM.B.y
        TM_B_y *= np.sign(TM_B_y[0])
        ax1.plot(xvals, TE_E_y, color="xkcd:nasty green", linestyle=ls, linewidth=lw)
        ax2.plot(xvals, TM_B_y, color="xkcd:burple", linestyle=ls, linewidth=lw)

    ax1.plot([], [], color="xkcd:nasty green", linestyle=ls, linewidth=lw, label=r"$E_y$ for TE")
    ax2.plot([], [], color="xkcd:burple", linestyle=ls, linewidth=lw, label=r"$B_y$ for TM")
    customline(ax1, x=0);       customline(ax2, x=0)
    customline(ax1, x=b);       customline(ax2, x=b)
    customline(ax1, y=0);       customline(ax2, y=0)
    ax1.set_xticks([0, b]);              ax2.set_xticks([0, b])
    ax1.set_xticklabels(["0", "b"]);     ax2.set_xticklabels(["0", "b"])
    ax1.legend();                        ax2.legend()
    fig.show()