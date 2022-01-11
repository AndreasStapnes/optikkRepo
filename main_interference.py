import matplotlib.pyplot as plt
import numpy as np
from cosmetics import customline

class UniformFreqDist:
    w0: float
    deltaw: float

    def __init__(self, w0:float, deltaw:float):
        self.w0 = w0
        self.deltaw = deltaw

    def Gamma(self, W_0, t, real=False):
        val = 2*W_0*self.deltaw*np.sinc(t*self.deltaw/2/np.pi)*np.exp(-1j*self.w0*t)
        return val if not real else np.real(val)

    def gamma(self, t, real=False):
        val = np.sinc(t*self.deltaw/2/np.pi)*np.exp(-1j*self.w0*t)
        return val if not real else np.real(val)


if __name__ == "__main__":
    s = np.linspace(-0.6e-5, 0.6e-5, 1000)
    tau = s / 3e8

    fig, ax = plt.subplots(1,1)
    w0 = 3.25e15
    wavelen0 = 2*np.pi*3e8/w0
    print(f"{wavelen0=}")
    delta_w = 0.84e15
    uni_dist = UniformFreqDist(w0, delta_w)
    I = uni_dist.Gamma(1, 0*s, True)*(1+uni_dist.gamma(tau, True)) / delta_w
    ax.plot(s, I, label=r"$I(\Delta x)$")
    ax.set_ylabel("I_0")
    ax.set_xlabel(r"$\Delta x$ i meter")
    null_len = wavelen0*w0/delta_w
    customline(ax, x=null_len)
    customline(ax, x=-null_len)
    customline(ax, y=2)
    print(null_len)
    #customline(ax, y=0)
    ax.legend()
    plt.show()
    fig.savefig("eksempelLongitudalInterference.svg")

    fig, ax = plt.subplots(1,1)
    I = 2*(1+np.cos(w0*s/3e8))
    ax.plot(s, I, label=r"$I(\Delta x)$")
    ax.set_ylabel(r"$I_0$")
    ax.set_xlabel(r"$\Delta x$ i meter")
    ax.legend(loc="upper right")
    plt.show()
    fig.savefig("simpleLongitudalInterference.svg")

    fig, ax = plt.subplots(1,1)
    mod_uni_dist = UniformFreqDist(0, delta_w)
    sinc_modulation = mod_uni_dist.gamma(tau)
    cos_modulation = np.cos(w0*tau)
    ax.plot(s, sinc_modulation, label=r"$sinc\left(\frac{\Delta w\Delta x}{2c}\right)$")
    ax.plot(s, cos_modulation, label=r"$cos\left(\frac{w_0\Delta x}{c}\right)$", color="r", alpha=0.3)
    ax.legend(loc="upper right")
    customline(ax, y=0)
    customline(ax, x=null_len)
    customline(ax, x=-null_len)
    ax.set_xlabel(r"$\Delta x$ i meter")
    plt.show()
    fig.savefig("eksempelLongitudalInterferenceDelkomp.svg")

