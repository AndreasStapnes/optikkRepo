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
    fig, ax = plt.subplots(1,1)
    w0 = 3.25e15
    wavelen0 = 2*np.pi*3e8/w0
    delta_w = 0.5e15
    uni_dist = UniformFreqDist(w0, delta_w)
    s = np.linspace(-1.2e-5, 1.2e-5, 1000)
    tau = s/3e8
    I = uni_dist.Gamma(1, 0*s, True)*(1+uni_dist.gamma(tau, True)) / delta_w
    ax.plot(s, I, label="I")
    ax.set_ylabel("I_0")
    ax.set_xlabel(r"$\Delta x$")
    null_len = wavelen0*w0/delta_w
    customline(ax, x=null_len)
    customline(ax, x=-null_len)
    customline(ax, y=2)
    print(null_len)
    #customline(ax, y=0)
    ax.legend()
    plt.show()
    fig.savefig("eksempelLongitudalInterference.svg")

