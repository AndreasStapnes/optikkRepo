# Ønsker her å løse problemer av formen
# 2*pi/lambda * b * sqrt(n_1^2 - N^2) = arctan(sqrt(N^2-n_2^2/n_1^2-N^2)) + m*pi for naturlig m.
# Grunnet problemets stygge natur, vil jeg i hovedsak basere meg på bisection-method.
from typing import Callable, Tuple, List, Dict, Optional, NamedTuple
from enum import Enum
from numba import njit
from scipy import linalg
import numpy as np
from time import time

c = 2.99e8

'''denne filen er ikke laget for å kunne tydes av andre en meg.
Jeg har ikke tatt meg betydelig tid til å kommentere ting. For å forstå
metodikken er en nødt til å ha satt seg inn i problemet jeg ønsker å løse,
hvilket (litt vel kortfattet) er gjennomgått i en Optikk-forelesning.

Uansett er funksjonaliteten til programmet oppsummert i .pyi-filen, hvilket
danner grunnlaget for importering av denne løsningsmetoden.
'''


def carryBisect(end_points: Tuple[float, float],
                function: Callable[[float], float],
                tolerance: float,
                prev_endpoints_eval: Tuple[bool, bool]) -> float:
    mid_point = (end_points[0] + end_points[1])/2
    if end_points[1]-end_points[0] > 2*tolerance:
        mid_val = function(mid_point)
        if mid_val == 0:
            return mid_point
        mid_eval = mid_val > 0
        if mid_eval ^ prev_endpoints_eval[1]:
            next_points = (mid_point, end_points[1])
            next_evals = (mid_eval, prev_endpoints_eval[1])
        else:
            next_points = (end_points[0], mid_point)
            next_evals = (prev_endpoints_eval[0], mid_eval)

        return carryBisect(next_points, function, tolerance, next_evals)
    else:
        return mid_point


def bisection(end_points: Tuple[float, float],
              function: Callable[[float], float],
              tolerance: float = 1e-10) -> float:
    end_evals = (function(end_points[0]) > 0, function(end_points[1]) > 0)
    if not (end_evals[0] ^ end_evals[1]):
        raise Exception("Function does not change sign through interval")
    return carryBisect(end_points, function, tolerance, end_evals)


class Mode(Enum):
    TE: int = 0
    TM: int = 1


@njit()
def lhs_basis(n_1: float, n_2: float, wavelen: float, b: float, N: float, m: int) -> float:
    return 2*np.pi*b/wavelen*np.sqrt(n_1**2-N**2)


@njit()
def rhs_basis_TE(n_1, n_2, wavelen, b, N, m) -> float:
    if N == n_1: return (m+1)*np.pi
    else: return 2*np.arctan(np.sqrt((N**2-n_2**2)/(n_1**2-N**2))) + m*np.pi


@njit()
def rhs_basis_TM(n_1, n_2, wavelen, b, N, m) -> float:
    if N == n_1: return (m+1)*np.pi
    else: return 2*np.arctan(n_1**2/n_2**2*np.sqrt((N**2-n_2**2)/(n_1**2-N**2))) + m*np.pi


def refractive_index_bisection_method(n_1: float, n_2: float, wavelen: float,
                                      b: float, m: int, mode: Mode = Mode.TE) -> float:
    rhs_basis: Callable[[float, float, float, float, float, int], float] = \
        rhs_basis_TE if mode is Mode.TE else rhs_basis_TM

    def method(N: float):
        return rhs_basis(n_1, n_2, wavelen, b, N, m) - lhs_basis(n_1, n_2, wavelen, b, N, m)
    return bisection((n_2, n_1), method)


def find_effective_refractive_indexes(n_1: float, n_2: float, wavelen: float, b: float, mode:Mode=Mode.TE) -> List[float]:
    m_max = int(np.floor(lhs_basis(n_1, n_2, wavelen, b, n_2, 0)/np.pi))
    m_points = m_max + 1
    N_vals = np.zeros(m_points)
    for m in range(m_points):
        N_vals[m] = refractive_index_bisection_method(n_1, n_2, wavelen, b, m, mode)
    return N_vals




def find_ABCD(n_1: float, n_2: float, wavelen: float, b: float,
              mode: Mode = Mode.TE, N_sols: Optional[np.ndarray]=None):
    if N_sols is None:
        N_sols = find_effective_refractive_indexes(n_1, n_2, wavelen, b, mode)
    N_sols = np.array(N_sols)
    K_sols = np.sqrt(n_1**2-N_sols**2)
    gamma_sols = np.sqrt(N_sols**2-n_2**2)  #En faktor er utelatt her (2*np.pi/wavelen)
    Kb_sols = K_sols*b*2*np.pi/wavelen      #men dette har (overraskende nok) ingenting å si for videre beregninger
    m_points = len(N_sols)
    As, Bs, Cs, Ds = [], [], [], []
    for m in range(m_points):
        K = K_sols[m]
        gamma = gamma_sols[m]
        Kb = Kb_sols[m]
        qm_matrix = \
            np.array([[1,                   0,              -1,      0],
                      [0,                   K,              -gamma,  0],
                      [np.cos(Kb),          np.sin(Kb),      0,     -1],
                      [-K * np.sin(Kb),     K * np.cos(Kb),  0,      gamma]])\
            if mode is Mode.TE else \
            np.array([[1,                     0,                   -1,             0],
                      [0,                     K/n_1**2,            -gamma/n_2**2,  0],
                      [np.cos(Kb),            np.sin(Kb),           0,            -1],
                      [-K*np.sin(Kb)/n_1**2,  K*np.cos(Kb)/n_1**2,  0,             gamma/n_2**2]])
        # print(qm_matrix)
        # print(np.linalg.det(qm_matrix))
        sol = np.array(linalg.null_space(qm_matrix, rcond=1e-8))
        if np.shape(sol)[1] != 1:
            print(sol)
            print(np.shape(sol))
            raise Exception("Something went wrong during search for null-vectors")
        As.append(sol[0, 0])
        Bs.append(sol[1, 0])
        Cs.append(sol[2, 0])
        Ds.append(sol[3, 0])
    return As, Bs, Cs, Ds


class VectorArray:
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray

    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray):
        self.x, self.y, self.z = x, y, z
        assert len(self.x) == len(self.y) == len(self.z)

    def __getitem__(self, item:int):
        return self.x[item], self.y[item], self.z[item]

    def __mul__(self, other):
        return VectorArray(self.x*other, self.y*other, self.z*other)

class ElmagTensor(NamedTuple):
    E: VectorArray
    B: VectorArray

    def __mul__(self, other):
        return ElmagTensor(*[elem*other for elem in self])

    def __truediv__(self, other):
        return self*(1/other)


class ModeSolutionObject(NamedTuple):
    TE: ElmagTensor
    TM: ElmagTensor





class system:
    def __init__(self, n_1: float, n_2: float, wavelen: float, b: float):

        self.n_1 = n_1
        self.n_2 = n_2
        self.wavelen = wavelen
        self.b = b
        TE_N_sols = find_effective_refractive_indexes(n_1, n_2, wavelen, b, Mode.TE)
        TE_N_sols = np.array(TE_N_sols)
        self.TE_N_sols = TE_N_sols
        self.TE_K_sols = np.sqrt(n_1**2-TE_N_sols**2)*2*np.pi/wavelen
        self.TE_beta_sols = self.TE_N_sols*2*np.pi/wavelen
        self.TE_gamma_sols = np.sqrt(TE_N_sols**2-n_2**2)*2*np.pi/wavelen

        self.AsTE, self.BsTE, self.CsTE, self.DsTE = find_ABCD(n_1, n_2, wavelen, b, Mode.TE, TE_N_sols)

        TM_N_sols = find_effective_refractive_indexes(n_1, n_2, wavelen, b, Mode.TM)
        TM_N_sols = np.array(TM_N_sols)
        self.TM_N_sols = TM_N_sols
        self.TM_K_sols = np.sqrt(n_1 ** 2 - TM_N_sols ** 2) * 2 * np.pi / wavelen
        self.TM_beta_sols = self.TM_N_sols*2*np.pi/wavelen
        self.TM_gamma_sols = np.sqrt(TM_N_sols ** 2 - n_2 ** 2) * 2 * np.pi / wavelen

        self.AsTM, self.BsTM, self.CsTM, self.DsTM = find_ABCD(n_1, n_2, wavelen, b, Mode.TM, TM_N_sols)

    def n(self, x: np.ndarray):
        val = np.where((x <= 0) | (x >= self.b), self.n_1, self.n_2)
        return val

    def _evalABCD(self, A: float, B: float, C: float, D: float,
                  gamma: float, K: float) -> Callable[[np.ndarray], np.ndarray]:
        def eval(x: np.ndarray):
            vals = np.zeros_like(x)
            vals += np.where(x <= 0, C * np.exp(gamma * x), 0)
            vals += np.where((0 < x) & (x < self.b), A * np.cos(K * x) + B * np.sin(K * x), 0)
            vals += np.where(x >= self.b, D * np.exp(-gamma * (x - self.b)), 0)
            return vals
        return eval

    def _derivativeEvalABCD(self, A: float, B: float, C: float, D: float,
                            gamma: float, K: float) -> Callable[[np.ndarray], np.ndarray]:
        def eval(x: np.ndarray):
            vals = np.zeros_like(x)
            vals += np.where(x<=0, gamma * C * np.exp(gamma * x), 0)
            vals += np.where((0 < x) & (x < self.b), K * (B * np.cos(K*x) - A * np.sin(K*x)), 0)
            vals += np.where(x >= self.b, -gamma*D*np.exp(-gamma*(x-self.b)), 0)
            return vals
        return eval

    def getRank(self, m: int) -> Callable[[np.ndarray, Optional[float]], ModeSolutionObject]:
        A, B, C, D = self.AsTE[m], self.BsTE[m], self.CsTE[m], self.DsTE[m]
        gamma, K = self.TE_gamma_sols[m], self.TE_K_sols[m]
        TE_beta = self.TE_beta_sols[m]
        TE_E_y_method = self._evalABCD(A, B, C, D, gamma, K)
        TE_dE_y_dx_method = self._derivativeEvalABCD(A, B, C, D, gamma, K)
        A, B, C, D = self.AsTM[m], self.BsTM[m], self.CsTM[m], self.DsTM[m]
        gamma, K = self.TM_gamma_sols[m], self.TM_K_sols[m]
        TM_beta = self.TM_beta_sols[m]
        TM_B_y_method = self._evalABCD(A, B, C, D, gamma, K)
        TM_dB_y_dx_method = self._derivativeEvalABCD(A, B, C, D, gamma, K)

        w = 2*np.pi*c/self.wavelen

        def method(x: np.ndarray, t: Optional[float] = None):
            TE_E = VectorArray(np.zeros_like(x), TE_E_y_method(x), np.zeros_like(x))
            TE_B_x = -TE_beta/w * TE_E.y
            TE_B_y = np.zeros_like(x)
            TE_B_z = -1j/w*TE_dE_y_dx_method(x)
            TE_B = VectorArray(TE_B_x, TE_B_y, TE_B_z)
            TE = ElmagTensor(TE_E, TE_B) / np.max(np.abs(TE_E.y))

            TM_B = VectorArray(np.zeros_like(x), TM_B_y_method(x), np.zeros_like(x))
            TM_E_x = TM_beta*c**2/(w*self.n(x)**2) * TM_B.y
            TM_E_y = np.zeros_like(x)
            TM_E_z = 1j*c**2/(w*self.n(x)**2) * TM_dB_y_dx_method(x)
            TM_E = VectorArray(TM_E_x, TM_E_y, TM_E_z)
            TM = ElmagTensor(TM_E, TM_B) / np.max(np.abs(TM_B.y))

            if t is not None:
                TE = TE*np.exp(1j*w*t)
                TM = TM*np.exp(1j*w*t)

            return ModeSolutionObject(TE, TM)
        return method

    def __getitem__(self, m: int):
        return self.getRank(m)

