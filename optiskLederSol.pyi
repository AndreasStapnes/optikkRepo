from enum import Enum
from typing import List, Callable, NamedTuple
import numpy as np

class Mode(Enum):
    TE: int
    TM: int

class VectorArray:
    '''
    En grunnleggende array som beskriver et felt (kun) definert langs en 1D-akse
    '''
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    def __init__(self, x: np.ndarray, y: np.ndarray, z: np.ndarray): ...
    def __getitem__(self, item: int): ...
    def __mul__(self, other): ...


class ElmagTensor(NamedTuple):
    '''
    Et tensor-objekt hvilket beskriver en fysisk mulig mode som består
    av et E-felt og et B-felt
    '''
    E: VectorArray
    B: VectorArray
    def __mul__(self, other): ...
    def __truediv__(self, other): ...

class ModeSolutionObject(NamedTuple):
    '''
    Et generelt objekt som angir både TE og TM-løsningene for en bestemt orden
    til et gitt system (orden kan være 0, 1, ..., som for QM partikkel i brønn)
    '''
    TE: ElmagTensor
    TM: ElmagTensor
    def __mul__(self, other): ...



def find_effective_refractive_indexes(n_1: float, n_2: float, wavelen: float, b: float, mode:Mode=Mode.TE) -> List[float]:
    '''
    metode som returnerer en list effektive refraktive indekser som løser optisk-leder-problemet med de
    aktuelle parametrene. Resultatet er en liste indekser.

    :param n_1: Refraktiv indeks inni den optiske leder-flaten
    :param n_2: Refraktiv indeks utfor den optiske leder-flaten
    :param wavelen: Bølgelengden den aktuelle elektromagnetiske bølgen ville hatt i vakuum
    :param b: Tykkelsen på leder-flaten i x-retning
    :param mode: Mode.TM eller Mode.TE avhengig av hvilken mode en ønsker å se på
    :return: Liste av floats som angir approksimative løsninger for effektive refraktive indekser som oppfyller problemets kriterier
    '''
    ...

class system:
    '''
    Et objekt som fullstendig beskriver et optisk leder-system slik behandlet i
    optikk-faget. Et system kan initialiseres med et sett beskrivende parametre.
    Objektet vil deretter beregne alle mulige effektive refraktive indekser for TE og TM
    (som kan hentes i system.TE_N_sols og .TM_N_sols)

    Man vil i tillegg kunne hente ut felt-informasjon om alle ordener av løsninger for
    både TE og TM. Operasjonen system[index](x-np-array) vil returnere et
    ModeSolutionObject som beskriver både TE og TM-moden i de gitte x-posisjonene
    for orden 'index'.
    '''
    TE_N_sols: np.ndarray
    TM_N_sols: np.ndarray
    def __init__(self, n_1: float, n_2: float, wavelen: float, b: float):
        ...

    def __getitem__(self, m: int) -> Callable[[np.ndarray], ModeSolutionObject]:
        ...