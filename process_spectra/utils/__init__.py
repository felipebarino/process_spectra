"""Esse módulo contêm algumas funções auxiliares"""


import numpy as np
import os
from scipy.interpolate import interp1d


def get_power(y, x=None, unit='dBm', noise=None):
    """
    Calcula a potência de um espectro (integral com np.trapz)

    :param y: O y da função (mais detalhes em np.trapz)
    :type y: np array

    :param x: O x da função (opcional, mais detalhes em np.trapz)
    :type x: np array

    :param unit: A unidade do espectro. dBm por padrão
    :type unit: string

    :param noise: O ruído elétrico na leitura da potência. O valor é usado
        como variância na geração de um aleatório
    :type noise: float

    :return: A potência calculada
    """

    if unit.upper() in ['DBM', 'DBMW']:
        y = dBmW_to_W(y)
    elif unit != 'W':
        raise ValueError(f'{unit} é uma unidade inválida, '
                         f'as implementadas são "dBm" e "W"')

    power = np.trapz(np.abs(y), x)          # A potência fica bem baixa,
                                            # na casa dos 1e-18

    if noise:
        power += np.random.randn()*noise

    return power


def dBmW_to_W(dBm):
    """
    Transforma de dBmW para W

    :param dBm: O valor em dBm (float ou np.ndarray)
    :type dBm: float ou np array

    :return: O valor em W
    """

    return 1e-3 * 10 ** (dBm/10)


def remove_extension(filename):
    """
    Retorna o nome base do arquivo sem extensão

    :param filename: O nome do arquivo
    :type filename: str

    :return: O nome extraído
    :rtype: str
    """
    return os.path.splitext(os.path.basename(filename))[0]


def interpolate(original, wl_limits, wl_step, kind='cubic'):
    """
    Interpola o array original entre os limites, com o passo especificado

    :param original: O array original
    :type original: np.array

    :param wl_limits: Os limites da interpolação
    :type wl_limits: tuple(float, float)

    :param wl_step: O passo de interpolação
    :type wl_step: float

    :param kind: O tipo de interpolação. Para mais detalhes, ver o kind do
        interp1d do scipy
    :type kind: str

    :return: O array interpolado
    :rtype: np.array
    """
    wl = [x for x in np.arange(wl_limits[0], wl_limits[1], wl_step)]
    f = interp1d(original[::, 0], original[::, 1], kind)

    final = np.array([x for x in zip(wl, f(wl))],
                     dtype=np.float64)

    return final


def gauss(x, a, x0, sigma, bias):
    """
    Retorna o valor da função gaussiana em x para os parâmetros

    :param x: Entrada (x)
    :type x: float, np.array

    :param a: Escala da função
    :type a: float

    :param x0: A média (Mu)
    :type x0: float

    :param sigma: O desvio padrão (Sigma)
    :type sigma: float

    :param bias: O deslocamento vertical
    :type bias: float

    :return: O valor no ponto específicado
    :rtype: Mesmo de x
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) + bias


def lorentz(x, a, x0, w, bias):
    """
    Retorna o valor da função Lorentziana em x para os parâmetros

    :param x: Entrada (x)
    :type x: float, np.array

    :param a: Escala da função
    :type a: float

    :param x0: A média (Mu)
    :type x0: float

    :param w: O FWHM (full width at half maximum) (Gama)
    :type w: float

    :param bias: O deslocamento vertical
    :type bias: float

    :return: O valor no ponto específicado
    :rtype: Mesmo de x
    """
    return a*(1 + ((x - x0)/(w/2))**2)**(-1) + bias
