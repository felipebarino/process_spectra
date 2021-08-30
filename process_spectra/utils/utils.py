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

    if unit == 'dBm' or unit == 'dBmW':
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
    """Função para interpolar um np array. Deve ser passado pro pacote"""
    wl = [x for x in np.arange(wl_limits[0], wl_limits[1], wl_step)]
    f = interp1d(original[::, 0], original[::, 1], kind)

    final = np.array([x for x in zip(wl, f(wl))],
                     dtype=np.float64)

    return final
