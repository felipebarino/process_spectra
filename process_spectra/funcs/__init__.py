"""
Esse módulo contêm as funções usadas no processo de análise, seja as de
carregamento, de extração ou filtragem do espectro
"""

import os
import numpy as np
from process_spectra.utils import lorentz
from scipy import signal as sg
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from process_spectra import utils


def load_from_optisystem(filename,
                         quiet=False,
                         delimiter=';',
                         wl_multiplier=1,
                         dtype=np.float64,
                         ignore_errors=False,
                         name_function=None):
    """
    Carrega o espectro e retorna o espectro (np array, comprimentos de onda e
    potência) e o dicionário com o nome extraído

    :param filename: O nome do arquivo do espectro
    :type filename: str

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :param delimiter: O delimitador que separa os dois valores no txt. ;
        por padrão (usado pelo optisystem)
    :type delimiter: str

    :param wl_multiplier: Um valor que multiplica todos os valores de
        comprimentos de onda. Deve ser um valor que faz os comprimentos serem
        vistos em metros
        Ex: para comprimentos salvos como nm (1e-9), esse
        valor deve ser 1e9.
        1 por padrão.
    :type wl_multiplier: float

    :param dtype: o tipo da variável a ser carregada dos arquivos. np.float64
        por padrão
    :type dtype: np.dtype

    :param ignore_errors: Se for False, vai parar o programa se não
        conseguir carregar algum espectro. Se for True, coloca os que não
        abriram como spectrum None, que é ignorado pelo programa
    :type ignore_errors: bool

    :param name_function: A função que vai identificar o nome. Deve receber
        o filename como argumento e retornar o nome. Deve tomar cuidado para
        que dois arquivos com nomes diferentes não gerem o mesmo nome na saída.
        Por padrão, o nome é extraído como o nome do arquivo sem extensões.
    :type name_function: function

    :return: um np array 2d com os valores (comprimentos de onda e
        potência) e o dicionário com o nome extraído
    :rtype: (np.ndarray, dict)
    """

    name_function = name_function or utils.remove_extension

    info = {'name': name_function(filename)}

    if not quiet:
        print(f'Carregando {info["name"]}')

    def conv(text):
        text = text.replace(b'E', b'e')
        text = text.replace(b',', b'.')
        return float(text)

    if ignore_errors:
        try:
            spectrum = np.loadtxt(filename, dtype=dtype,
                                  converters={0: conv, 1: conv},
                                  delimiter=delimiter)
        except ValueError:
            return None, info
    else:
        spectrum = np.loadtxt(filename, dtype=dtype,
                              converters={0: conv, 1: conv},
                              delimiter=delimiter)

    if spectrum[0, 0] > spectrum[-1, 0]:
        spectrum = spectrum[::-1]

    spectrum[:, 0] *= wl_multiplier

    return spectrum, info


def load_from_optigrating(filename, quiet=False):
    """Função para carregar o espectro de transmitância do arquivo do
    optigrating.

    :param filename: O nome do arquivo do espectro
    :type filename: str

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :return: um np array 2d com os valores (comprimentos de onda e
        potência) e o dicionário com o nome extraído
    :rtype: (np.ndarray, dict)
    """

    info = {'name': utils.remove_extension(filename)}

    if not quiet:
        print(f'Carregando {info["name"]}')

    complex_spectrum = np.genfromtxt(filename, dtype=np.float64, delimiter=' ')
    spectrum = np.zeros((complex_spectrum.shape[0], 2))

    spectrum[::, 0] = complex_spectrum[::, 0] * 1e-6
    spectrum[::, 1] = \
        (complex_spectrum[::, 1]**2 + complex_spectrum[::, 2]**2)**(1/2)

    spectrum[::, 1] = np.array(list(map(
        lambda x: 10 * np.log10(x), spectrum[::, 1]         # escalar pra dB
    )))

    return spectrum, info


def mask_spectrum(spectrum, _, wl_limits, quiet=False):
    """
    Corta o espectro com base em um intervalo de comprimento de onda

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param _: Ignorado. O programa entrega o info, e essa função não usa
    :type _: Any

    :param wl_limits: Uma tupla com 2 valores, contendo os limites do corte
    :type wl_limits: (float, float)

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :return: O espectro cortado e um dicionário vazio
    :rtype: (np.ndarray, dict)
    """

    if not quiet:
        print(f'Cortando')

    region = np.where((min(wl_limits) <= spectrum[::, 0]) &
                      (spectrum[::, 0] <= max(wl_limits)))
    masked = spectrum[region]

    return masked, dict()


def filter_spectrum(spectrum, _, window_length, polyorder, quiet=False):
    """
    Filtra o espectro utilizando o filtro de Savitzky–Golay

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param _: Ignorado. O programa entrega o info, e essa função não usa
    :type _: Any

    :param window_length: O tamanho da 'janela do filtro' (o número de
        coeficientes usados nos cálculos do filtro, deve ser ímpar)
    :type window_length: int

    :param polyorder: A ordem dos polinômios usados nos cálculos do filtro
        Deve ser menor do que o 'window_length'
    :type polyorder: int

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :return: O espectro após aplicação do filtro e um dicionátio vazio
    :rtype: (np.ndarray, dict)
    """

    if not quiet:
        print(f'Filtrando')

    filtered = spectrum * 1  # só pra copiar rápido
    filtered[::, 1] = sg.savgol_filter(spectrum[::, 1],
                                       window_length, polyorder)

    return filtered, dict()


def interpolate_spectrum(spectrum, _, wl_step, wl_limits=None,
                         kind='cubic', quiet=False):
    """
    Interpola o espectro para os valores de comprimento de onda dentro dos
    limites, considerando incrementos discretos de 'wl_step'

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param _: Ignorado. O programa entrega o info, e essa função não usa
    :type _: Any

    :param wl_step: O tamanho do incremento
    :type wl_step: float

    :param wl_limits: Os limites da região para interpolar. Usa os limites
        do espectro original como padrão
    :type wl_limits: (float, float)

    :param kind: Tipo de interpolação. Deve ser um desses: (‘linear’,
        ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’,
        ‘next’)
    :type kind: str

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :return: O espectro gerado pela interpolação e um dicionário vazio
    :rtype: (np.ndarray, dict)
    """
    if not quiet:
        print(f'Interpolando')

    wl_limits = wl_limits or (spectrum[0, 0], spectrum[-1, 0])
    wl = [x for x in np.arange(wl_limits[0], wl_limits[1], wl_step)]

    xs = spectrum[::, 0]
    ys = spectrum[::, 1]

    f = interp1d(xs, ys, kind=kind)
    interpolated = np.array(
        [x for x in zip(wl, f(wl))],
        dtype=np.float64)

    return interpolated, dict()


def simulate_gain(spectrum, _, other, unit='dB'):
    """
    Simula o efeito da associação de dois espectros. Pode ser usado para
    obter o ganho / atenuação total de dois filtros em série (como uma LPG
    e uma FBG), ou pra simular o espectro resultante após efeito de um
    filtro sobre a fonte, por exemplo

    :param spectrum: Um dos espectros
    :type spectrum: np.ndarray

    :param other: O outro espectro
    :type other: np.ndarray

    :param unit: A unidade dos espectros. dB por padrão
    :type unit: string

    :return: O espectro resultante
    :rtype: np.ndarray
    """

    if (spectrum[::, 0] != other[::, 0]).any():
        raise ValueError('Os espectros devem ter os wavelengths iguais. Por '
                         'favor interpole os dois com os mesmos parâmetros '
                         'para equalizar isso.')

    final_spectrum = spectrum * 1          # só pra copiar

    if unit == 'dB':
        final_spectrum[::, 1] += other[::, 1]
    elif unit == 'scalar':
        final_spectrum[::, 1] *= other[::, 1]
    else:
        raise ValueError('Unidade inválida. As implementadas são "dB" e '
                         '"scalar"')
    return final_spectrum, dict()


def find_valley(spectrum, _, prominence=5, ignore_errors=False, quiet=False):
    """
    Tenta achar o vale ressonante do espectro. Se achar, retorna as
    coordenadas no dicionário

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param _: Ignorado. O programa entrega o info, e essa função não usa
    :type _: Any

    :param prominence: A 'proeminência' ('altura do pico') mínima do vale.
    :type prominence: float

    :param ignore_errors: Se for False, vai parar o programa se não achar
        o vale. Se for True ignora e continua.
    :type ignore_errors: bool

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :return: O espectro original o dicionário com as coordenadas
    :rtype: (np.ndarray, dict)
    """

    if not quiet:
        print(f'Tentando achar o vale')

    xs = spectrum[::, 0]
    ys = spectrum[::, 1]
    valleys, properties = sg.find_peaks(-ys, prominence=prominence)

    info = dict()

    if len(valleys) < 1:
        if ignore_errors:
            info['resonant_wl'] = None
            info['resonant_wl_power'] = None
        else:
            raise Exception(f'Nenhum vale encontrado!')

    if len(valleys) > 1:
        print(f'Achou {len(valleys)} vales pro espectro, '
              f'retornando o com maior proeminência')

    best_match = np.argmax(properties['prominences'])
    x, y = xs[valleys[best_match]], ys[valleys[best_match]]

    info['resonant_wl'] = x
    info['resonant_wl_power'] = y

    return spectrum, info


def get_approximate_valley(spectrum, info, approx_func=lorentz, prominence=5,
                           resolution_proximity=2, p0=None, valley_samples=100,
                           plot=False):
    """
    Aproxima a região do vale como uma curva determinada na 'approx_func',
        depois extrai o comprimento de onda ressonante a partir da curva
        aproximada (funciona também para espectros com mais de um vale)

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param info: O dicionário com as informações (não é usado)
    :type info: dict

    :param approx_func: A função de approximação (lorentziana por padrão)
    :type approx_func: function

    :param prominence: A prominência mínima dos vales
    :type prominence: float

    :param resolution_proximity: A proximidade que o vale aproximado deve
        estar do observado sem aproximação (multiplicado pela resolução)
    :type resolution_proximity: float

    :param p0: Os parâmetros iniciais da aproximação. Deve ser uma lista com
        os parâmetros da função de aproximação (em ordem, ignorando x). Devem
        ser ajustados se a aproximação consistentemente falhar.
    :type p0: list

    :param valley_samples: A quantidade máxima de samples em um vale. Deve ser
        ajustado de acordo com a resolução do medidor
    :type valley_samples: int

    :param plot: Se deve ou não plotar os gráficos com aproximações. Se for
        True, o script deve estar rodando em um caminho que possui uma pasta
        chamada 'plots/'. Esse parâmetro serve para ajudar a ajustar os p0 e
        valley_samples
    :type plot: bool

    :return: O espectro original e o dicionário com os valores de comprimento
        de onda e potência extraídos
    :rtype: (np.ndarray, dict)
    """

    wl = spectrum[::, 0]
    power = spectrum[::, 1]
    resolution = np.mean(np.diff(wl))

    peaks, peak_info = sg.find_peaks(-power, prominence=prominence,
                                  plateau_size=0, wlen=valley_samples)

    _info = dict()

    if plot:
        fig, ax = plt.subplots()
        ax.plot(wl, power)

    for i in range(len(peaks)):
        valley = spectrum[peak_info['left_bases'][i]: peak_info['right_bases'][i], ::]

        p0 = p0 or [-peak_info['prominences'][i],
              wl[peaks[i]] * 1e6,
              1,
              power[peaks[i]]]

        try:
            popt, _ = curve_fit(approx_func, valley[::, 0] * 1e6, valley[::, 1],
                                p0=p0, max_nfev=10000,
                                bounds=(-np.inf, (0, np.inf, np.inf, np.inf)))

            resonant_wl = popt[1] * 1e-6
            resonant_power = approx_func(popt[1], *popt)

        except RuntimeError:
            resonant_wl = wl[peaks[i]]
            resonant_power = power[peaks[i]]

        if abs(resonant_wl - wl[peaks[i]]) > resolution_proximity * resolution:
            resonant_wl = wl[peaks[i]]
            resonant_power = power[peaks[i]]

        if len(peaks) == 1:
            _info['resonant_wl'] = resonant_wl
            _info['resonant_wl_power'] = resonant_power
        else:
            _info[f'resonant_wl_{i}'] = resonant_wl
            _info[f'resonant_wl_power_{i}'] = resonant_power

        if plot:
            ax.plot(valley[::, 0], approx_func(valley[::, 0] * 1e6, *popt))  # debug

    # Get the index to the best prominence. Can be useful for interrogation
    best_index = np.argmax(peak_info['prominences'])
    _info['best_index'] = best_index

    if plot:
        fig.savefig(f"plots/{info['name']}")  # debug

    return spectrum, _info


def get_max_power(spectrum, _):
    """
    Essa função pega o maior valor de potência do espectro

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param _: Ignorado. O programa entrega o info, e essa função não usa
    :type _: Any

    :return: O espectro e o dicionário com a potência máxima
    :rtype: (np.ndarray, dict)
    """
    max_value = spectrum[:, 1].max()
    info = {'max_power': max_value}

    return spectrum, info


def fill_name_zeros(spectrum, info, zeros=4):
    """
    Adiciona zeros na parte numérica do nome pra facilitar a ordenação

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param info: O dicionário com info original
    :type info: dict

    :param zeros: A quantidade de zeros pra preencher
    :type zeros:int

    :return: O espectro original e o dicionário com o novo nome
    :rtype: (np.ndarray, dict)
    """
    name = info['name']
    numbers = ''.join(x for x in name if x.isdigit())
    non_numbers = ''.join([x for x in name if not x.isdigit()])
    _info = {'name': non_numbers + numbers.zfill(zeros)}

    return spectrum, _info


def set_name_numbers(spectrum, info):
    """
    Coloca o nome como apenas os números no nome do arquivo

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param info: O dicionário como nome extraído
    :type info: dict

    :return: O espectro original e o dicionário com o novo nome
    :rtype: (np.ndarray, dict)
    """

    name = info['name']
    numbers = ''.join(x for x in name if x.isdigit())
    _info = {'name': numbers}

    return spectrum, info


def plot_spectrum(spectrum, info, plot_opts=None, subplots=None, save_folder=None,
                  quiet=False):
    """
    Plota o espectro em um gráfico, com o vale marcado com um x, se tiver
    sido encontrado

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param info: O dicionário com as informações
    :type info: dict

    :param plot_opts: Um dicionário com opções para plotar, se nenhum for
        passado, será usado um padrão definido na função.
    :type plot_opts: dict

    :param subplots: Os objetos criados pelo matplotlib para fazer um plot,
        retornados pelo 'plt.subplots()'. Geralmente salvos como (fig, ax).
    :type subplots: tuple

    :param save_folder: A pasta para salvar as imagens dos plots, se for
        None, não salva
    :type save_folder: str, path

    :param quiet: Se o programa deve printar sobre o progresso
    :type quiet: bool

    :return: None
    """
    if not quiet:
        print(f'Plotando')

    plot_opts = plot_opts or {
        'xlim': (spectrum[0, 0], spectrum[-1, 0]),
        'ylim': (-100, 0), 'animated': False, 'interval': 1e-6}

    if subplots is None:
        fig, ax = plt.subplots()
    else:
        fig, ax = subplots

    ax.set_xlim(plot_opts['xlim'])
    ax.set_ylim(plot_opts['ylim'])
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Optical power (dBm)')

    xs = spectrum[::, 0]
    ys = spectrum[::, 1]

    ax.plot(xs, ys)
    if 'resonant_wl' in info.keys():
        ax.plot(info['resonant_wl'], info['resonant_wl_power'], 'xk')

    if save_folder is not None:
        if not quiet:
            print(f'Salvando o plot')
        full_path = os.path.join(save_folder,
                                 f'_{info["name"]}.png')
        fig.savefig(full_path, transparent=False)

    return spectrum, dict()
