# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 2021

@author: Felipe Barino, Guilherme Sampaio

Esse pacote adiciona funções para simular fbgs acopladas a componentes
piezoelétricos, de forma que a vibração do componente modula a refletância
das fbgs.
"""

import numpy as np
from math import sin, pi
from process_spectra.utils import get_power
from process_spectra.utils.fbg import get_fbg_reflectance
from process_spectra.funcs import simulate_gain


def simulate_piezo_fbg(spectrum,
                       _,
                       frequency,
                       amplitude,
                       fbg_wl_bragg,
                       fbg_fwhm,
                       fbg_label='fbg',
                       total_time=None,
                       sample_rate=None,
                       electrical_noise=None,
                       sensitivity=6.475e-12,
                       max_harmonic_index=2,
                       quiet=False):
    """
    Essa função simula a atuação de uma fbg acoplada com um piezoelétrico ao
        longo do tempo. Vale notar que por enquanto a função é extremamente
        lenta. Se for usar ela em tempo real, deve pensar em otimizar

        |Artigo de onde foram pegos os dados do piezoelétrico:
    @article{dante2016temperature,
    title={A temperature-independent interrogation technique for FBG sensors using monolithic multilayer piezoelectric actuators},
    author={Dante, Alex and Bacurau, Rodrigo Moreira and Spengler, Anderson Wedderhoff and Ferreira, Elnatan Chagas and Dias, Jos{\'e} Ant{\^o}nio Siqueira},
    journal={IEEE Transactions on Instrumentation and Measurement},
    volume={65},
    number={11},
    pages={2476--2484},
    year={2016},
    publisher={IEEE}
    }

    :param spectrum: O espectro
    :type spectrum: np.ndarray

    :param _: Ignorado. O programa entrega o info, e essa função não usa
    :type _: Any

    :param frequency: A frequência da tensão no piezoelétrico
    :type frequency: float

    :param amplitude: A amplitude da tensão no piezoelétrico
    :type amplitude: float

    :param fbg_wl_bragg: O comprimento de onda ressonante da fbg
    :type fbg_wl_bragg: float

    :param fbg_fwhm: (full wifth at half maximum) A largura da fbg
    :type fbg_fwhm: float

    :param fbg_label: O nome do fbg no csv final, se tiver usando apenas um,
        pode deixar o padrão
    :type fbg_label: str

    :param total_time: A duração total da simulação em segundos. É 3 ciclos do
        sistema por padrão
    :type total_time: float

    :param sample_rate: A taxa de amostragem temporal, em Hz. É de
        (max_harmonic_index + 1) * frequency por padrão
    :type sample_rate: float

    :param electrical_noise: O ruído elétrico na leitura de potência da fbg.
        Se for float será usado como variância na geração do ruído com
        distribuição normal. Se for um np array vai ser a distribuição. 0 por
        padrão
    :type electrical_noise: float

    :param sensitivity: A sensitividade do piezoelétrico (tensão para
        deslocamento do vale ressonante da fbg acoplada). Apenas mude se souber
        exatamente o que está fazendo. O valor padrão foi pego de um artigo
        (menção acima)
    :type sensitivity: float

    :param max_harmonic_index: O harmônico máximo que vai ser salvo no
        espectro. 2 por padrão
    :type max_harmonic_index: int

    :param quiet: Se o programa deve printar o progresso
    :type quiet: bool

    :return: O espectro de entrada e o dicionário com as informações de
        magnitude e fase dos harmônicos pedidos
    :rtype: np.ndarray
    """

    if not quiet:
        print(f'Simulando fbg piezoelétrico: {fbg_label}')

    total_time = total_time or 4/frequency
    sample_rate = sample_rate or (max_harmonic_index + 1) * frequency * 2

    sim_time = np.arange(0, total_time, 1/sample_rate)

    wls = spectrum[::, 0]

    fbg_powers = []
    for c, t in enumerate(sim_time):
        piezo_voltage = amplitude * sin(2*pi*frequency*t)
        piezo_modulation = piezo_voltage * sensitivity

        reflectance = get_fbg_reflectance(fbg_wl_bragg + piezo_modulation,
                                          fbg_fwhm, wls)

        reflected_spectrum, _ = simulate_gain(spectrum, '', reflectance)
        reflected_spectrum = reflected_spectrum[::, 1]

        fbg_powers.append(
            get_power(reflected_spectrum, wls, noise=electrical_noise))

    fourier_frequencies = np.fft.fftfreq(len(fbg_powers), d=1/sample_rate)
    fourier_frequencies = np.fft.fftshift(fourier_frequencies)

    fourier_transform = np.fft.fft(np.array(fbg_powers))
    fourier_transform = np.fft.fftshift(fourier_transform) / len(fbg_powers)
    fourier_transform_mag = np.abs(fourier_transform)
    fourier_transform_phase = np.angle(fourier_transform, deg=True)

    info = dict()

    for i in range(0, max_harmonic_index + 1):
        harmonic_freq = i*frequency
        info[f'{fbg_label}_harmonic{i}_mag'] = \
            fourier_transform_mag[np.where(fourier_frequencies == harmonic_freq)][0]

        info[f'{fbg_label}_harmonic{i}_pahse'] = \
            fourier_transform_phase[np.where(fourier_frequencies == harmonic_freq)][0]

    return spectrum, info
