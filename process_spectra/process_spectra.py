# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:58:09 2020

@author: Felipe Barino
"""

import os
import numpy as np
import pandas as pd
from scipy import signal as sg
from copy import deepcopy
from scipy.interpolate import interp1d
from matplotlib import pyplot as plt


class SpectrumData:
    """
    Uma classe usada para processar as informações de um espectro, com o
    objetivo final de extrair o comprimento de onda e atenuação do vale
    ressonante.
    """

    def __init__(self, filename, quiet=False, name=None):
        """
        Inicia o objeto, criando umas variáveis necessárias e carregando o
        espectro.

        :param filename: O nome do arquivo do espectro (saída do optisystem)
        :type filename: str

        :param quiet: Se True, o programa deixa de dar mensagens a cada etapa,
            False por padrão.
        :type quiet: bool, optional

        :param name: O nome do espectro. Usado para motivos de salvar e debug.
            É o filename, removendo a extensão, por padrão
        :type name: str, optional
        """
        self.quiet = quiet
        self.info = {
            'name': name or filename[0:-4],
        }

        self.resonant_wl = None
        self.resonant_wl_power = None
        self.spectrum = self._load_spectrum(filename)

    def _load_spectrum(self, filename,
                       delimiter=';',
                       wl_multiplier=1,
                       dtype=np.float64):
        """
        Carrega o espectro e retorna um np array 2d com os valores
        (comprimentos de onda e potência)

        :param filename: O nome do arquivo do espectro
        :type filename: str

        :param delimiter: O delimitador que separa os dois valores no txt. ;
            por padrão (usado pelo optisystem)
        :type delimiter: str, optional

        :param wl_multiplier: Um valor que multiplica todos os valores de
            comprimentos de onda. Deve ser um valor que faz os comprimentos serem
            vistos em metros.
            Ex: para comprimentos salvos como nm (10e-9), esse
            valor deve ser 10e9.
            1 por padrão.
        :param dtype: float, optional

        :return: um np array 2d com os valores (comprimentos de onda e
            potência)
        """
        if not self.quiet:
            print(f'Carregando {self.info["name"]}')

        def conv(text):
            text = text.replace(b'E', b'e')
            text = text.replace(b',', b'.')
            return float(text)

        spectrum = np.loadtxt(filename, dtype=dtype,
                              converters={0: conv, 1: conv}, delimiter=';')

        if spectrum[0, 0] > spectrum[-1, 0]:
            spectrum = spectrum[::-1]

        return spectrum

    def set_additional_info(self, info):
        """
        Adiciona informações adicionais ao espectro, como medidas, por exemplo

        :param info: Um dicionário com as informações a serem adicionadas.
            Ex: { 'periodo': 400 , 'temp': 20, 'strain': 200, 'mode': 4 }
        :type info: dict

        :return: None
        """
        for key in info.keys():
            self.info[key] = info[key]

    def mask_spectrum(self, wl_limits):
        """
        Corta o espectro com base em um intervalo de comprimento de onda

        :param wl_limits: Uma tupla com 2 valores, contendo os limites do corte
        :type wl_limits: tuple

        :return: O espectro (objeto dessa classe) cortado
        """
        if not self.quiet:
            print(f'Cortando {self.info["name"]}')

        region = np.where(min(wl_limits) <= self.spectrum[0] <= max(wl_limits))

        masked = deepcopy(self)
        masked.spectrum = masked.spectrum[region]

        return masked

    def filter_spectrum(self, window_length, polyorder):
        """
        Filtra o espectro utilizando o filtro de Savitzky–Golay

        :param window_length: O tamanho da 'janela do filtro' (o número de
            coeficientes usados nos cálculos do filtro.
        :type window_length: int, deve ser ímpar

        :param polyorder: A ordem dos polinômios usados nos cálculos do filtro
        :type polyorder: int, deve ser menor do que o 'window_length'

        :return: O espectro (objeto dessa classe) após aplicação do filtro.
        """
        if not self.quiet:
            print(f'Filtrando {self.info["name"]}')

        filtered = deepcopy(self)

        print(window_length)

        filtered.spectrum[::, 1] = sg.savgol_filter(
            filtered.spectrum[::, 1], window_length, polyorder)

        return filtered

    def interpolate_spectrum(self, wl_step, wl_limits=None,
                             kind='cubic'):
        """
        Interpola o espectro para os valores de comprimento de onda dentro dos
        limites, considerando incrementos discretos de 'wl_step'.

        :param wl_step: O tamanho do incremento.
        :type wl_step: float

        :param wl_limits: Os limites da região para interpolar. Usa os limites
            do espectro original como padrão.
        :type wl_limits: tuple with 2 floats

        :param kind: Tipo de interpolação. Deve ser um desses: (‘linear’,
            ‘nearest’, ‘zero’, ‘slinear’, ‘quadratic’, ‘cubic’, ‘previous’,
            ‘next’).
        :type kind: str

        :return: O espectro (objeto dessa classe) gerado pela interpolação.
        """
        if not self.quiet:
            print(f'Interpolando {self.info["name"]}')

        wl_limits = wl_limits or (self.spectrum[0, 0], self.spectrum[-1, 0])
        wl = [x for x in np.arange(wl_limits[0], wl_limits[1], wl_step)]

        xs = self.spectrum[::, 0]
        ys = self.spectrum[::, 1]

        f = interp1d(xs, ys, kind=kind)

        interpolated = deepcopy(self)
        interpolated.spectrum = np.array(
            [x for x in zip(wl, f(wl))],
            dtype=np.float64)

        return interpolated

    def find_valley(self, prominence=5):
        """
        Tenta achar o vale ressonante do espectro. Se achar, salva no
        dicionário 'info' objeto.

        :param prominence: A 'proeminência' ('altura do pico') mínima do vale.
        :type prominence: float

        :return: As coordenadas do vale (comprimento de onda, potência)
        """
        if not self.quiet:
            print(f'Tentando achar o vale do {self.info["name"]}')

        xs = self.spectrum[::, 0]
        ys = self.spectrum[::, 1]
        valleys, properties = sg.find_peaks(-ys, prominence=prominence)

        if len(valleys) < 1:
            raise Exception(f'No valleys found for {self.info["name"]}.')

        if len(valleys) > 1:
            print(f'Found {len(valleys)} valleys for the spectrum'
                  f' \"{self.info["name"]}\", returning the one'
                  f' with higher prominence.')

        best_match = np.argmax(properties['prominences'])
        x, y = xs[valleys[best_match]], ys[valleys[best_match]]

        self.info['resonant_wl'] = x
        self.info['resonant_wl_power'] = y
        return x, y

    def plot_spectrum(self, plot_opts=None, subplots=None):
        """
        plota o espectro em um gráfico, com o vale marcado com um x, se tiver
        sido encontrado

        :param plot_opts: Um dicionário com opções para plotar, se nenhum for
            passado, será usado um padrão definido na função.
        :type plot_opts: dict

        :param subplots: Os objetos criados pelo matplotlib para fazer um plot,
            retornados pelo 'plt.subplots()'. Geralmente salvos como (fig, ax).
        :type subplots: tuple

        :return: None
        """
        if not self.quiet:
            print(f'Plotando {self.info["name"]}')

        plot_opts = plot_opts or {
            'xlim': (self.spectrum[0, 0], self.spectrum[-1, 0]),
            'ylim': (-100, 0), 'animated': False, 'interval': 1e-6}

        if subplots is None:
            fig, ax = plt.subplots()
        else:
            fig, ax = subplots

        ax.set_xlim(plot_opts['xlim'])
        ax.set_ylim(plot_opts['ylim'])
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Optical power (dBm)')

        xs = self.spectrum[::, 0]
        ys = self.spectrum[::, 1]

        ax.plot(xs, ys)
        if 'resonant_wl' in self.info.keys() is not None:
            ax.plot(self.info['resonant_wl'], self.info['resonant_wl_power'],
                    'xk')


class MassSpectraData:
    """
    Uma classe usada para processar vários espectros de uma vez.
    """
    def __init__(self, filenames, output_path=None, add_info=None,
                 batch_size=50):
        """
        Inicia o objeto, criando umas variáveis necessárias..

        :param filenames: Um array like com os nomes dos arquivos dos
            espectros a serem abertos.
        :type filenames: array like

        :param output_path: O caminho da pasta onde os arquivos vão ser salvos.
        :type output_path: str

        :param add_info: Um array like contendo os nomes das informações a
            serem adicionadas no arquivo de saída (como temp, strain, etc...)
        :type add_info: array like

        :param batch_size: A quantidade de espectros a serem analisados entre
            cada salvamento de segurança (checkpoint.csv)
        :type batch_size: int
        """
        self.filenames = list(filenames)
        self.out_path = output_path or 'output'
        if not os.path.isdir(self.out_path):
            os.mkdir(self.out_path)

        self.steps = list()
        self.kwargs = list()
        self.add_info = add_info
        self.batch_size = batch_size

        self.columns = ['name', 'resonant_wl', 'resonant_wl_power']
        if add_info is not None:
            for info in add_info:
                self.columns = self.columns.append(info)

        self.df = pd.DataFrame(columns=self.columns)

    def add_step(self, step, kwargs=None):
        """
        Adiciona uma função para ser aplicada à todos os espectros, com os
        argumentos definidos no dicionário kwargs. Se a função tiver
        argumentos obrigatórios, esses devem estar inclusos no kwargs.

        :param step: A função a ser adicionada. Deve aceitar um objeto da
            classe SpectrumData como primeiro argumento ou ser uma função da
            classe.
        :type step: Uma função

        :param kwargs: Um dicionário contendo pares com a key sendo uma str
            com o nome de um argumento da função em step e value o valor desse
            argumento.
        :type kwargs: dict, optional

        :return: None
        """
        self.steps.append(step)
        kwargs = kwargs or dict()
        self.kwargs.append(kwargs)

    def run(self):
        """
        Aplica todas as funções em self.steps a todos os espectros (um por
        vez).  Ao finalizar as funções de um espectro, tenta passar as
        informações definidas no self.columns para o dataframe do objeto.
        Salva checkpoints a cada intervalo de self.batch_size e um último
        no final.

        :return: None
        """
        try:
            index = self.steps.index(SpectrumData.plot_spectrum)
            self.kwargs[index]['subplots'] = plt.subplots()
        except ValueError:
            pass

        for i, filename in enumerate(self.filenames):
            print(f'{5 * "-"}calculando {i + 1}/{len(self.filenames)}{5 * "-"}')
            spectrum = SpectrumData(filename)
            for c, step in enumerate(self.steps):
                ret = step(spectrum, **self.kwargs[c])
                if type(ret) == SpectrumData:
                    spectrum = ret

            self.df = self.df.append(spectrum.info, ignore_index=True)

            if (i % self.batch_size) == 0 and i != 0:
                self.export_csv(os.path.join(self.out_path, 'checkpoint.csv'))

        self.export_csv(os.path.join(self.out_path, 'checkpoint.csv'))

    def export_csv(self, export_name):
        """
        Exporta o dataframe atual como um csv

        :param export_name: O nome do arquivo de saída
        :type export_name: str

        :return: None
        """
        if export_name[-4:] != '.csv':
            export_name += '.csv'

        self.df.to_csv(export_name, index=False, sep=',', decimal='.')
