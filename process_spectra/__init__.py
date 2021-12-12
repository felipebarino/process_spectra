# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:58:09 2020

@author: Felipe Barino, Guilherme Sampaio

Nesse pacote está a biblioteca criada com a intenção de analisar os espectros.

Esse módulo é o básico para simulações e extração de informações a partir de
espectros ópticos. É feito para trabalhar com grandes conjuntos de espectros
com conjuntos de funções específicas para tais simulações e extrações de
informação.
"""

import pandas as pd
from process_spectra.funcs import load_from_optisystem


class MassSpectraData:
    """
    Uma classe usada para processar vários espectros de uma vez.
    """
    def __init__(self, filenames, out_filename=None, load_function=None,
                 quiet=False):
        """
        Inicia o objeto, criando umas variáveis necessárias

        :param filenames: Um array like com os nomes dos arquivos dos
            espectros a serem abertos
        :type filenames: list

        :param out_filename: O arquivo de saída
        :type out_filename: str
        """
        self.filenames = list(filenames)
        self.load_spectrum = load_function or load_from_optisystem
        self.out_filename = out_filename
        self.quiet = quiet

        self.steps = list()
        self.kwargs = list()

        self.df = pd.DataFrame(columns=['name', ])

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

    def run(self, quiet=False, **load_kwargs):
        """
        Aplica todas as funções em self.steps a todos os espectros (um por
        vez).  Ao finalizar as funções de um espectro, tenta passar as
        informações definidas no self.columns para o dataframe do objeto.
        Salva checkpoints a cada intervalo de self.batch_size e um último
        no final.

        :return: None
        """

        for i, filename in enumerate(self.filenames):
            if not quiet:
                padding = 5 * "-"
                print(f'\n{padding}calculando {i + 1}/'
                      f'{len(self.filenames)}{padding}')

            spectrum, info = self.load_spectrum(self.filenames[i], **load_kwargs)

            if spectrum is not None:
                for c, step in enumerate(self.steps):
                    spectrum, _info = step(spectrum, info, **self.kwargs[c])
                    info = {**info, **_info}

            self.df = self.df.append(info, ignore_index=True)

        if self.out_filename:
            self.export_csv(self.out_filename)

    def export_csv(self, export_name):
        """
        Exporta o dataframe atual como um csv

        :param export_name: O nome do arquivo de saída
        :type export_name: str

        :return: None
        """
        if export_name[-4:] != '.csv':
            export_name += '.csv'

        self.df.sort_values('name').\
            to_csv(export_name, index=False, sep=',', decimal='.')

