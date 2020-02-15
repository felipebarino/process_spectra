# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 17:58:09 2020

@author: Felipe Barino
"""

import sys
import codecs
import os, os.path
import numpy as np
from scipy import signal as sg

class SpectraData:
    spectra = list()
    wl_res = list()
    wl_res_att = list()
    files = list()
    measurements = np.array(list())
    
    wl = list()
    optical_powers = list()
    
    """ Construtor:
    Ao iniciar o objeto passar o caminho da pasta que contém os espectros
    
    OBS.: a) Tipo: str; b) Usar '/' no caminho
    Exemplo:
        C:/Users/Usuario/Documents/My_Experiment/Spectra/
    """
    def __init__(self, path):
        if path[-1] != '/':
            path = path + '/'
        self.path = path
        self.updateFiles()
    
    """ updateFiles
    Atualiza o nome dos arquivos da pasta onde estão os espectros
    """
    def updateFiles(self):
        self.files = list()
        for name in os.listdir(self.path):
            if os.path.isfile(os.path.join(self.path, name)):
                self.files.append(self.path + name)
    
    """ getFiles
    Retorna o caminho dos arquivos que contém os espectros
    """
    def getFiles(self):
        return self.files
    
    """ fixOptisystem
    Arruma o formato do texto que foi salvo os espectros pelo Optisystem
    """
    def fixOptisystem(self, file):
        f = codecs.open(file, encoding='utf-8')
        contents = f.read()
        newcontents = contents.replace('E', 'e')
        newcontents = newcontents.replace(',', '.')
        f = open(file, 'w')
        f.write(newcontents)
        f.close()
    
    """ fixOptisystemFormats
    Aplica fixOptisystem para todos os arquivos contidos em files
    """
    def fixOptisystemFormats(self):
        for file in self.files:
            self.fixOptisystem(file)
    
    """ setMeasurements
    Seta quais foram os valores que geraram os respectivos espectros da pasta
    
    OBS.: a) Tipo: np.array() b) Usar a mesma ordem da variável files  
    """
    def setMeasurements(self, measurements):
        self.measurements = measurements
    
    """ updateOptisystemSpectra
    Atualiza a variável spectra com os espectros que têm o caminho descrito 
    por files
    
    separator: char que separa os valores, se os arquivos são no formato:
    1.2e-6, -80
    1.3e-6, -81
    ...
    1.8e-6, -85
    então deve utilizar: separator=','
    
    OBS.: caso não utilize o parâmetro, o padrão ';' será utilizado
    """
    def importSpectra(self, wl_multiplier=1, separator=';'):
        self.spectra = list()
        for file in self.files:
            try:
                spectrum = np.loadtxt(fname=file, delimiter=separator)
            except:
                self.fixOptisystem(file)
                spectrum = np.loadtxt(fname=file, delimiter=separator)
            if spectrum.shape[0] < spectrum.shape[1]:
                spectrum = spectrum.transpose()
            spectrum[::,0] = wl_multiplier*spectrum[::,0]
            self.spectra.append(spectrum)
    
    """ maskSpectra
    Corta o espectro com base em um intervalo de comprimento de onda
    
    wl_limits: é uma lista de dois valores com os comprimentos de onda máximo
    e mínimo que se deseja obter no espectro final
    """
    def maskSpectra(self, wl_limits):
        aux = list()
        for spectrum in self.spectra:
            wl = spectrum[::, 0]
            a = wl >= min(wl_limits)
            b =  wl <= max(wl_limits)
            mask = [all(tup) for tup in zip(a, b)]
            aux.append(spectrum[mask, ::])
        self.spectra = aux
    
    """ filterSpectra
    Filtra os espectros utilizando filtro de Savitz-Golay
    https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter
    
    size: tamanho do filtro
            window_length : (int) The length of the filter window 
                            (i.e. the number of coefficients). 
                            window_length must be a positive odd
    order: ordem do polinômio
            polyorder :     (int) The order of the polynomial used 
                            to fit the samples. polyorder must be less 
                            than window_length.
    
    keep_old: se mantém ou não os espectros não filtrados
                se mantém: apenas retorna o filtrado
                se não: só atualiza os espectros que estão no objeto da classe
            Por padrão, os espectros antigos não são mantidos
    """
    def filterSpectra(self, size, order, keep_old=False):       
        aux = list()
        for spectrum in self.spectra:
            spectrum[::, 1] =  sg.savgol_filter(spectrum[::, 1], size, order)
            aux.append(spectrum)
        if keep_old:
            return aux
        else:
            self.spectra = aux
    
    """ interpolateSpectra
    Interpola os espectros com base em um vetor pré-fixado de comprimento de onda
    
    wl: vetor de comprimentos de onda com a resolução desejada
    kind: tipo de função interpolante (‘linear’, ‘nearest’, ‘zero’, ‘slinear’, 
                                       ‘quadratic’, ‘cubic’, ‘previous’, ‘next’)
    """
    def interpolateSpectra(self, wl, kind):
        self.optical_powers = list()
        from scipy import interpolate
        self.wl = wl
        for spectrum in self.spectra:
            interpolant = interpolate.interp1d(spectrum[::,0], spectrum[::,1], kind=kind)
            self.optical_powers.append(interpolant(wl))
    
    """ detectWlRes
    Detecta o comprimento de onda ressonante
    
    wl: vetor dos compriemntos de onda ressonante do espectro
    power: vetor com potência optica do espectro
    threshold: limiar para detecção do vale
    """
    def detectWlRes(self, wl, power, threshold):
        peaks = sg.find_peaks(-power, prominence=threshold)
        peaks = peaks[0]
        if len(peaks) < 1:
            print('Resonant wavelength detection error!')
            return -1
        else:
            return wl[peaks], power[peaks]
    
    """ updateWlRes
    Encontra o comprimento de onda resonante de cada espectro
    
    threshold: limiar para detecção do vale
    use_interp_data: se usa os dados interpolados ou não, por padrão: não usa
    """
    def updateWlRes(self, threshold, use_interp_data=True):
        self.wl_res = list()
        if use_interp_data and len(self.optical_powers) < 1:
            print('Spectra has not been interpolated\n' + 
                  'Use interpolateSpectra before using use_interp_data as True')
            return -1

        if use_interp_data:
            for optical_power in self.optical_powers:
                wl, power = self.detectWlRes(self.wl, optical_power, threshold)
                self.wl_res.append(wl)
                self.wl_res_att.append(power)
        else:
            for spectrum in self.spectra:
                wl, power = self.detectWlRes(spectrum[::,0], spectrum[::,1], threshold)
                self.wl_res.append(wl)   
                self.wl_res_att.append(power)
    
    """ getWlRes
    Retorna os comprimentos de onda ressonante
    """
    def getWlRes(self):
        return np.array(self.wl_res)
    
    """ plotSpectra
    Mostra os espectros
    
    plot_opts: dicionário com opções do plot
    """
    def plotSpectra(self, plot_opts):
        from matplotlib import pyplot as plt
        
        fig, ax = plt.subplots()
        ax.set_xlim(plot_opts['xlim'])
        ax.set_ylim(plot_opts['ylim'])
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('Optical power (dBm)')
        i = 0
        for spectrum in self.spectra:
            if len(self.wl_res) > 0:
                ax.plot(self.wl_res[i]*1e9, self.wl_res_att[i], 'xk')  
            ax.plot(spectrum[::,0]*1e9, spectrum[::,1])
            i+=1
    
    """ getSpectra
    Retorna uma lista com os espectros
    
    interpolated: se retorna os espectros interpolados junto com o wl
    """
    def getSpectra(self, interpolated=False):
        if interpolated:
            return self.wl, self.optical_powers
        else:
            return self.spectra
    
    """ exportData
    Exporta o comprimento de onda ressonante e os mensurandos
    
    path: caminho de saída para o arquivo ser salvo, incluindo nome do arquivo 
    com extensão do mesmo (csv ou xls)
    """
    def exportData(self, path):
        import pandas as pd
        data = [self.measurements, np.array(self.wl_res)[::,0]]
        cols = ['Measure', 'lambda_res']
        df = pd.DataFrame(data, index=cols).transpose()
        
        _, extension = path.split(sep='.')
        if extension == 'csv':
            df.to_csv(path, index=False)
        elif extension == 'xls':
            df.to_excel(path, index=False)
        else:
            print('Sorry, ' + extension + ' is an unsuported format.')
