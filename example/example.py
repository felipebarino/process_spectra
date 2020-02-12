import process_spectra as ps
import numpy as np
from matplotlib import pyplot as plt
from scipy import stats as st
import os 

path = os.path.dirname(os.path.realpath(__file__)) + '/'

sd = ps.SpectraData(path + 'data/')
# o OSA exporta o compriemnto de onda em nm, devemos multiplicar por 1e-9 na importação
sd.importSpectra(wl_multiplier=1e-9)
# dicionário com opções de plotar (é usado pelo pacote)
plot_opts = {'xlim': (1470, 1530), 'ylim': (-70, -50), 'animated': False, 'interval': 1e-6}
sd.plotSpectra(plot_opts)

# filtrar o espectro com savitz-golay
window_size = 45
order = 3
sd.filterSpectra(window_size, order)

# interpola os espectros para os compriemntos de onda abaixo
wl = np.arange(1470e-9, 1530e-9, 0.5e-12)
sd.interpolateSpectra(wl, 'cubic')
#detecta os comprimentos de onda ressonantes
threshold = 5
sd.updateWlRes(threshold)
# seta os valores de mensurando utilizados no experimento
temps = np.array([21.9,23,24.6,29.8,32.3,37.4])
sd.setMeasurements(temps)
# exporta os comprimentos de onda ressonantes e mensurandos para um excel
sd.exportData(path + "export.xls")

# curva de calibração
wl_res = sd.getWlRes()
plt.figure()
plt.plot(wl_res*1e9, temps, 'or', label='Experimental data')
plt.xlabel('$\lambda_{res}$ (nm)')
plt.ylabel('Temperature (°C)')

slope, intercept, r_value, p_value, std_err = st.linregress(wl_res[::,0], temps)

wl = np.linspace(min(wl_res), max(wl_res), 10)

plt.plot(wl*1e9, wl*slope + intercept, '-k', label='Fitted curve')
print('Fitting results:')
print('temp = ' + str(slope) + '*lambda_res + ' + str(intercept))
print('R²: ' + str(r_value**2))
