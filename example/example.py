path = 'example/'

sd = SpectraData(path + 'data/')
# o OSA exporta o compriemnto de onda em nm, devemos multiplicar por 1e-9 na importação
sd.importSpectra(wl_multiplier=1e-9)
# dicionário com opções de plotar (é usado pelo pacote)
plot_opts = {'xlim': (1470, 1530), 'ylim': (-80, -40), 'animated': False, 'interval': 1e-6}
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
sd.setMeasurements(np.array([21.9,23,24.6,29.8,32.3,37.4]))
# exporta os comprimentos de onda ressonantes e mensurandos para um excel
sd.exportData(path + "export.xls")