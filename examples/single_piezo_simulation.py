"""
Esse script simula a potência óptica refletida por uma fbg, percebida por
um fotodetector. Essa fbg é acopladas a um componente piezoelétrico, o que a
faz vibrar e, consequentemente, modular seu vale ressonante.

A leitura do sensor é simulada ao longo do tempo para obter os
componentes harmônicos da frequência da tensão aplicada ao piezoelétrico.

Esse script está com os valores usados para um estudo específico, e não deve
funcionar com os espectros da pasta que está aqui. Pra funcionar basta ajustar
os valores.
"""
import os
import numpy as np
from process_spectra import MassSpectraData, utils, funcs
from process_spectra.funcs.piezo_fbg import simulate_piezo_fbg


files = os.listdir('spectra')

spectra = MassSpectraData([os.path.join('data/spectra', x) for x in files])

spectra.add_step(funcs.set_name_numbers)

step = funcs.filter_spectrum
kwargs = {
    'window_length': 45,
    'polyorder': 3
}
spectra.add_step(step, kwargs)

step = funcs.interpolate_spectrum
kwargs = {'wl_step': 0.5e-12, 'wl_limits':
          (1.500000000000000039e-06, 1.599999999999999927e-06)}
spectra.add_step(step, kwargs)

led = np.loadtxt('data/S5FC1550S-A2.csv', delimiter=', ')
led[:, 0] = led[:, 0] * 1e-9
led = utils.interpolate(original=led, **kwargs)

kwargs = {'other': led}
spectra.add_step(funcs.simulate_gain, kwargs)

step = funcs.find_valley
kwargs = {'prominence': 5}
spectra.add_step(step, kwargs)

step = simulate_piezo_fbg
kwargs = {'frequency': 60,
          'amplitude': 127,
          'fbg_wl_bragg': 1.5424e-6,
          'fbg_fwhm': 1e-9,
          'fbg_label': 'fbg'}
spectra.add_step(step, kwargs)


spectra.run()
spectra.export_csv('final')

