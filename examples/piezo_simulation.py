"""
Esse script simula as potências ópticas refletidas por fbgs, percebidas por
fotodetectores (um em cada fbg). Essas fbgs são acopladas a componentes
piezoelétricos, o que as faz vibrar e, consequentemente, modular seus vales
ressonantes.

A leitura do sensor é simulada ao longo do tempo para obter os
componentes harmônicos da frequência da tensão aplicada ao piezoelétrico.
"""
import os
import numpy as np
from process_spectra import MassSpectraData, utils, funcs
from process_spectra.funcs.piezo_fbg import simulate_piezo_fbg


files = os.listdir('data/spectra')

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
          (1.47411538420759E-06, 1.6397432121049e-06)}
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
          'fbg_wl_bragg': 1537.5e-9,
          'fbg_fwhm': 1e-9,
          'fbg_label': 'fbg1'}
spectra.add_step(step, kwargs)

kwargs = {'frequency': 60,
          'amplitude': 127,
          'fbg_wl_bragg': 1547.5e-9,
          'fbg_fwhm': 1e-9,
          'fbg_label': 'fbg2'}
spectra.add_step(step, kwargs)

spectra.run()
spectra.export_csv('csvs/piezo.csv')
