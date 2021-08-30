"""
Esse script extrai o comprimento de onda ressonante e a potência do vale.
Também pega a potência máxima para contexto e plota os gráficos juntos.
"""

import os
from process_spectra import MassSpectraData
from process_spectra import funcs
import timeit
import matplotlib.pyplot as plt


files = os.listdir('data/spectra')
files_complete = [os.path.join('data/spectra', x) for x in files]

spectra = MassSpectraData(files_complete, 'csvs/resonant_wavelengths.csv')

spectra.add_step(funcs.fill_name_zeros)

step = funcs.filter_spectrum
kwargs = {
    'window_length': 45,
    'polyorder': 3
}
spectra.add_step(step, kwargs)

step = funcs.interpolate_spectrum
kwargs = {'wl_step': 0.5e-12, 'wl_limits': (1.49e-6, 1.62e-6)}
spectra.add_step(step, kwargs)

step = funcs.find_valley
kwargs = {'prominence': 5, 'ignore_errors': True}
spectra.add_step(step, kwargs)

spectra.add_step(funcs.get_max_power)


step = funcs.plot_spectrum
fig, ax = plt.subplots()
kwargs = {'subplots': (fig, ax)}
spectra.add_step(step, kwargs)


start = timeit.default_timer()
spectra.run(ignore_errors=True)
stop = timeit.default_timer()

fig.savefig('plots/resonant_wavelengths.png', transparent=False)

print('')
print('-'*50)
print(f'Time: {stop - start}s')
