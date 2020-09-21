"""
Um exemplo mostrando uma utilização do process_spectra.
Nesse exemplo, um grupo de espectros na pasta local '/spectrums' é lido,
filtrado, interpolado e tem seus valores do vale ressonante encontrados. Ao
final, são plotados em um único gráfico e os comprimentos de onda ressonante
são salvos em um arquivo 'final.csv'.
"""

import os
import sys
from matplotlib import pyplot as plt

sys.path.insert(0, os.path.abspath('..'))
from process_spectra.process_spectra import SpectrumData, MassSpectraData


files = os.listdir('spectrums')

spectra = MassSpectraData([os.path.join('spectrums', x) for x in files])

step = SpectrumData.filter_spectrum
kwargs = {
    'window_length': 45,
    'polyorder': 3
}
spectra.add_step(step, kwargs)

step = SpectrumData.interpolate_spectrum
kwargs = {'wl_step': 0.5e-12, }
spectra.add_step(step, kwargs)

step = SpectrumData.find_valley
kwargs = {'prominence': 5}
spectra.add_step(step, kwargs)

step = SpectrumData.plot_spectrum
spectra.add_step(step)

spectra.run()
spectra.export_csv('final')

plt.waitforbuttonpress()
