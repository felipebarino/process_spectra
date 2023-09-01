import numpy as np
import pprint
from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join
from process_spectra.funcs import get_approximate_valley, filter_spectrum, mask_spectrum
from process_spectra.utils import lorentz

spectra_path = '../examples/data/spectra/'

spectra_files = [spectra_path + f for f in listdir(spectra_path) if isfile(join(spectra_path, f))]

lims = [1480, 1600]

plt.figure()
for spectrum_file in spectra_files:
    print(f'\n{spectrum_file}: ')
    spectrum = np.loadtxt(spectrum_file, delimiter=';')

    spectrum, info = filter_spectrum(spectrum, {}, 15, 2, quiet=True)

    spectrum, info = mask_spectrum(spectrum, {}, lims, quiet=True)

    spectrum, info = get_approximate_valley(spectrum, info, approx_func=lorentz,
                                            prominence=2, resolution_proximity=2,
                                            p0=None, dwl=2)

    pprint.pprint(info)

    plt.plot(spectrum[:, 0], spectrum[:, 1])

    for i in range(info['valley_count']):
        if info['valley_count'] > 1:
            x, y = info[f'resonant_wl_{i}'], info[f'resonant_wl_power_{i}']
        else:
            x, y = info[f'resonant_wl'], info[f'resonant_wl_power']
        plt.plot(x, y, 'ok')
plt.show()

