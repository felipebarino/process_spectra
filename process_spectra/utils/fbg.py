"""Esse módulo tem a função para calcular a refletância da fbg"""

import numpy as np


def get_fbg_reflectance(wl_bragg, fwhm, wls, unit='dB'):
    """
    Aproximação do espectro de reflexão de uma FBG

    Referência da eq.:
    @article{peternella2017,
      title={Interrogation of a ring-resonator ultrasound sensor using a fiber Mach-Zehnder interferometer},
      author={Peternella, Fellipe Grillo and Ouyang, Boling and Horsten, Roland and Haverdings, Michael and Kat, Pim and Caro, Jacob},
      journal={Optics express},
      volume={25},
      number={25},
      pages={31622--31639},
      year={2017},
      publisher={Optical Society of America}
    }

    :param wl_bragg: (bragg wavelengths) comprimento de onda de bragg da FBG
    :type wl_bragg: float

    :param fwhm: (full width at half maximum) largura de banda da FBG
    :type fwhm: float

    :param wls: (wavelengths) comprimento de onda ressonante da simulação
    :type wls: np array

    :param unit: unidade da saída:
        * dB (standard)
        * linear

    :return: Retorna um np array com o espectro de refleção da FBG

    """
    r = (1 + ((wls - wl_bragg)/(fwhm/2))**8)**(-1)
    if unit == 'dB':
        r = 10*np.log10(r)
    elif unit != 'scalar':
        raise ValueError('Unidade inválida. As implementadas são "dB" e '
                         '"scalar"')

    return np.array([x for x in zip(wls, r)],
                    dtype=np.float64)

