""" package implementing useful rayoptics graphics using matplotlib

		The :mod:`~.mpl` subpackage provides useful basic optical graphics
		using the matplotlib plotting package. Particular features include:

				- 2D lens layout, :mod:`~.interactivelayout`
				- |ybar| and |nubar| paraxial ray diagrams, :mod:`~.interactivediagram`
				- ray aberration and wavefront pupil/field plots,
					:mod:`~.analysisfigure`, :mod:`~.axisarrayfigure` and
					:mod:`~.analysisplots`
				- base class to manage light and dark UI styles, :mod:`~.styledfigure`
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


def plot_z_enp(flds, obj_angs, img_hts, z_enps, plot_vs='object angle'):
    """ Plot z_enp vs `object angle` or `image height` """
    if plot_vs == 'object angle':
        vs = obj_angs
    elif plot_vs == 'image height':
        vs = img_hts
    tck = interpolate.splrep(vs, z_enps)
    polyline = np.linspace(vs[0], vs[-1], 50)
    z_enps_interps = interpolate.splev(polyline, tck, der=0)
    
    fig, ax = plt.subplots()
    ax.scatter(z_enps, vs)
    ax.plot(z_enps_interps, polyline)
    ax.set_title("Pupil Aberration")
    ax.set_ylabel(plot_vs)
    ax.set_xlabel("z_enp")
    
    return fig
