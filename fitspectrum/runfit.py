from fitspectrum import fitspectrum
from spectra import *
import numpy as np

const = get_spectrum_constants()
test1 = freefree(const, 100.0, 1.0, 4000.0, 1.0)
test2 = freefree(const, 101.0, 1.0, 4000.0, 1.0)
alpha = np.log(test1/test2) / np.log((100.0) / (101.0))
print alpha

# fitspectrum('ngc253_spectrum.txt')