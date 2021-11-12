import numpy as np
import sys

from dark_emulator import darkemu
from lsst.cp.pipe.utils import (funcPolynomial, irlsFit)


# Fit Delta sigma measurement using emulator

emu = darkemu.base_class()

cparam = np.array([0.02225,0.1198,0.6844,3.094,0.9645,-1.])
emu.set_cosmology(cparam)


def funcDeltaSigma (M, rs):
    return emu.get_DeltaSigma_mass(rs, M, 0.0)

file1 = np.genfromtxt(sys.argv[1])

r1, y1, y1err = file1[1:,7], file1[1:,5], file1[1:,11]


print (r1)

mask= (r1>0.5) & (r1<10)

r1=r1[mask]
y1=y1[mask]
y1err=y1err[mask]

polyFit, polyFitErr, chiSq, weights = irlsFit([1e13], r1, y1, funcDeltaSigma) #, weightsY=y1err)

print (polyFit, polyFitErr, chiSq, weights)

