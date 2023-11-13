"""@package collinear parton-factorised process objects definition

A collection of useful objects for the definition of a
general collinear parton momentum-factorised process steering card
"""

from .containers_cfi import Module, Parameters
from .collinearFluxes_cff import *


process = Module('collinearProcess',
    ktFactorised = False,
    outKinematics = Parameters(
        #--- cuts on initial-state partons
        q2 = (0., 10.),
        #--- cuts on individual particles defining the central system
        rapidity = (-6., 6.),
    ),
)
