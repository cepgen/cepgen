"""@package kT-factorised process objects definition

A collection of useful objects for the definition of a
general kT-factorised process steering card
"""

from math import pi
from .containers_cfi import Module, Parameters
from .ktFluxes_cff import *


process = Module('ktProcess',
    ktFactorised = True,
    outKinematics = Parameters(
        #--- cuts on initial-state partons
        qt = (0., 50.),
        phiqt = (0., 2.*pi),
        #--- cuts on individual particles defining the central system
        rapidity = (-6., 6.),
        #--- cuts on the pt(outgoing system) (hyper-)plane
        ptdiff = (0., 500.),
        phidiff = (0., 2.*pi),
    ),
)
