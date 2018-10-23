"""@package kT-factorised process objects definition

A collection of useful objects for the definition of a
general kT-factorised process steering card
"""

from containers_cfi import Module, Parameters
from math import pi

class ProtonFlux:
    '''Type of parton (from proton) flux modelling'''
    PhotonElastic         = 0
    PhotonInelastic       = 1
    PhotonInelasticBudnev = 11
    GluonKMR              = 20
class HeavyIonFlux:
    '''Type of parton (from heavy ion) flux modelling'''
    PhotonElastic         = 100

process = Module('ktProcess',
    outKinematics = Parameters(
        #--- cuts on initial-state partons
        qt = (0., 50.),
        phiqt = (0., 2.*pi),
        #--- cuts on individual particles defining the central system
        rapidity = (-6., 6.),
        #--- cuts on the pt(outgoing system) (hyper-)plane
        ptdiff = (0., 500.),
        phiptdiff = (0., 2.*pi),
    ),
)
