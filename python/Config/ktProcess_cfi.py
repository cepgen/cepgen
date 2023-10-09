"""@package kT-factorised process objects definition

A collection of useful objects for the definition of a
general kT-factorised process steering card
"""

from math import pi

from .containers_cfi import Module, Parameters

class ProtonFlux:
    """Type of parton (from proton) flux modelling"""
    PhotonElastic         = Module('Elastic')
    PhotonInelastic       = Module('Inelastic')
    PhotonElasticBudnev   = Module('BudnevElastic')
    PhotonInelasticBudnev = Module('BudnevInelastic')
    GluonKMR              = Module('KMR')


class HeavyIonFlux:
    """Type of parton (from heavy ion) flux modelling"""
    PhotonElastic         = Module('ElasticHeavyIon')


class ElectronFlux:
    """Type of parton (from electron) flux modelling"""
    PhotonElasticBudnev   = Module('BudnevElasticLepton', pdgId = 11)


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
