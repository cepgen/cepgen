"""@package kT-dependent fluxes definition

A collection of beam parton fluxes f(x, kt^2, Q^2)
"""

from .containers_cfi import Module, Parameters


class ProtonFlux:
    """Type of parton (from proton) flux modelling"""
    PhotonElastic = Module('Elastic')
    PhotonInelastic = Module('Inelastic')
    PhotonElasticBudnev = Module('BudnevElastic')
    PhotonInelasticBudnev = Module('BudnevInelastic')
    GluonKMR = Module('KMR')


class HeavyIonFlux:
    """Type of parton (from heavy ion) flux modelling"""
    PhotonElastic = Module('ElasticHeavyIon')


class ElectronFlux:
    """Type of parton (from electron) flux modelling"""
    PhotonElasticBudnev = Module('BudnevElasticLepton',
        pdgId = 11)
