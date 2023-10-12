"""@package collinear parton-factorised fluxes definition

A collection of beam-collinear parton fluxes f(x, Q^2)
"""

from .containers_cfi import Module, Parameters


class ProtonFlux:
    """Type of parton (from proton) flux modelling"""
    PhotonElastic = Module('EPAFlux',
        formFactors = Module('StandardDipole')
    )
    PhotonElasticDZ = Module('DreesZeppenfeld')
    PhotonInelastic = Module('EPAFlux',
        formFactors = Module('InelasticNucleon')
    )
    IntegratedPhotonElastic = Module('KTIntegrated',
        ktFlux = Module('BudnevElastic')
    )
    IntegratedPhotonInelastic = Module('KTIntegrated',
        ktFlux = Module('BudnevElastic',
            formFactors = Module('InelasticNucleon')
        )
    )
    def LHAPDF(pdfset='', fromRemnant: bool=False):
        return Module('LHAPDF',
            set = pdfset,
            fromRemnant = fromRemnant)
    LHAPDFLUXlep = LHAPDF(
        pdfset = 'LUXlep-NNPDF31_nlo_as_0118_luxqed',
        fromRemnant = False
    )


class HeavyIonFlux:
    """Type of parton (from heavy ion) flux modelling"""
    PhotonElastic = Module('KTIntegrated',
        ktFlux = Module('BudnevElastic',
            formFactors = Module('HeavyIonDipole')
        )
    )


class ElectronFlux:
    """Type of parton (from electron) flux modelling"""
    PhotonElastic = Module('KTIntegrated',
        ktFlux = Module('BudnevElastic',
            formFactors = Module('PointLikeFermion')
        )
    )
