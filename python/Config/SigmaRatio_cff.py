##
# \file
# \ingroup python
#
# Enumeration of \f$\sigma_R/\sigma_L\f$ computation algorithms

from .containers_cff import Module


class SigmaRatio:
    """R-ratio computation method"""
    E143             = 1
    R1990            = 2
    CLAS             = 3
    SibirtsevBlunden = 4
