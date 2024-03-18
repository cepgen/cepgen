##
# \file
# \defgroup kt kt-factorised processes
# \ingroup python kt
#
# A collection of objects for the definition of a general \f$k_{\rm T}\f$-factorised process

from math import pi
from .containers_cff import Module, Parameters
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
