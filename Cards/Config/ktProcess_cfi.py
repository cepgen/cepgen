import Config.Core as cepgen
from math import pi

class ProtonFlux:
    PhotonElastic         = 0
    PhotonInelastic       = 1
    PhotonInelasticBudnev = 11

process = cepgen.Module('ktProcess',
    outKinematics = cepgen.Parameters(
        qt = (0., 50.),
        phiqt = (0., 2.*pi),
        #--- cuts on individual particles defining the central system
        rapidity = (-6., 6.),
        #--- cuts on the pt(outgoing system) (hyper-)plane
        ptdiff = (0., 500.),
        phiptdiff = (0., 2.*pi),
    ),
)
