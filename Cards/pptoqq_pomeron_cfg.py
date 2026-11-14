import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator

#--- photon-pomeron (and pomeron-pomeron) production of a quark pair, with pomeron-initiated legs modelled as a
#--- C=+1 isoscalar photon (Donnachie-Landshoff)
#--- switch to (kt.ProtonFlux.PomeronElastic, kt.ProtonFlux.PomeronElastic) for pomeron-pomeron exchange
process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.charm,
        method = 1,  # only the off-shell matrix element is defined for pomeron-initiated modes
        pomeronBeta2 = 3.5,  # GeV^-2, to be kept consistent with the pomeron flux 'beta2' parameter
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        partonFluxes = (kt.ProtonFlux.PhotonElasticBudnev, kt.ProtonFlux.PomeronElastic),
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (3., 500.),
        qt = (0., 50.),
        eta = (-2.5, 2.5),
    ),
)

generator.numEvents = 10000
