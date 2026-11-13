import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.generator_cfi import generator as _gen
from Config.PDG_cfi import PDG


# gammaUPC single-beam photon flux (proton beam at 6.5 TeV)
gammaUPC_flux = cepgen.Module('gammaUPC', pdgId = PDG.proton, gamma = 6500. / 0.938272)

process = coll.process.clone('pptoff',
    partonsGenerator = 'coll',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
        method = 0,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        partonFluxes = (gammaUPC_flux, gammaUPC_flux),
        structureFunctions = cepgen.StructureFunctions.suriYennie,
    ),
    outKinematics = coll.process.outKinematics.clone(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

generator = _gen.clone(
    numEvents = 10000
)
