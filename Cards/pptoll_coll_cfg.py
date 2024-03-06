import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.generator_cfi import generator as _gen
from Config.PDG_cfi import PDG


process = coll.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        pair = PDG.muon,
        method = 0,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        partonFluxes = (coll.ProtonFlux.PhotonElastic, coll.ProtonFlux.PhotonElastic),
        #partonFluxes = (coll.ProtonFlux.PhotonElasticDZ, coll.ProtonFlux.PhotonElasticDZ),
        #partonFluxes = (coll.ProtonFlux.IntegratedPhotonElastic, coll.ProtonFlux.IntegratedPhotonElastic),
        #partonFluxes = (coll.ProtonFlux.LHAPDF(pdfset='cteq6l1', extrapolatePDF=True),
        #                coll.ProtonFlux.LHAPDF(pdfset='cteq6l1', extrapolatePDF=True)),
        #partonFluxes = (coll.ProtonFlux.LHAPDFLUXlep, coll.ProtonFlux.LHAPDFLUXlep),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
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
