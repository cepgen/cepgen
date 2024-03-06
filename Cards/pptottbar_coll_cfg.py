import Config.Core as cepgen
import Config.collinearProcess_cfi as coll
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator


process = coll.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.top,
        method = 0,  # only on-shell method supported with collinear emission
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = coll.process.outKinematics.clone(
        #eta = (-2.5, 2.5),
        mx = (1.07, 2000.),
        ptdiff = (0., 2000.),
    ),
)

generator.numEvents = 25000
