import Config.Core as cepgen
from Config.integrators_cff import miser as integrator
#from Config.pythia8_cff import pythia8 as hadroniser
from Config.ktProcess_cfi import ktProcess

process = ktProcess.clone('pptoll',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
    ),
    outKinematics = ktProcess.outKinematics.clone(
        pair = 13,
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 10000
