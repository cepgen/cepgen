import Cards.utils_cfi as cg
from Cards.integrators_cff import miser as integrator
#from Cards.pythia8_cff import pythia8 as hadroniser

process = cg.Parameters(
    name = 'pptoll',
    mode = cg.ProcessMode.ElasticElastic,
    inKinematics = cg.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cg.StructureFunctions('Suri-Yennie'),
        #structureFunctions = cg.StructureFunctions('Fiore'),
    ),
    outKinematics = cg.Parameters(
        pair = 13,
        pt = (25.0, -1.),
        energy = (0., -1.),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- either use the default generation (100k events)
#from Cards.generator_cff import generator

#--- or let the user specify the run conditions
generator = cg.Parameters(
    numEvents = 100000,
    printEvery = 10000,
)
