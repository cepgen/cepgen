import Cards.Core as cepgen
from Cards.integrators_cff import miser as integrator
#from Cards.pythia8_cff import pythia8 as hadroniser

process = cepgen.Module('pptoll',
    mode = cepgen.ProcessMode.ElasticElastic,
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions(
            'Suri-Yennie'
            #'Fiore'
        ),
    ),
    outKinematics = cepgen.Parameters(
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

#--- events generation
from Cards.generator_cff import generator
generator.numEvents = 1e4

