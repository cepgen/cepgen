import Cards.utils_cfi as cg
from Cards.integrators_cff import miser as integrator
#from Cards.pythia8_cff import pythia8 as hadroniser

process = cg.Parameters(
    name = 'lpair',
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
    ),
    tamingFunctions = cg.Parameters(
        #{ variable = "m_central", expression = "(m_central>80.) ? exp(-(m_central-80)/10) : 1.0" } // example of a complex taming function
        #{ variable = "q2", expression = "exp(-q2)" }
    ),
)

#--- either use the default generation (100k events)
#from Cards.generator_cff import generator

#--- or let the user specify the run conditions
generator = cg.Parameters(
    numEvents = 100000,
    printEvery = 10000,
)
