import Cards.Core as cepgen
from Cards.integrators_cff import vegas as integrator
#from Cards.pythia8_cff import pythia8 as hadroniser

process = cepgen.Module('lpair',
    mode = cepgen.ProcessMode.ElasticElastic,
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        pair = 13,
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
    #tamingFunctions = cepgen.Parameters(
    #    # example of a complex taming function
    #    cepgen.Parameters(
    #        variable = "m_central",
    #        expression = "(m_central>80.) ? exp(-(m_central-80)/10) : 1.0",
    #    ),
    #),
)

#--- either use the default generation (100k events)
#from Cards.generator_cff import generator

#--- or let the user specify the run conditions
generator = cepgen.Parameters(
    numEvents = 100000,
    printEvery = 10000,
)

print integrator.numFunctionCalls
