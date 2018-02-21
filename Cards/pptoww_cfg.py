import Cards.Core as cepgen
from Cards.integrators_cff import miser as integrator
#from Cards.pythia8_cff import pythia8 as hadroniser

process = cepgen.Module('pptoww',
    mode = cepgen.ProcessMode.ElasticElastic,
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
        #structureFunctions = cepgen.StructureFunctions.ALLM91,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        pt = (0., -1.),
        energy = (0., -1.),
        rapidity = (-6.0, 6.0),
        mx = (1.07, 1000.),
        #--- extra cuts on the pt(W+) and pt(W-) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between W+ and W-
        #rapiditydiff = (4., 5.),
        cuts = {
            11: cepgen.Parameters(
                pt = (20., -1),
                eta = (-2.5, 2.5),
            ),
            13: cepgen.Parameters(
                pt = (20., -1),
                eta = (-2.5, 2.5),
            )
        },
    )
)

#--- either use the default generation (100k events)
from Cards.generator_cff import generator

#--- or let the user specify the run conditions
#generator = cepgen.Parameters(
#    numEvents = 1000,
#    printEvery = 500,
#)

integrator.numPoints = 10000
generator.numEvents = 100000

#print(process)
#print(integrator)
#print(hadroniser)
