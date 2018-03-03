import Cards.Core as cepgen
from Cards.integrators_cff import miser as integrator
from Cards.ktProcess_cfi import ktProcess
'''from Cards.pythia8_cff import pythia8 as hadroniser

hadroniser.pythiaProcessConfiguration = (
    # process-specific
    '13:onMode = off', # disable muon decays
    '24:onMode = off', # disable all W decays, but...
    '24:onIfAny = 11 13' # enable e-nue + mu-numu final states
)'''

process = ktProcess.clone('pptoww',
    mode = cepgen.ProcessMode.InelasticElastic,
    inKinematics = cepgen.Parameters(
        cmEnergy = 13.e3,
        #structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.SzczurekUleshchenko,
        #structureFunctions = cepgen.StructureFunctions.ALLM91,
        #structureFunctions = cepgen.StructureFunctions.ALLM97,
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = ktProcess.outKinematics.clone(
        mx = (1.07, 1000.),
        qt = (0., 1000.),
        #--- extra cuts on the pt(W+) and pt(W-) plane
        ptdiff = (0., 2000.),
        #--- distance in rapidity between W+ and W-
        #rapiditydiff = (4., 5.),
        cuts = {
            # cuts on the single W level
            24: cepgen.Parameters(
                pt = (0.,),
            ),
            # cuts on the W decay products
            11: cepgen.Parameters(
                pt = (20.,),
                eta = (-2.5, 2.5),
            ),
            13: cepgen.Parameters(
                pt = (20.,),
                eta = (-2.5, 2.5),
            )
        },
    )
)

#integrator.numPoints = 10000

#--- import the default generation parameters
from Cards.generator_cff import generator
generator.numEvents = 10000
#generator.printEvery = 1

#print(process)
#print(integrator)
#print(hadroniser)
