import Cards.utils_cfi as cg
from Cards.integrators_cff import miser as integrator
from Cards.pythia8_cff import pythia8 as hadroniser

process = cg.Parameters(
    name = "pptoww",
    mode = cg.ProcessMode.InelasticInelastic,
    #'mode': "elastic/elastic",
    inKinematics = cg.Parameters(
        cmEnergy = 13000.,
        #structureFunctions = cg.StructureFunctions('Szczurek-Uleshchenko'),
        #structureFunctions = cg.StructureFunctions('ALLM', '91'),
        #structureFunctions = cg.StructureFunctions('ALLM', '97'),
        structureFunctions = cg.StructureFunctions('LUXlike'),
        #structureFunctions = cg.StructureFunctions('Suri-Yennie'),
        #structureFunctions = cg.StructureFunctions('Fiore'),
    ),
    outKinematics = cg.Parameters(
        pt = (0., -1.),
        energy = (0., -1.),
        rapidity = (-6.0, 6.0),
        mx = (1.07, 1000.),
        #--- extra cuts on the pt(W+) and pt(W-) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between W+ and W-
        #rapiditydiff = (4., 5.),
        cuts = {
            11: dict(
                pt = (20., -1),
                eta = (-2.5, 2.5),
            ),
            13: dict(
                pt = (20., -1),
                eta = (-2.5, 2.5),
            )
        },
    )
)

#--- either use the default generation (100k events)
from Cards.generator_cfi import generator

#--- or let the user specify the run conditions
#generator = cg.Parameters(
#    numEvents = 1000,
#    printEvery = 500,
#)

#integrator['numPoints'] = 1000

#print integrator
#print hadroniser
