import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator as _gen


process = kt.process.clone('pptozz',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        polarisationStates = 0,  # full
        model = cepgen.Module("dim6EFTzz",
            a0z = 1.e-6,  # in GeV^-2
        ),
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        cmEnergy = 13.e3,
        structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = kt.process.outKinematics.clone(
        mx = (1.07, 1000.),
        qt = (0., 1000.),
        #ptdiff = (0., 2000.),
        invmass = (0.,),
        ptsum = (0.,),
        #xi = (0.02, 0.15),
    )
)

generator = _gen.clone(  # generation parameters
    numEvents = 50000,
    printEvery = 5000,
)
text = cepgen.Module('text',  # histogramming/ASCII output capability
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xrange=(50., 1000.), nbins=19),
    }
)
output = cepgen.Sequence(text)
