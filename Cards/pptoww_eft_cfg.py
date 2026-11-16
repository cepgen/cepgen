import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator as _gen


process = kt.process.clone('pptoww',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        method = 2,  # on-shell (0) or off-shell (2, helicity amplitudes) formula
        polarisationStates = 0,  # full
        #model = cepgen.Module("W", s1 = 1.e-2),
        model = cepgen.Module("dim6EFTww",
            #a0w = 1.e-5,
            a0w = 1.e-6,
            includeSM = False,
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
        #--- extra cuts on the pt(W+) and pt(W-) plane
        ptdiff = (0., 2000.),
        #--- extra cuts on the W+W- system
        invmass = (0.,),
        ptsum = (0.,),
        #--- cuts on single particles' level
        #minFinalState = [11, 13],
        cuts = {
            # cuts on the single W level
            #24: cepgen.Parameters(pt = (0.,)), # no pt cut on Ws
            # cuts on the W decay products
            # (mimicking LHC-like experimental cuts)
            #11: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5)),
            #13: cepgen.Parameters(pt = (20.,), eta = (-2.5, 2.5))
        },
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
