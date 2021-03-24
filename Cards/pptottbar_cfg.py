import Config.Core as cepgen
import Config.ktProcess_cfi as kt
#from Config.Hadronisation.pythia6_cff import pythia6
#from Config.Hadronisation.pythia8_cff import pythia8
from Config.PDG_cfi import PDG

#--- redefinition of top to modify its bare mass
#registerParticle(6, 'top', mass=174., charge=2./3., fermion=True)

#--- process definition
process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.top,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = kt.process.outKinematics.clone(
        #eta = (-2.5, 2.5),
        qt = (0., 2000.),
        mx = (1.07, 2000.),
        #--- extra cuts on the p1t(t) and p2t(t) plane
        ptdiff = (0., 2000.),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)

#--- events modification sequence
eventSequence = cepgen.Sequence(
    #pythia6,
)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 25000
