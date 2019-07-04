import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.Integration.vegas_cff import integrator
#from Config.Hadronisation.pythia6_cff import pythia6 as hadroniser
#from Config.Hadronisation.pythia8_cff import pythia8 as hadroniser
from Config.PDG_cfi import PDG, registerParticle

#--- auxiliary particles definition
#registerParticle(1000001, 'sd_l', mass=100., charge=1., fermion=True)

process = kt.process.clone('pptoll',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = PDG.muon,
        #pair = PDG.up,
        #pair = PDG.sd_l, # whatever was defined above as "new" particle
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #pdgIds = (PDG.proton, PDG.electron),
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #--- extra cuts on the p1t(l) and p2t(l) plane
        #ptdiff = (0., 2.5),
        #--- distance in rapidity between l^+ and l^-
        #dely = (4., 5.),
    ),
)
#print PDG

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 10000
