import Config.Core as cepgen
import Config.ktProcess_cfi as kt
#from Config.Hadronisation.pythia6_cfi import pythia6 as hadroniser
#from Config.Hadronisation.pythia8_cfi import pythia8 as hadroniser
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator
from OutputModules.dump_cfi import dump as _dump_output # periodic event printout
from OutputModules.text_cfi import text as _text_output # ASCII histograms

#--- example of an auxiliary particles definition
#registerParticle(1000001, 'sd_l', mass=100., charge=1., fermion=True) # right now, only fermionic coupling handled

process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = PDG.muon,
        #pair = PDG.up,
        #pair = PDG.sd_l, # whatever was defined above as "new" particle
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        #pdgIds = (PDG.proton, PDG.electron),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.suriYennie,
        #structureFunctions = cepgen.StructureFunctions.fioreBrasse,
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

#--- events generation
generator.numEvents = 10000
text = _text_output.clone(
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 250, 10)]),
        'm(ob1)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
dump = _dump_output.clone(printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
