import Config.Core as cepgen
import Config.ktProcess_cfi as kt
#from Config.Hadronisation.pythia6_cfi import pythia6 as hadroniser
#from Config.Hadronisation.pythia8_cfi import pythia8 as hadroniser
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator
from FormFactors.heavyIonDipole_cfi import heavyIonDipole
from OutputModules.dump_cfi import dump as _dump_output # periodic event printout
from OutputModules.text_cfi import text as _text_output # ASCII histograms


process = kt.process.clone('pptoff',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        heavyIon1 = (208, 82),
        heavyIon2 = (208, 82),
        formFactors = [heavyIonDipole, heavyIonDipole],
        pz = (2562.2, 2562.2),
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (5.,),
        eta = (-2.5, 2.5)
    )
)

# events generation
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
