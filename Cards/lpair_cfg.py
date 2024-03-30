import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator as _gen
#from Config.timer_cfi import timer # enable timing framework
from OutputModules.dump_cfi import dump as _dump_output # periodic event printout
from OutputModules.text_cfi import text as _text_output # ASCII histograms
#from OutputModules.rootTree_cfi import rootTree # dump everything into a flat ROOT tree


process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.suriYennie,
        #--- other structure functions parameterisation are also available
        #structureFunctions = cepgen.StructureFunctions.fioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = cepgen.Parameters(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
    #--- example of a taming function definition
    #tamingFunctions = [
    #    cepgen.Parameters( # central dilepton mass reweighting
    #        variables = ["m(4)"],
    #        expression = "(m(4)>80.) ? exp(-(m(4)-80)/10) : 1.0"
    #    )
    #],
)

#--- let the user specify the events generation parameters
generator = _gen.clone(
    numEvents = 100000,
    printEvery = 10000,
)

#--- example of an event modification procedure
#from Config.Hadronisation.pythia8_cfi import pythia8
#eventSequence = cepgen.Sequence(pythia8)

#--- example of an output module(s) procedure
text = _text_output.clone(
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 250, 10)]),
        'm(ob1)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
#lhef = cepgen.Module('lhef', filename='test.lhe')
#hepmc = cepgen.Module('hepmc', filename='test.hepmc')
dump = _dump_output.clone(printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
