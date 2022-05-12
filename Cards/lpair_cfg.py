import Config.Core as cepgen
from Config.PDG_cfi import PDG
#--- enable timing framework
#from Config.Timer_cfi import timer

process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #--- other structure functions parameterisation are also available
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.LUXlike,
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
    #        variable = "m(4)",
    #        expression = "(m(4)>80.) ? exp(-(m(4)-80)/10) : 1.0"
    #    )
    #],
)

#--- let the user specify the events generation parameters
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 100000,
    printEvery = 10000,
)

#--- example of an events modification procedure
#from Config.Hadronisation.pythia8_cff import pythia8
#eventSequence = cepgen.Sequence(pythia8)

#--- example of an output module(s) procedure
#... dump everything into a flat ROOT tree (if CepGenRoot was built and loaded)
#from Config.OutputModule.ROOTTree_cfi import rootTree
#... accumulate and dump everything into an ASCII text output
text = cepgen.Module('text',
    #variables = ['nev', 'm(4)', 'tgen'],
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 250, 10)]),
        'm(ob1)': cepgen.Parameters(xrange=(0., 250.), nbins=10, log=True),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 250.), yrange=(0., 250.), log=True)
    }
)
#... or write onto one of various event formats handled
#lhef = cepgen.Module('lhef', filename='test.lhe')
#hepmc = cepgen.Module('hepmc', filename='test.hepmc')
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(
    #rootTree,
    text,
    #lhef,
    #hepmc,
    dump,
)
