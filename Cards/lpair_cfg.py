import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator
from Config.PDG_cfi import PDG

#--- example of a hadronisation algorithm steering
#from Config.Hadronisation.pythia6_cff import pythia6 as hadroniser
#from Config.Hadronisation.pythia8_cff import pythia8 as hadroniser

process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticElastic,
        pair = PDG.muon,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
        #structureFunctions = cepgen.StructureFunctions.FioreBrasse,
        #structureFunctions = cepgen.StructureFunctions.LUXlike,
    ),
    outKinematics = cepgen.Parameters(
        pt = (25.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
    #--- example of a complex taming function definition
    #tamingFunctions = cepgen.Parameters(
    #    cepgen.Parameters(variable = "m_ll", expression = "(m_ll>80.) ? exp(-(m_ll-80)/10) : 1.0"),
    #),
)

#--- example of an output module parameterisation
#output = cepgen.Module('text', variables = ['nev', 'm(4)', 'tgen'], histVariables={'m(4)': cepgen.Parameters(low=0., high=250.)})
#output = cepgen.Module('lhef', filename='test.lhe')
#output = cepgen.Module('hepmc', filename='test.hepmc')

#--- let the user specify the events generation parameters
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 100000,
    printEvery = 10000,
)

