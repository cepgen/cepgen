import Config.Core as cepgen
from Config.Integration.vegas_cff import integrator
#from Config.Hadronisation.pythia6_cff import pythia6 as hadroniser
#from Config.Hadronisation.pythia8_cff import pythia8 as hadroniser
from Config.PDG_cfi import PDG

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
    #tamingFunctions = cepgen.Parameters(
    #    # example of a complex taming function
    #    cepgen.Parameters(variable = "m_ll", expression = "(m_ll>80.) ? exp(-(m_ll-80)/10) : 1.0"),
    #),
)

#--- example of an output module parameterisation
#output = cepgen.Module('text', variables = ['nev', 'm(4)', 'tgen'])
#output = cepgen.Module('lhef', output = 'test.lhe')
#output = cepgen.Module('hepmc', output = 'test.hepmc')

#--- let the user specify the run conditions
from Config.generator_cff import generator
generator = generator.clone(
    numEvents = 100000,
    printEvery = 10000,
)

