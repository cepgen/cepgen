import Config.Core as cepgen
import Config.ktProcess_cfi as kt
from Config.Integration.vegas_cff import integrator
from Config.PDG_cfi import PDG, registerParticle

#--- auxiliary SUSY particles definition
registerParticle(1000001, 'sd_l', mass=100., charge=1., fermion=False)
registerParticle(1000024, 'chi_1', mass=100., charge=1., fermion=True)
registerParticle(37, 'H_ch', mass=100., charge=1., fermion=True)

process = kt.process.clone('pptosusy',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.ElasticElastic,
        pair = PDG.sd_l,
        #pair = PDG.chi_1,
        #pair = PDG.H_ch,
    ),
    inKinematics = cepgen.Parameters(
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.SuriYennie,
    ),
    outKinematics = kt.process.outKinematics.clone(
        pt = (0.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
        #ptdiff = (0., 2.5),
        #dely = (4., 5.),
    ),
)

#root = cepgen.Module('root_tree')
text = cepgen.Module('text',
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 500, 20)]),
        'pt(7)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 500, 25)]),
        'acop(7,8)': cepgen.Parameters(xrange=(0., 0.25), nbins=15),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 500.), yrange=(0., 500.), log=True)
    }
)
output = cepgen.Sequence(text)

#--- events generation
from Config.generator_cff import generator
generator.numEvents = 50000
