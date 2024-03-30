import Config.Core as cepgen
from Config.PDG_cfi import PDG
from Config.generator_cfi import generator as _gen
#from Config.timer_cfi import timer  # enable timing framework


process = cepgen.Module('lpair',
    processParameters = cepgen.Parameters(
        mode = cepgen.ProcessMode.InelasticInelastic,
        pair = PDG.top,
    ),
    inKinematics = cepgen.Parameters(
        pdgIds = (PDG.proton, PDG.proton),
        pz = (6500., 6500.),
        structureFunctions = cepgen.StructureFunctions.luxLike,
    ),
    outKinematics = cepgen.Parameters(
        pt = (0.,),
        energy = (0.,),
        eta = (-2.5, 2.5),
        mx = (1.07, 1000.),
    ),
)

generator = _gen.clone(
    numEvents = 100000,
    printEvery = 10000,
)

text = cepgen.Module('text',
    histVariables={
        'm(4)': cepgen.Parameters(xbins=[float(bin) for bin in range(200, 1500, 50)]),
        'pt(7):pt(8)': cepgen.Parameters(xrange=(0., 750.), yrange=(0., 750.), log=True)
    }
)
dump = cepgen.Module('dump', printEvery = generator.printEvery)
output = cepgen.Sequence(text, dump)
