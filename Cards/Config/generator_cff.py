import Config.Core as cepgen

generator = cepgen.Parameters(
    numEvents = 100000,
    numPoints = 100,
    printEvery = 10000,
    numThreads = 2,
    treat = True,
)
