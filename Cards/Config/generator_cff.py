from containers_cfi import Parameters

generator = Parameters(
    numEvents = 100000,
    numPoints = 100,
    printEvery = 10000,
    numThreads = 2,
    treat = True, # smoothing of the integrand
)
