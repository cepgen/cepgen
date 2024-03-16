##
# \file
# \ingroup python
#
# Collection of generator parameters

from .containers_cff import Parameters

generator = Parameters(
    numEvents = 100000,
    numPoints = 100,
    printEvery = 10000,
    numThreads = 2,
)
