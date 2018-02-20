import Cards.Core as cepgen

vegas = cepgen.Module('Vegas',
    numIntegrationCalls = 500000,
    numIntegrationIterations = 10,
    numPoints = 100,
)

miser = cepgen.Module('MISER',
    numIntegrationCalls = 1000000,
    numPoints = 100,
)
