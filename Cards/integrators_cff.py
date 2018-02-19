import Cards.utils_cfi as cg

vegas = cg.Parameters(
  algorithm = "Vegas",
  numIntegrationCalls = 500000,
  numIntegrationIterations = 10,
  numPoints = 100,
)

miser = cg.Parameters(
  algorithm = "MISER",
  numIntegrationCalls = 1000000,
  numPoints = 100,
)
