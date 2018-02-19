import Cards.utils_cfi as cg

vegas = dict(
  algorithm = "Vegas",
  numIntegrationCalls = 500000,
  numIntegrationIterations = 10,
  numPoints = 100,
)

miser = dict(
  algorithm = "MISER",
  numIntegrationCalls = 1000000,
  numPoints = 100,
)
