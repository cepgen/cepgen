# Monte Carlo numerical integrators

The modularity of CepGen also allows for multiple integration algorithms to be steered and profit for some increased numerical stability of a wise choice of parameters.
Several interfaces to external algorithms are provided in the core and `CepGenAddOns` libraries.

In the [Python](/cards-python.md) cards parsing these can be steered through the `integrator` keyword.
All modules are derived from a common {cpp:class}`cepgen::Integrator` object, described below:

A full list of the algorithms and their parameters can be found [here](/raw-modules.md#integr).

```{doxygenclass} cepgen::Integrator
:outline:
```

Detailed description

````{toggle}
```{doxygenclass} cepgen::Integrator
:members:
:no-link:
```
````

______________________________________________________________________

## GSL Monte Carlo integrator algorithms

Historical integrator algorithms, relying on the GSL implementation of `gsl_monte_xxx` integration routines (see [the GSL manual](https://www.gnu.org/software/gsl/doc/html/montecarlo.html) for a description of all relevent parameters).

```{doxygenclass} cepgen::GSLIntegrator
:private-members:
```

### Vegas integration algorithm

```{doxygenclass} cepgen::VegasIntegrator
:private-members:
```

### Miser integration algorithm

```{doxygenclass} cepgen::MISERIntegrator
```

### Plain integration algorithm

```{doxygenclass} cepgen::PlainIntegrator
:private-members:
```

## ROOT integration algorithms

Interfacing to the two ROOT numerical integration algorithms:

- a one-dimensional integrator (see [ROOT::Math::IntegratorOneDim](https://root.cern.ch/doc/master/classROOT_1_1Math_1_1IntegratorOneDim.html)) ;
- a multidimensional integrator (see [ROOT::Math::IntegratorMultiDim](https://root.cern.ch/doc/master/classROOT_1_1Math_1_1IntegratorMultiDim.html)).

According to the user's request, either of the two objects is populated and configured. The following parameters are to be steered by the end user:

- `type`: integrator algorithm type:
  : - `gauss`, `legendre`, `adaptive`, `adaptiveSingular`, `nonAdaptive`, for one-dimensional integration,
    - `adaptive`, `plain`, `miser`, `vegas` for multidimensional integration. The last three are one-to-one equivalent to the GSL Monte Carlo algorithms described above (except for the interface and parameters definition).
- `absToL`: desired absolute error ;
- `relToL`: desired relative error ;
- `size`: maximum number of sub-intervals.

```{versionadded} 0.9.10
```

```{doxygenclass} cepgen::ROOTIntegrator
:private-members:
```

## Python integration algorithms

Among the newest features is the interfacing between CepGen and any Python numerical integrator, provided the definition of a wrapper function (in Python) called by an interfacing object:

```{versionadded} 1.2.0
```

```{doxygenclass} cepgen::PythonIntegrator
:private-members:
```

The Python interfacing/wrapper function for CepGen numerical integration has the form:

```python
def integrate(f,  # CepGen function to integrate (a [](vector<double>) -> double function)
              num_dim: int,  # number of dimensions to integrate
              num_iter: int,  # number of iterations for integration (if used)
              num_warmup: int,  # number of calls to perform for an integrator warmup (if any)
              num_calls: int,  # number of function calls for each iteration
              limits: list[tuple[float]]=[]  # list of (min, max) limits for each dimension
              ):
    # definition of integration procedure
    # ...
    return (average, standard_deviation)
```
