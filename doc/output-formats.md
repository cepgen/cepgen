# Output formats

In CepGen, the event (particles, their parentage and kinematics) is handled through the `cepgen::Event` object described [in the dedicated event format page](/event.md).

To ease the user interaction with this object, a few output writers (defined here as "handlers") are given as examples.
All handlers are defined as modules derivating from the following abstract base class:

```{doxygenclass} cepgen::EventExporter
:outline:
```

````{toggle}
```{doxygenclass} cepgen::EventExporter
:members:
:no-link:
```
````

A full list of the output modules currently supported in CepGen addons, along with their user-steerable parameters, can be found [here](/raw-modules.md#evtout).

______________________________________________________________________

## Event harvester

```{versionadded} 1.0.0
```

```{doxygenclass} cepgen::EventHarvester
:outline:
```

This simplest case of an output module allows to generate **integrated histograms** of kinematic variables, fully configurable by the user.
Using the Python steering cards definition, a dictionary `histVariables` of variable-indexed `cepgen.Parameters` objects is fed to the `output` module.

A valid implementation of such objects requires a set of attributes depending on the type of distribution requested by the user.
If the variable string contains one `:`, a 2D distribution is automatically booked for the two variables around it.
Otherwise, a 1D distribution is assumed.

These attributes are, namely for 1-dimensional histograms:

- a number of bins `nbins`, or `nbinsX`, and
- a range (`low` and `high`) of interest for the variable, or a set of bins in a `xbins` Python list

and for 2-dimensional distributions:

- the two `nbinsX` and `nbinsY` number of bins, and
- the two ranges (`lowX` and `highX`, and `lowY` and `highY`) of interest for the variables, or equivalently one or two sets of bins in `xbins`/`ybins` lists.

As an example, equivalently to `vars` output defined above, the following `output` block may be used for a `text` output histogram with kinematics equivalent to the `lpair` process:

```python
cepgen.Module('text',
    histVariables = {
       # 1D histogram (pt of central system)
       'pt(4)': cepgen.Parameters(nbins=10, xrange=(0., 20.)),
       # 1D histogram (outgoing proton mass)
       'm(5)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 20, 1)]),
       # 2D histogram (central system rapidity vs. mass)
       'y(cs):m(cs)': cepgen.Parameters(nbinsX=10, xrange=(-10., 10.),
                                        nbinxY=10, yrange=(0., 400.)),
    },
    save = False,
    show = True
)
```

To quote a few introduced since v1 of CepGen, the following plotters are handled for harvesters:

- `gnuplot`, on systems which handle a working version of the `gnuplot` executable documented [here](http://www.gnuplot.info/), and for which the command line piper works,
- `matplotlib`, in case the [C++ wrapper](https://matplotlib-cpp.readthedocs.io) of the [matplotlib](https://matplotlib.org/),
- `text`, for the CepGen-specific terminal-based text histogram drawer, introduced as an overlay of the [GSL histogramming](https://www.gnu.org/software/gsl/doc/html/histogram.html) capability.

## Other types of output modules

In this page you will find a list of all currently supported output formats, covering a broad spectrum of usages, both in the phenomenological and experimental communities.
Please note that this list is under constant evolution, you may contact us with requests for additional interfacing capabilities.

### `dump`

A simple text-based event dumper, useful for debugging the process and its kinematics, is steered using the {cpp:class}`cepgen::TextEventHandler` module.

```{doxygenclass} cepgen::TextEventHandler
:outline:
```

### `lhef`

This output format handles the conversion into the [Les Houches standard definition](https://en.wikipedia.org/wiki/Les_Houches_Accords).
Currently, two implementations of this export module exist:

- a `Pythia 8` LHEF output module (described [here](http://home.thep.lu.se/~torbjorn/pythia82html/LesHouchesAccord.html)) as the default handler, {cpp:class}`cepgen::LHEFPythiaHandler`,

  ```{doxygenclass} cepgen::LHEFPythiaHandler
  :outline:
  ```

- a `HepMC (v≥3)` implementation, {cpp:class}`cepgen::LHEFHepMCHandler`.

  ```{doxygenclass} cepgen::LHEFHepMCHandler
  :outline:
  ```

### `hepmc2`, `hepmc2_ascii`, ...

```{doxygenclass} cepgen::HepMC2Handler
:outline:
```

This handler allows to translate the CepGen event record into one (or multiple) implementation(s) of the version 2 of the [HepMC](http://hepmc.web.cern.ch/hepmc) {cite}`Dobbs:2001ck` ASCII output format.
By default, this version is used in older releases. It allows a `hepmc2` output format to be supported.

### `hepmc`, `hepmc_root`, `hepevt`, ...

```{doxygenclass} cepgen::HepMC3Handler
:outline:
```

This handler allows to translate the CepGen event record into one (or multiple) implementation(s) of the version 3 of the  [HepMC](http://hepmc.web.cern.ch/hepmc) {cite}`Dobbs:2001ck` ASCII output format.

By default, the version 3 of the file format is chosen for versions of `HepMC` starting from `v3.1.0`.
It may be updated with future derivatives of [the HepMC writer base class](http://hepmc.web.cern.ch/hepmc/classHepMC3_1_1Writer.html).

Alternatively, as from this version `3.1.0` of `HepMC`, the following output formats are also handled:

- a `hepevt` ASCII format using the {cpp:class}`HepMC3::WriterHEPEVT` handler,
- a `hepmc_root` format using the {cpp:class}`HepMC3::WriterRoot` export module,
- a `hepmc_root_tree` using the {cpp:class}`HepMC3::WriterRootTree` module.

### `promc`

```{versionadded} 0.9.8
```

```{doxygenclass} cepgen::ProMCHandler
:outline:
```

The support has been added for the [ProMC](http://jwork.org/wiki/PROMC) highly compressed output format.

### `vars`

```{versionadded} 1.0.0
```

```{doxygenclass} cepgen::TextVariablesHandler
:outline:
```

This simplest case of an output module allows to generate a **generic (ASCII) output format** along with **raw text histograms** of kinematic variables, fully configurable by the user.
Using the Python steering cards definition, a list of variables to be stored is defined through the `variables` list/array of string-typed definition.

For this **text output format**, the default behaviour is storing one event per line with variables separated with an user-parameterisable separator (`separator` string parameter, default is the standard tabulation `\t`).

The variable (here, `var` is used as an example) may be defined using the three following conventions:

- `var` for event-level information (e.g. diffractive outgoing proton state multiplicity)

- `var(role)` for the retrieval of a single particle with a given role

  These particle roles may be one of the followings:
  : - `ib1` and `ib2` (resp. `ob1` and `ob2`) for the incoming (resp. outgoing) beam kinematics,
    - `pa1` and `pa2` for the parton/initiator particle emitted from the first/second incoming beam particle respectively,
    - `cs` for the two-parton/initators system, and
    - `int` for any intermediate $s$-channel particle exchange (depending on the process),

- `var(id)` for the retrieval of a single particle with a given integer identifier

- `var(role1,role2)` or `var(id1,id2)` for the multi-particles kinematics correlation.

As from version `0.9.7` of CepGen, the following variables are handled for the particles (single, or combined) momentum definition:

```{doxygenvariable} cepgen::utils::EventBrowser::m_mom_str_
```

In addition, particles momenta's correlations can be accessed through the following keywords:

```{doxygenvariable} cepgen::utils::EventBrowser::m_two_mom_str_
```

Two extra boolean parameters may also be fed to the module configuration:

- `saveBanner`, to enable/disable the CepGen banner printout (containing useful information about the process and cuts definition), and
- `saveVariables`, to show/hide the list of variables used in this file.

As an example, the following `output` block may be used for a 2-to-4 process such as `lpair`:

```python
cepgen.Module('vars',
    filename = 'test.txt',
    variables = [
        'm(4)',    # central two-photon/central system mass
        'pt(cs)',  # central two-photon/central system transverse momentum
        'pt(6)'    # first outgoing central system particle transverse momentum
    ],
    saveBanner = False,
    saveVariables = True,
    separator = ' ',  # single space
)
```

### `root_hist`, `root_tree`

```{versionadded} 0.9.7
```

```{note}
Previously used in dedicated test executables, resp. `test_distributions` and `cepgen-root`.
```

These two modules module allow to produce a **ROOT** {cite}`Brun:1997pa` **file** containing either:

- a list of histograms (stored as ROOT {cpp:class}`TH1D` objects) provided as an input for the earlier:

  ```{doxygenclass} cepgen::ROOTHistsHandler
  :outline:
  ```

- or a set of **events** and **run information** (stored as ROOT {cpp:class}`TTree` objects) for the latter:

  ```{doxygenclass} cepgen::ROOTTreeHandler
  :outline:
  ```

The histogramming utilitary follows the same procedure as introduced for the {cpp:class}`cepgen::TextHandler` module above to define the histograms list.

As an example, the following `output` block may be used:

```python
cepgen.Module('root_hist',
    filename = 'output.hists.root',
    variables = {
       'pt(4)': cepgen.Parameters(nbins=10, xrange=(0., 20.)),
       'm(5)': cepgen.Parameters(nbins=10, xrange=(0., 100.)),
       'y(cs)': cepgen.Parameters(nbins=10, xrange=(-10., 10.)),
    },
)
```

The tree handler may be used in parallel to the two {cpp:class}`ROOT::CepGenRun` and {cpp:class}`ROOT::CepGenEvent` helper reader objects for a compact analysis workflow:

````{toggle}
```{doxygenclass} ROOT::CepGenRun
:members:
```
```{doxygenclass} ROOT::CepGenEvent
:members:
```
````

### `delphes`

```{versionadded} 0.9.7
```

An interface to the [Delphes](https://cp3.irmp.ucl.ac.be/projects/delphes) {cite}`deFavereau:2013fsa` fast simulation framework is provided through the `CepGenDelphes` add-on implemented [here](https://github.com/cepgen/cepgen/blob/master/CepGenAddOns/ROOTWrapper/DelphesHandler.cpp).

Beside the usual `filename` flag specifying the file name Delphes will use for its output, a path to the [Tcl](https://www.tcl.tk/) configuration card is also required to steer the output module through the `inputCard` string parameter.

Please refer to the Delphes manual and comprehensive list of examples for more information on the steering of the detector simulation.
