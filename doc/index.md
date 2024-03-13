% CepGen documentation master file, created by
% sphinx-quickstart on Thu Jun 13 15:50:31 2019.
% You can adapt this file completely to your liking, but it should at least
% contain the root `toctree` directive.

# CepGen, a generic central exclusive processes events generator

Welcome to the home page of the CepGen Monte Carlo event generator, released under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.en.html).

This documentation page attempts to summarise the various processes and configuration parameters currently, and soon to be handled by this piece of software.

Please use the reference {cite:p}`Forthomme:2018ecc`: for citation:

> [Comput. Phys. Commun. 271 (2022) 108225](https://doi.org/10.1016/j.cpc.2021.108225). [arXiv:1807.06059 \[hep-ph\]](https://arxiv.org/abs/1808.06059).

This page is constantly evolving, but do not hesitate to contact us through the project’s [GitHub issues page](https://github.com/cepgen/cepgen/issues).

```{toctree}
:caption: 'Contents:'
:hidden: true
:maxdepth: 1

Downloads <https://github.com/cepgen/cepgen/releases>
changelog
install
usage
zz-bibliography
```

```{toctree}
:caption: 'Configuration:'
:hidden: true
:maxdepth: 1

cards-lpair
cards-python
event-modifiers
output-formats
```

```{toctree}
:caption: 'Processes:'
:hidden: true
:maxdepth: 1

lpair <processes-lpair>
pptoff/pptoll <processes-pptoff>
pptoww <processes-pptoww>
```

```{toctree}
:caption: 'Advanced usage:'
:hidden: true
:maxdepth: 1

processes-devel
structure-functions
event
integrators
random-generators
raw-modules
```
