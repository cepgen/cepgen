[![build status](https://gitlab.cern.ch/lforthom/cepgen/badges/master/build.svg)](https://gitlab.cern.ch/lforthom/cepgen/commits/master) [![coverage report](https://gitlab.cern.ch/lforthom/cepgen/badges/master/coverage.svg)](https://gitlab.cern.ch/lforthom/cepgen/commits/master)

Release notes:
--------------

V0.8 (24 Aug 2017)
* Major refactoring of the code
* Dropped the hadronisers and other unfinished processes definition for clarity
* Several bug fixes and memory management improvements through the usage of smart pointers / C++11 migration

V0.7 (3 May 2016)
* Introduced a generic templating class for the kT factorisation method
* Added a placeholder for the yy -> W+W- process computation
