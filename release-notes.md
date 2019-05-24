# Release notes

## v0.9.4 (9 May 2019)
* pptoff polarisation terms may be steered for off-shell ME
* Functional is now handling exprtk in addition to &mu;Parser
* Refactored exceptions and tests
* Improved accessors for integrator/grid definition
* Better memory management for parameters/processes/hadronisers definition

## v0.9.1-3 (29-30 Nov 2018)
* Fixes in includes chain
* Fix in LPAIR cards parser, simplified cards handler interface
* External Fortran processes now defined through a function instead of a subroutine
* Fix in Pythia8 interface parton colours definition and resonances decay
* Structure functions interface modified

## v0.9 (28-29 Nov 2018)
* Fixed memory leaks in Python configuration files handler
* First HI process! (and first Fortran process interfaced)
* Polarisation states may be specified for the off-shell &gamma;&gamma; &rarr; W<sup>+</sup>W<sup>-</sup> ME
* New structure functions definitions: ALLM (+subsequent), MSTW grid, CLAS, LUXlike, LHAPDF
* Introduced R-ratio definitions for F<sub>L</sub> calculation from F<sub>2</sub>
* New cuts definition
* Better memory management

## v0.8 (24 Aug 2017)
* Major refactoring of the code
* Dropped the hadronisers and other unfinished processes definition for clarity
* Several bug fixes and memory management improvements through the usage of smart pointers / C++11 migration
* Taming functions introduced

## v0.7 (3 May 2016)
* Introduced a generic templating class for the k<sub>T</sub> factorisation method
* Added a placeholder for the &gamma;&gamma; &rarr; W<sup>+</sup>W<sup>-</sup> process computation
