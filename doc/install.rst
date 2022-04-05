======================
Installation procedure
======================

The recipe detailed below is intended for “regular users”, who downloaded a version of the code from a current stable release.
You may check in parallel the latest upstream version from `the CepGen repository on GitHub <https://github.com/cepgen/cepgen>`_.

.. note::
   For developers, please find a collection of useful procedures in `this page <install-dev>`_.

General and usage-specific dependencies
---------------------------------------

CepGen was designed to require only a limited set of external dependencies for its core features.

Among the mandatory dependencies,

* `CMake <https://cmake.org/>`_ for the management of the build process,
* `GSL <https://www.gnu.org/software/gsl/>`_, for the integrators definition and general utilitaries.
  A version greater than or equal to ``2.1`` is recommended for the bilinear spline interpolation capability.
  This latter is used in e.g. the MSTW structure functions and KMR :math:`\kt`-factorised gluon flux definitions.

Depending on your system architecture and distribution, these dependencies may be installed the following way:

.. note:: Debian/Ubuntu

   .. code:: sh

      sudo apt-get install cmake
      sudo apt-get install g++ gfortran
      sudo apt-get install libgsl2 libgsl-dev

.. note:: RHEL/Fedora

   .. code:: sh

      sudo dnf install cmake
      sudo dnf install gcc-c++ gcc-gfortran
      sudo dnf install gsl gsl-devel

It can further be extended through a set of add-ons loaded through shared libraries into the runtime environment.

For instance, it is highly recommended to install the following dependencies to ease your first steps:

* ``CepGenPython``: complex Python version 2 or 3 steering cards parsing ; it requires the ``python-devel`` or ``libpython-dev`` development headers to be built,
* ``CepGenLHAPDF`` with `LHAPDF <https://lhapdf.hepforge.org/>`_ version 5 or above, for definition of partonic-level nucleon structure functions, and its :math:`\alpha_S(Q)` strong coupling constant evolution algorithm.

A set of add-ons specifically designed to ease your workflow when generating and storing events is also provided, given the following libraries are available on your system:

* ``CepGenPythia8``: `Pythia <https://pythia.org>`_ version 8.1 or above, for the various particles decays and excited proton fragmentation, but also the LHEF event format generation,
* ``CepGenPythia6``: `Pythia <https://pythia.org/pythia6>`_ version 6, for the legacy steering of particles decays and excited proton fragmentation,
* ``CepGenHepMC2`` and ``CepGenHepMC3`` for the `HepMC <https://hepmc.web.cern.ch/hepmc/>`_ version ≥ 2, to handle its various ASCII output formats,

  .. note::
     From version 3 on, a LHEF output is also supported.
     If not present, this latter format can be handled through the Pythia 8 interface instead.

* ``CepGenROOT``: for the `ROOT <https://root.cern.ch/>`_ ntuples building/histogramming/plotting capabilities widely used in HEP.

CepGen installation
-------------------

Start by downloading the latest release on `our downloads page </downloads>`_.
Unpack the sources in a location referred here as: ``$CEPGEN_SOURCES``.

For instance, for a CepGen version ``x.y.z`` package:

.. code:: sh

   tar xvfz cepgen-x.y.z.tar.gz
   cd cepgen
   export CEPGEN_SOURCES=`pwd -P`

Once you have set up the sources and downloaded/installed the required dependencies, create your building environment (here, we will use the general ``CMake`` convention ``$CEPGEN_SOURCES/build``).
The compilation is hence done with:

.. code:: sh

   mkdir $CEPGEN_SOURCES/build && cd $CEPGEN_SOURCES/build
   cmake $CEPGEN_SOURCES # or `cmake ..`, path to the sources make
   # optionally, add -jN (N=number of parallel threads for compilation)

This compilation will build a collection of required sub-libraries to be linked against any executable built on top of CepGen.:

* ``libCepGenCore`` contains all physics constants, calculators, helpers, along with "non-physics" standard objects implementation, and nucleon structure function calculators objects,
* ``libCepGenProcesses`` contains all processes definitions and implementations,
* ``libCepGenEvent`` holds the definition of events and subleading particles objects (useful for analyses of CepGen outputs),
* ``libCepGenAddOns`` provides a set of helper tools for the interfacing with external applications,
* ``libCepGenCards`` for the input cards definition and handling part.

.. note::
   If your usage requires the import of CepGen libraries and includes in your standard ``PATH``, e.g. for the purpose of interfacing library development, run (as root)

   .. code:: sh

      make install

   This will copy all required headers into the local includes directory (e.g. ``/usr/local/include``), and copy the shared objects into the library path (e.g. ``/usr/local/lib64`` or ``/usr/local/lib``).

As described `here </usage>`_, several test executables can be linked against the CepGen libraries.

.. note::
   You may build these executables using the ``make`` command.
   For instance, ``make cepgen``. The test executable will then be located in the ``test/`` directory.
   You may run it using, for instance (in ``cepgen-dev/build/``):

   .. code:: sh

      ./test/cepgen <path to your steering card>

.. .. doxygengroup:: Executables
