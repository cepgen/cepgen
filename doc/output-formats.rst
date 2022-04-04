Output formats
==============

In CepGen, the event (particles, their parentage and kinematics) is handled through the ``cepgen::Event`` object described `in the dedicated event format page </event>`_.

To ease the user interaction with this object, a few output writers (defined here as "handlers") are given as examples.
All handlers are defined as modules derivating from the following abstract base class:

.. doxygenclass:: cepgen::io::ExportModule
   :outline:

.. container:: toggle

   .. container:: header

      Detailed description

   .. doxygenclass:: cepgen::io::ExportModule
      :members:
      :no-link:

----

In this page you will find a list of all currently supported output formats, covering a broad spectrum of usages, both in the phenomenological and experimental communities.
Please note that this list is under constant evolution, you may contact us with requests for additional interfacing capabilities.

``dump``
--------

A simple text-based event dumper, useful for debugging the process and its kinematics, is steered using the :cpp:class:`cepgen::io::TextEventHandler` module.

  .. doxygenclass:: cepgen::io::TextEventHandler
     :outline:

``lhef``
--------

This output format handles the conversion into the `Les Houches standard definition <https://en.wikipedia.org/wiki/Les_Houches_Accords>`_.
Currently, two implementations of this export module exist:

- a ``Pythia 8`` LHEF output module (described `here <http://home.thep.lu.se/~torbjorn/pythia82html/LesHouchesAccord.html>`_) as the default handler, :cpp:class:`cepgen::io::LHEFPythiaHandler`,

  .. doxygenclass:: cepgen::io::LHEFPythiaHandler
     :outline:

- a ``HepMC (v≥3)`` implementation, if the earlier is not found in the standard libraries path: :cpp:class:`cepgen::io::LHEFHepMC3Handler`.

``hepmc2``, ``hepmc2_ascii``, ...
---------------------------------

.. doxygenclass:: cepgen::io::HepMC2Handler
   :outline:

This handler allows to translate the CepGen event record into one (or multiple) implementation(s) of the version 2 of the `HepMC <http://hepmc.web.cern.ch/hepmc>`_ :cite:`Dobbs:2001ck` ASCII output format.
By default, this version is used in older releases. It allows a ``hepmc2`` output format to be supported.

``hepmc``, ``hepmc_root``, ``hepevt``, ...
------------------------------------------

.. doxygenclass:: cepgen::io::HepMC3Handler
   :outline:

This handler allows to translate the CepGen event record into one (or multiple) implementation(s) of the version 3 of the  `HepMC <http://hepmc.web.cern.ch/hepmc>`_ :cite:`Dobbs:2001ck` ASCII output format.

By default, the version 3 of the file format is chosen for versions of ``HepMC`` starting from ``v3.1.0``.
It may be updated with future derivatives of `the HepMC writer base class <http://hepmc.web.cern.ch/hepmc/classHepMC3_1_1Writer.html>`_.

Alternatively, as from this version ``3.1.0`` of ``HepMC``, the following output formats are also handled:

- a ``hepevt`` ASCII format using the :cpp:class:`HepMC3::WriterHEPEVT` handler,
- a ``hepmc_root`` format using the :cpp:class:`HepMC3::WriterRoot` export module,
- a ``hepmc_root_tree`` using the :cpp:class:`HepMC3::WriterRootTree` module.

``promc``
---------

.. versionadded:: 0.9.8

.. doxygenclass:: cepgen::io::ProMCHandler
   :outline:

The support has been added for the `ProMC <http://jwork.org/wiki/PROMC>`_ highly compressed output format.

``vars``
--------

.. versionadded:: 1.0.0

.. doxygenclass:: cepgen::io::TextVariablesHandler
   :outline:

This simplest case of an output module allows to generate a **generic (ASCII) output format** along with **raw text histograms** of kinematic variables, fully configurable by the user.
Using the Python steering cards definition, a list of variables to be stored is defined through the ``variables`` list/array of string-typed definition.

For this **text output format**, the default behaviour is storing one event per line with variables separated with an user-parameterisable separator (``separator`` string parameter, default is the standard tabulation ``\t``).

The variable (here, ``var`` is used as an example) may be defined using the three following conventions:

- ``var`` for event-level information (e.g. diffractive outgoing proton state multiplicity)
- ``var(role)`` for the retrieval of a single particle with a given role.

  This latter may be one of the followings:
   - ``ib1`` and ``ib2`` (resp. ``ob1`` and ``ob2``) for the incoming (resp. outgoing) beam kinematics,
   - ``pa1`` and ``pa2`` for the parton/initiator particle emitted from the first/second incoming beam particle respectively,
   - ``cs`` for the two-parton/initators system, and
   - ``int`` for any intermediate :math:`s`-channel particle exchange (depending on the process),
- ``var(id)`` for the retrieval of a single particle with a given integer identifier.

As from version ``0.9.7`` of CepGen, the following variables are handled for the particles momentum definition:

.. doxygenvariable:: cepgen::utils::EventBrowser::m_mom_str_

Two extra boolean parameters may also be fed to the module configuration:

- ``saveBanner``, to enable/disable the CepGen banner printout (containing useful information about the process and cuts definition), and
- ``saveVariables``, to show/hide the list of variables used in this file.

As an example, the following ``output`` block may be used for the ``lpair`` process:

.. code:: python

   output = cepgen.Module('text',
       filename = 'test.txt',
       variables = [
           'm(4)', 'pt(cs)', 'pt(6)'
       ],
       saveBanner = False,
       saveVariables = True,
       separator = ' ', # single space
   )

``text``
--------

.. versionadded:: 1.0.0

.. doxygenclass:: cepgen::io::IntegratedEventVariablesHandler
   :outline:

This simplest case of an output module allows to generate **raw text histograms** of kinematic variables, fully configurable by the user.
Using the Python steering cards definition, a dictionary ``histVariables`` of variable-indexed ``cepgen.Parameters`` objects is fed to the ``output`` module.

A valid implementation of such objects requires a set of attributes depending on the type of distribution requested by the user.
If the variable string contains one ``:``, a 2D distribution is automatically booked for the two variables around it.
Otherwise, a 1D distribution is assumed.

These attributes are, namely for 1-dimensional histograms:

- a number of bins ``nbins``, or ``nbinsX``, and
- a range (``low`` and ``high``) of interest for the variable, or a set of bins in a ``xbins`` Python list

and for 2-dimensional distributions:

- the two ``nbinsX`` and ``nbinsY`` number of bins, and
- the two ranges (``lowX`` and ``highX``, and ``lowY`` and ``highY``) of interest for the variables, or equivalently one or two sets of bins in ``xbins``/``ybins`` lists.

As an example, equivalently to ``vars`` output defined above, the following ``output`` block may be used for the ``lpair`` process:

.. code:: python

   output = cepgen.Module('text',
       histVariables = {
          # 1D histogram (pt of central system)
          'pt(4)': cepgen.Parameters(nbins=10, low=0., high=20.),
          # 1D histogram (outgoing proton mass)
          'm(5)': cepgen.Parameters(xbins=[float(bin) for bin in range(0, 20, 1)]),
          # 2D histogram (central system rapidity vs. mass)
          'y(cs):m(cs)': cepgen.Parameters(nbinsX=10, lowX=-10., highX=10.,
                                           nbinxY=10, lowY=0., highY=400.),
          # 1D histogram (event generation time)
          'tgen': cepgen.Parameters(nbins=100, low=0., high=1.e-5),
       },
       save = False,
       show = True
   )

``root``, ``root_tree``
-----------------------

.. versionadded:: 0.9.7
.. note:: Previously used in dedicated test executables, resp. ``test_distributions`` and ``cepgen-root``.

These two modules module allow to produce a **ROOT** :cite:`Brun:1997pa` **file** containing either:

- a list of histograms (stored as ROOT :cpp:class:`TH1D` objects) provided as an input for the earlier:

  .. doxygenclass:: cepgen::io::ROOTHistsHandler
     :outline:

- or a set of **events** and **run information** (stored as ROOT :cpp:class:`TTree` objects) for the latter:

  .. doxygenclass:: cepgen::io::ROOTTreeHandler
     :outline:

The histogramming utilitary follows the same procedure as introduced for the :cpp:class:`cepgen::io::TextHandler` module above to define the histograms list.

As an example, the following ``output`` block may be used:

.. code:: python

   output = cepgen.Module('root',
       filename = 'output.hists.root',
       variables = {
          'pt(4)': cepgen.Parameters(nbins=10, low=0., high=20.),
          'm(5)': cepgen.Parameters(nbins=10, low=0., high=100.),
          'y(cs)': cepgen.Parameters(nbins=10, low=-10., high=10.),
          'tgen': cepgen.Parameters(nbins=100, low=0., high=1.e-5),
       },
   )

The tree handler may be used in parallel to the two :cpp:class:`ROOT::CepGenRun` and :cpp:class:`ROOT::CepGenEvent` helper reader objects for a compact analysis workflow:

  .. container:: toggle

     .. container:: header

        Detailed description

     .. doxygenclass:: ROOT::CepGenRun
        :members:
     .. doxygenclass:: ROOT::CepGenEvent
        :members:

``delphes``
-----------

.. versionadded:: 0.9.7
.. doxygenclass:: cepgen::io::DelphesHandler
   :outline:

An interface to the `Delphes <https://cp3.irmp.ucl.ac.be/projects/delphes>`_ :cite:`deFavereau:2013fsa` fast simulation framework is provided if the library is installed on the system.

Beside the usual ``filename`` flag specifying the file name Delphes will use for its output, a path to the `Tcl <https://www.tcl.tk/>`_ configuration card is also required to steer the output module through the ``inputCard`` string parameter.

Please refer to the Delphes manual and comprehensive list of examples for more information on the steering of the detector simulation.

