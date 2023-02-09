Changelog
=========

Latest revisions
^^^^^^^^^^^^^^^^

.. changelog::
    :version: 1.0.2
    :released: Mon Aug 22 2022

    .. change::
        :tags: core

        An analytic integrators base class was introduced, including the few GSL-based implementations already present since the earlier version, along with Boost- and ROOT-based integrators.

.. changelog::
    :version: 1.0.1
    :released: Tue May 24 2022

    .. change::
        :tags: core

        This version fixes a few coding flaws and allows build using Clang versions 12+.
        Build on LXPLUS/environment capable of accessing LCG releases is made safer through the use of standard LCG v101.

    .. change::
        :tags: utils

        A wrapper to a good fraction of the GSL one-dimensional integration algorithms is introduced. This paves the ground for the future (integrated) collinear fluxes computation and several underlying utilities.
        The ROOT drawer is also made safer through the use of standard colours if the number of subplots exceeds the pool of CepGen-themed colours.
        Additionally, the command lines arguments parser now allows limits to be specified through a min,max couple.

.. changelog::
    :version: 1.0.0
    :released: Thu May 12 2022

    .. change::
        :tags: core

        Several utilities are now converted into the post-0.9.X modules schema.
        This allows to delegate a few definitions in the runtime loading of all libraries compiled against CepGen, thus reducing the overhead of dependencies for the CepGen core library.

        A new parameters documentation system was introduced to list all possible keys and their default/expected values and ease the user-interaction with all module parameters.

    .. change::
        :tags: processes

        ``PPtoWW`` now includes more (incl. anomalous) matrix element implementations listed in `Eur.Phys.J.C45:679-691,2006 <https://doi.org/10.1140/epjc/s2005-02450-3>`_.

    .. change::
        :tags: strfun

        New hybrid Kulagin-Barinov structure functions, as implemented in `Phys. Rev. C 105 (2022) 045204 <https://doi.org/10.1103/PhysRevC.105.045204>`_.

    .. change::
        :tags: utils

        Added a set of utilities to ease the drawing of 1- and 2-dimensional graphs and histograms.
        In addition to the "standard" text-based renderer, several libraries are interfaced to generate their output (ROOT, YODA, Gnuplot, Matplotlib, Topdrawer).

    .. change::
        :tags: external

        Python cards steering et al. interface is now stripped off the core ``CepGen`` library into a dedicated ``CepGenPython`` library.
        It now includes a functional parser and output configuration producer.

        HepMC inteface is now further splitted between its pre3 and 3+ versions.
        This allows to ease the interfacing between CepGen event content and several libraries accepting a HepMC2 or HepMC3 event content.

        Added a Photos++ and a Tauola++ algorithms interface for event modification.
        Included a testing suite for e.g. Pythia 6 steering through its CepGen interface.

.. changelog::
    :version: 1.0.0alpha2
    :released: Fri Apr 23 2021

.. changelog::
    :version: 0.9.9
    :released: Tue Dec 31 2019

.. changelog::
    :version: 0.9.8
    :released: Wed Oct 16 2019

.. changelog::
    :version: 0.9.7
    :released: Thu Jul 25 2019

    .. change::
        :tags: processes

        Fortran processes can now be fed a generic set of parameters, thanks to additional getter functions

    .. change::
        :tags: output
        :changeset: b8e5927e52, 507f8ccdc8

        Output handlers may now be constructed directly from steering cards, thus enhancing overall modularity.

    .. change::
        :tags: output
        :changeset: d59f3702ca

        New text output handler (raw text output, and ASCII histograms)

    .. change::
        :tags: output
        :changeset: 7f982e3a3d

        New HepMC ASCII output handler (for HepMC v<3), refactored HepMC event builder in preparation for future developments

    .. change::
        :tags: output
        :changeset: e467dcf1a0, e3b10e3572

        New ROOT histogram collections and ntuple files writers. Dropped the support for the ``cepgen-root`` executable.

    .. change::
        :tags: output
        :changeset: 0f0e541a2f

        Interface to Delphes for the simulation of detectors effects

    .. change::
        :tags: core
        :changeset: 65ae85039c

        Added a helper for the retrieval of events properties through human-readable getters

.. changelog::
    :version: 0.9.6
    :released: Thu Jul 11 2019

    .. change::
        :tags: external
        :changeset: 06ebf75259

        Added support of Pythia6 hadronisation/fragmentation algorithm for legacy tests

    .. change::
        :tags: core
        :changeset: 7c57a24d31, 1c5e353895

        Structure functions parameterisation objects polished

    .. change::
        :tags: output

        New output modes handled for HepMC interfacing module

Latest developments
^^^^^^^^^^^^^^^^^^^

.. collapse:: Show

   .. git_changelog::
      :rev-list: 1.0.2..HEAD
      :filename_filter: CepGen*

