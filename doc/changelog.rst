Changelog
=========

.. changelog::
    :version: 1.0.0alpha2
    :released: Fri Apr 23 2021

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

