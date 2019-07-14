Output formats
==============

Internal event record
---------------------

In CepGen, the event (particles, their parentage and kinematics) is handled through the ``cepgen::Event`` object:

.. doxygenclass:: cepgen::Event
   :members:

To ease the user interaction with this object, a few output writers (defined here as "handlers") are given as examples.
All handlers are defined as modules derivating from the following abstract base class:

.. doxygenclass:: cepgen::io::GenericExportHandler
   :members:

In this page you will find a list of all currently supported output formats.

``lhef``
--------

This output format handles the conversion into the `Les Houches standard definition <https://en.wikipedia.org/wiki/Les_Houches_Accords>`_.
Currently, two implementations of this export module exist:

- a ``Pythia 8`` LHEF output module (described `here <http://home.thep.lu.se/~torbjorn/pythia82html/LesHouchesAccord.html>`_) as the default handler, :cpp:class:`cepgen::io::Pythia8LHEFHandler`,
- a ``HepMC (v≥3)`` implementation, if the earlier is not found in the standard libraries path: :cpp:class:`cepgen::io::HepMCLHEFHandler`.

``hepmc``
---------

.. doxygenclass:: cepgen::io::HepMCHandler
   :outline:

This handler allows to translate the CepGen event record into one (or multiple) implementation(s) of the `HepMC <http://hepmc.web.cern.ch/hepmc>`_ ASCII output format, as the ``hepmc2`` and/or ``hepmc3`` derivatives may be used according to the version installed on the system.

By default, the version 3 of the file format is chosen for versions of ``HepMC`` starting from ``v3.1.0``, and the version 2 is used in older releases.
For the earlier, a ``hepmc2`` output format is however also supported.
It may be updated with future derivatives of `the HepMC writer base class <http://hepmc.web.cern.ch/hepmc/classHepMC3_1_1Writer.html>`_.

Alternatively, as from this version ``3.1.0`` of ``HepMC``, the following output formats are also handled:

- a ``hepevt`` ASCII format using the :cpp:class:`HepMC3::WriterHEPEVT` handler,
- a ``hepmc_root`` format using the :cpp:class:`HepMC3::WriterRoot` export module,
- a ``hepmc_root_tree`` using the :cpp:class:`HepMC3::WriterRootTree` module.

``text``
--------

.. versionadded:: 0.9.7

.. doxygenclass:: cepgen::io::TextHandler
   :outline:

This last module allows to generate a generic (ASCII) output format, fully configurable by the user.
Using the Python steering cards definition, a list of variables to be stored is defined through the ``variables`` list/array of string-typed definition.

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

.. doxygenvariable:: cepgen::io::TextHandler::m_mom_str_

As an example, the following ``output`` block may be used for the ``lpair`` process:

.. code:: python

   output = cepgen.Module('text',
       filename = 'test.txt',
       variables = [
           'm(4)', 'pt(cs)', 'pt(6)'
       ],
   )

