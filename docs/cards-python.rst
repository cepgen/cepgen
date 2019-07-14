=====================
Python configurations
=====================

This format allows a better control over all run parameters through an increased granularity.
It also allows to reduce the configuration overhead through Python’s ``import`` statements.
Obviously, it requires the Python C API to be linked against the steering cards utils library.

Two core object types (both inheriting from Python’s ``dict`` containers) are imported in the scripting scope through the general statement:

.. code:: python

   import Config.Core as cepgen

These two types are ``Module``, and ``Parameters``.
The earlier is a subset of the latter, with the first string-type attribute defining the module name.

The common usage for such a module definition is, for instance:

.. code:: python

   import Config.Core as cepgen
   module = cepgen.Module('my_first_module', foo = 'bar')

In `this page <python-containers>`_, one can see an illustrated example of this ``Module``/``Parameters``/``dict`` relation.

Several blocks are required to be defined in a CepGen steering card.
In the following sections, you may find a nonexhaustive (and evolving) list of such attributes.

.. _pdg-block:

``PDG`` parameters block
------------------------

.. doxygennamespace:: Cards::Config::PDG_cfi
   :members:

This :mod:`PDG` ``Parameters`` container object holding the list of particles definitions is parsed by the CepGen core at the initialisation level.
It allows to specify new PDG members and propagate their properties for its usage in all parts of the framework.

The :func:`registerParticle` helper function allows to add such a particle in one single line before the ``process`` block definition.

``integrator`` module block
---------------------------

.. warning:: Under construction

As of version 0.8 of CepGen, the three GSL implementations of the following integration algorithms are supported:

* ``Vegas`` by Lepage :cite:`Lepage:1977sw`,
* ``MISER`` stratified sampling by Press et al. :cite:`Press:1989vk`,
* ``plain`` "hit-and-miss" algorithm.

``generator`` module block
--------------------------

.. warning:: Under construction

``hadroniser`` module block
---------------------------

.. warning:: Under construction

``process`` module block
------------------------

This block comes as a required ``Module`` object defined in the general scope.
Its first feature is to specify the process to account for in the user-defined run.
See the list of processes section of the left hand side menu to find your model of interest.

``process.inKinematics`` parameters block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A ``pz`` Python pair (or list) of floating point numbers allows to specify the two incoming protons’ longitudinal momentum (in GeV).
The ``cmEnergy`` keyword can also be used to define directly the centre of mass energy :math:`\sqrt{s}` of the two incoming beams for symmetric, head-on collisions.
In that latter case, :math:`p _ {z,1-2} = \pm \sqrt{s}/2`.

Equivalently, a ``pdgIds`` pair/list of `integer-type PDG identifiers <http://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf>`_ (complete list handled :ref:`here <pdg-block>`) may be used to control beam particles type.
A default ``pdgIds = (2212, 2212)`` initial state, or equivalently ``(PDG.proton, PDG.proton)``, is used.

.. doxygenclass:: Cards::Config::StructureFunctions_cfi::StructureFunctions
   :members:

The ``structureFunctions`` attribute specifies the :math:`F _ {2/L}(\xbj,Q^2)` structure function to use in the parameterisation of the incoming photon fluxes.
The name of the structure functions set (see `the complete list here </structure-functions>`_) has to be prepended by ``StructureFunctions``

For instance, the *Suri-Yennie* set may be selected through ``StructureFunctions.SuriYennie``.

``process.outKinematics`` parameters block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The kinematics phase space to be used in the integration and events production can be specified using a set of cuts applied on the matrix element level:

* ``pt``: single central particle transverse momentum range definition,
* ``energy``: single central particle energy range definition,
* ``eta``: single central particle pseudo-rapidity range definition,
* ``rapidity``: single central particle rapidity range definition,
* ``mx``: outgoing excited proton mass range definition,
* ``xi``: outgoing proton fractional longitudinal momentum loss :math:`\xi = \Delta p/p`.

  .. versionadded:: 0.9.2

``process.processParameters`` parameters block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This block is a generic placeholder for all process-dependent parameters.
See the description page of each process to get a list of supported parameters to include in this collection.

``output`` module block
---------------------------

.. warning:: Under construction

.. _configuration-card-example:

Configuration card example
--------------------------

The generation of 100k single-dissociative :math:`\gg{\mu^+\mu^-}` events at 13 TeV with the `LPAIR matrix element <processes-lpair>`_ implementation with the following phase space cuts:

* :math:`\pt(\mu^\pm)>` 25 GeV, :math:`\lvert\eta(\mu^\pm)\rvert<` 2.5
* 1.07 $< M_X <$ 1000 GeV

can be steered using the following card:

.. code:: python

   import Config.Core as cepgen
   from Config.Integration.vegas_cff import integrator
   from Config.generator_cff import generator as gentmpl
   from Config.PDG_cfi import PDG

   process = cepgen.Module('lpair',
       processParameters = cepgen.Parameters(
           mode = cepgen.ProcessMode.InelasticElastic, # single-dissociation
           pair = PDG.muon, # or, equivalently, 13
       ),
       inKinematics = cepgen.Parameters(
           pz = (6500., 6500.), # or cmEnergy = 13.e3,
           structureFunctions = cepgen.StructureFunctions.SuriYennie,
       ),
       outKinematics = cepgen.Parameters(
           pt = (25., ),
           energy = (0., ),
           eta = (-2.5, 2.5),
           mx = (1.07, 1000.),
       )
   )

   generator = gentmpl.clone(
       numEvents = 1e5,
   )

   output = cepgen.Module('lhef',
       filename = 'lpair-example.lhef',
   )

This configuration is equivalent to the *LPAIR card* shown `here <cards-lpair#configuration-card-example>`_.
