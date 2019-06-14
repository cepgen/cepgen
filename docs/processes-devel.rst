.. _processes-devel:

=====================
Processes development
=====================

In CepGen, all processes are defined as an object derivating from the following base class:

.. doxygenclass:: cepgen::proc::GenericProcess
   :members:

:math:`k_{\rm T}`-factorised processes
--------------------------------------

Since its version 0.9, CepGen handles the transverse-momentum dependent factorisation of two-photon processes.
For this purpose, a :class:`cepgen::proc::GenericKTProcess` helper derivated-class of the earlier is introduced to allow the photon fluxes part to be transparent to the process developper.

.. doxygenclass:: cepgen::proc::GenericKTProcess
   :outline:

Fortran interface
-----------------

For development or testing purposes (and increased flexibility in general), Fortran definitions of physics processes can be fed and used by CepGen for cross section computation and event generation.
In this page a summary and hands-on example of such a Fortran implementation and linkage is described.

Before other things, you will have to provide the definition of your process and its topology for CepGen to handle it properly.
This requires you to fill in a set of predefined common blocks shared between your process definition and the core CepGen instance.

All these new processes are then linked to CepGen through the construction and association of each matrix element definition and interface to the following object:

.. doxygenclass:: cepgen::proc::FortranKTProcess
   :outline:

In a later paragraph, this linking recipe will be described.

Output event kinematics
~~~~~~~~~~~~~~~~~~~~~~~

This common block lets CepGen access the whole event structure and kinematics information.
Beside the outgoing beam particles 4-momenta ``px`` and ``py``, it contains the PDG id and 4-momentum of all central system particles in ``pc``.

For this latter, a maximum multiplicity of 10 particles is handled by default in CepGen.
Should this maximal value not be sufficient for your implementation purposes, please get in touch with us.

Process definition
~~~~~~~~~~~~~~~~~~

This block specifies all useful run-level parameters, such as:

* ``icontri``, the kinematics mode considered (1 = elastic, 2-3 = single-dissociative, 4 = double-dissociative),
* ``iflux1``, ``iflux2``, the type of :math:`k_{\rm T}`-dependent flux used to parameterise the incoming partons’ kinematics,
* ``imethod``, the type of matrix element considered (for flexibility),
* ``sfmod``, the structure functions modelling used,
* ``pdg_l``, the PDG id of central particles,
* ``inp1``, ``inp2``, the longitudinal beam momenta (in GeV/c).

Additionally, for heavy ions initial states, the following parameters are used:

* ``a_nuc1``, ``a_nuc2``, the mass numbers of positive- and negative-z incoming beam particles,
* ``z_nuc1``, ``z_nuc2``, the atomic numbers of positive- and negative-z incoming beam particles,

:math:`k_{\rm T}`-factorisation kinematics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This common block holds the following :math:`k_{\rm T}`-factorisation processes kinematics variables (generated and filled by the core CepGen instance):

* ``q1t`` and ``q2t``, the 2-norm of transverse incoming parton virtualities,
* ``phiq1t`` and ``phiq2t``, the equivalent azimutal angle of these virtualities,
* ``y1`` and ``y2`` the outgoing system particles’ rapidities,
* ``ptdiff`` and ``phiptdiff`` the 2-norm and azimutal angle of the outgoing system’s transverse momentum difference,
* ``am_x`` and ``am_y`` the dissociative outgoing beams’ invariant masses (or :math:`m_p` for elastic parton emission from proton).

Phase space cuts
~~~~~~~~~~~~~~~~

A limited set of pre-registered phase space cuts (one boolean switch, a lower and an upper value to apply) may be found in this interfacing library, such as:

* single particles :math:`p_{\rm T}`, energy, pseudo-rapidity :math:`\eta`,
* central system invariant mass, :math:`p_{\rm T}`,
* central system correlations: :math:`\Delta y`.

Should this collection not be sufficient for your purposes, please contact us.

Physics constants
~~~~~~~~~~~~~~~~~

This block introduces all useful and standard physical constants in double precision to help the process definition: :math:`m_p`, GeV²/barn, :math:`\pi`, :math:`\alpha_{\rm em}`.

Matrix element definition
~~~~~~~~~~~~~~~~~~~~~~~~~

The core definition of your process is to be implemented within a new ``.f`` file you will store in the ``External/Processes/`` directory of your local CepGen install.

Currently, the only following Fortran process is defined: ``nucl_to_ff``

For illustration in this documentation page, we will be using a dummy process name, such as ``dummy_f77_process``.
For the sake of simplicity in your processes definition bookkeeping we suggest you to match the filename to the matrix element weighting function name.
In our example, this corresponds to ``External/Processes/dummy_f77_process.f``

We expect your weighting function to follow the C-equivalent signature:

.. code-block:: c

   double dummy_f77_process_(void);

This might get translated, for instance, into the following minimal working example (here using the free form Fortran coding convention):

.. code-block:: fortran

   function dummy_f77_process
   implicit none
   double precision dummy_f77_process
   !--------------------------------------------------------------------------
   ! CepGen overhead
   !--------------------------------------------------------------------------
   include 'KTBlocks.inc' ! mandatory, include the kinematics common blocks
   call CepGen_print      ! optional, display some run parameters information
   !--------------------------------------------------------------------------
   ! end of overhead, beginning of process definition
   !--------------------------------------------------------------------------
   dummy_f77_process = 1.D0 ! placeholder, your actual definition is to be
                            ! implemented here
   !--------------------------------------------------------------------------
   ! end of process definition
   !--------------------------------------------------------------------------
   return
   end

With the kinematics common blocks defined through the include statement described above.

Helper tools
~~~~~~~~~~~~

To ease the work of the process developer, we provide utilitaries for the evaluation of unintegrated parton fluxes and other physics quantities through the CepGen core C++ library methods.
The external methods exposed to the Fortran process through the include statement above are:

.. doxygenfunction:: cepgen_kt_flux_
.. doxygenfunction:: cepgen_kt_flux_hi_
.. doxygenfunction:: cepgen_particle_charge_
.. doxygenfunction:: cepgen_particle_mass_

These can be translated in the following Fortran subroutines/functions definitions:

.. code-block:: fortran

  external CepGen_kT_flux,CepGen_kT_flux_HI
  double precision CepGen_kT_flux,CepGen_kT_flux_HI
  external CepGen_particle_charge,CepGen_particle_mass
  double precision CepGen_particle_charge,CepGen_particle_mass

Overall linking
~~~~~~~~~~~~~~~

To interface your process to CepGen, edit the ``External/Processes/CepGenWrapper.cpp`` file to link your function implementation to the core processes module.

The following macros are used to declare the function name and link it to CepGen.

.. doxygendefine:: DECLARE_FORTRAN_FUNCTION
.. doxygendefine:: REGISTER_FORTRAN_PROCESS

.. note:: No single, nor double trailing underscore (_) is required for this name registration.

For the latter, the signature is the following:

.. code-block:: C

   REGISTER_FORTRAN_PROCESS( my_dummy_process, dummy_f77_process,
                             "this process is only for illustration purposes" )

Or, following this nomenclature:

* the first argument corresponds to the process name (i.e. may be used with the my_dummy_process string in steering cards),
* the second argument is the function name as defined above,
* the third one is a one-line, human-readable process description that will be shown at the run start.

For this version, the following process is already registered:

.. literalinclude:: ../External/Processes/CepGenWrapper.cpp
   :language: cpp
   :caption: List of registered external Fortran processes, as imported from ``External/Processes/CepGenWrapper.cpp``.
