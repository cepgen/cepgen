.. title:: kT-factorised two-photon production of fermion pair

:math:`\kt`-factorised :math:`\ggff` process
============================================

The photon transverse momentum-dependant description of this process was previously developed in ``PPtoLL``, and described in :cite:`daSilveira:2014jla`.
It is defined as the ``pptoll`` process in CepGen.

Process-specific options
------------------------

``method``
~~~~~~~~~~

This integer allows to switch between the two matrix element definitions:

* ``0``: on-shell amplitude,
* ``1``: off-shell amplitude.

``mode``
~~~~~~~~

This enumeration allows to specify the kinematic regime to generate and the size of the phase space to perform the integration.
It can take the following values:

* ``ProcessMode.ElasticElastic := 1``: elastic emission of photons from the incoming protons (default value if unspecified),
* ``ProcessMode.ElasticInelastic := 2 / ProcessMode.InelasticElastic := 3``: elastic scattering of one photon and an inelastic/semi-exclusive emission of the other photon, resulting in the excitation/fragmentation of the outgoing proton state,
* ``ProcessMode.InelasticInelastic := 4``: both protons fragmented in the final state.

``pair``
~~~~~~~~

This integer value allows the end-user to specify the PDG identifier of the fermion pair to be produced in the final state.

It can take the following values:

* ``PDG.electron := 11``: :math:`e^+e^-` pair production
* ``PDG.muon := 13``: :math:`\mu^+\mu^-` pair production
* ``PDG.tau := 15``: :math:`\tau^+\tau^-` pair production
* ``PDG.down``, ``PDG.up``, ``PDG.strange``, ``PDG.charm``, ``PDG.bottom``, ``PDG.top`` (or equivalently ``1-6``): quark pair production.

Full object reference
---------------------

.. doxygenclass:: PPtoFF
   :private-members:
   :undoc-members:
