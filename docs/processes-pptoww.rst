.. title:: kT-factorised two-photon production of gauge boson pair

:math:`\kt`-factorised :math:`\ggww`
====================================

The two-photon production of gauge boson pairs process, i.e. :math:`pp \rightarrow p^{(\ast)}(\ggww)p^{(\ast)}`, is only defined using the :math:`\kt`-factorisation approach.

A full review of this process may be found in :cite:`Luszczak:2018ntp`.

Process-specific options
------------------------

``method``
~~~~~~~~~~

This integer allows to switch between the two matrix element definitions:

* ``0``: on-shell amplitude :cite:`Denner:1995jv`,
* ``1``: off-shell amplitude :cite:`Nachtmann:2005en,Nachtmann:2005ep`.

``mode``
~~~~~~~~

This enumeration allows to specify the kinematic regime to generate and the size of the phase space to perform the integration.
It can take the following values:

* ``ProcessMode.ElasticElastic := 1``: elastic emission of photons from the incoming protons (default value if unspecified),
* ``ProcessMode.ElasticInelastic := 2 / ProcessMode.InelasticElastic := 3``: elastic scattering of one photon and an inelastic/semi-exclusive emission of the other photon, resulting in the excitation/fragmentation of the outgoing proton state,
* ``ProcessMode.InelasticInelastic := 4``: both protons fragmented in the final state.

``polarisationStates``
~~~~~~~~~~~~~~~~~~~~~~

This enumeration lets you switch between all combinations of polarisation states to be included in the matrix element.
It can take the following values:

* ``0``: all contributions,
* ``1``: longitudinal-longitudinal,
* ``2``: longitudinal-transverse,
* ``3``: transverse-longitudinal,
* ``4``: transverse-transverse.

.. doxygenclass:: cepgen::proc::PPtoWW
