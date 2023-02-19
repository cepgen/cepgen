Event modification algorithms
=============================

The user's needs may sometimes require the modification on an event-by-event basis of the process-generated particles kinematics and relationships,
for instance in the scope of a beam particles remnants (dissociative proton states, for instance) hadronisation/fragmentation or unstable particles decays.

Several interfaces to external algorithms are therefore provided in the ``CepGenAddOns`` library, and easily steerable through the ``eventSequence`` sequential block (in `Python </cards-python>`_ cards) or ``HADR`` variable (in `LPAIR-like </cards-lpair>`_ cards).
All modules are derived from a common :cpp:class:`cepgen::hadr::Hadroniser` class itself derived from a :cpp:class:`cepgen::EventModifier` object, both described below:

.. doxygenclass:: cepgen::EventModifier
   :outline:

.. toggle::

   .. container:: header

      Detailed description

   .. doxygenclass:: cepgen::EventModifier
      :members:
      :no-link:

----

.. doxygenclass:: cepgen::hadr::Hadroniser
   :outline:

.. toggle::

   .. container:: header

      Detailed description

   .. doxygenclass:: cepgen::hadr::Hadroniser
      :members:
      :no-link:

----

A full list of the algorithms and their parameters can be found `here <raw-modules#evtmod>`_.
In particular, a sub-collection of algorithms currently supported is:


``pythia6``
-----------

.. versionadded:: 0.9.6
.. doxygenclass:: cepgen::hadr::Pythia6Hadroniser
   :outline:

This legacy fragmentation module mimicks the original LPAIR Jetset interfacing.
Thus, in dissociative photon emission, this latter is approximated as emitted from a valence quark tied to a diquark system in a beam remnant.
The flavours mixing is performed randomly on an event-by-event basis (with values chosen in :math:`(u,ud_0)`, :math:`(u,ud_1)`, and :math:`(d,uu_1)`).

``pythia8``
-----------

.. warning:: Under construction

.. versionadded:: 0.9
.. doxygenclass:: cepgen::hadr::Pythia8Hadroniser
   :outline:

