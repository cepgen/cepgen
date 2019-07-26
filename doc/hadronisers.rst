Hadronisation/fragmentation
===========================

The decay of secondary products and fragmentation of outgoing dissociative proton states may be triggered directly at the CepGen level.

Several interfaces to external algorithms are therefore provided in the ``CepGenAddOns`` library, and easily steerable through the ``hadroniser`` block (in `Python </cards-python>`_ cards) or ``HADR`` variable (in `LPAIR-like </cards-lpair>`_ cards).
All modules are derived from a common :cpp:class:`cepgen::hadr::GenericHadroniser` class described below:

.. doxygenclass:: cepgen::hadr::GenericHadroniser
   :outline:

.. container:: toggle

   .. container:: header

      Detailed description

   .. doxygenclass:: cepgen::hadr::GenericHadroniser
      :members:
      :no-link:

----

A list of algorithms currently supported may be found below.


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

