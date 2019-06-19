LPAIR-like steering cards
=========================

The second (and simplest one), inherited from ``LPAIR`` and ``PPtoLL``, only allows to set a well-defined set of parameters through a *key → value* scheme.
All keys currently handled by this parser are listed below, along with their default value (if not retrieved from the steering card).

General parameters
------------------

* ``PROC``
    Process to generate (complete list `here </processes>`__)

* ``MODE`` (default: ``1`` = elastic-elastic)
    Subprocess’ mode

* ``NGEN`` (default: ``0``)
    Number of unweighted events to generate

* ``NTRT`` (default: ``1``)
    Flag to specify if the integrant is required to be smoothed

* ``DEBG`` (default: ``2`` = warning)
    Debugging verbosity

Vegas integration parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``NCVG`` (default: ``1D5``)
    Number of function calls

* ``NCSG`` (default: ``100``)
    Number of points to probe

* ``ITVG`` (default: ``10``)
    Number of iterations for the integration

Kinematics parameters
---------------------

Incoming state
~~~~~~~~~~~~~~

* ``IN1P|INPP`` (default: ``6500.0``)
    First incoming particle’s momentum, in GeV

* ``IN2P|INPE`` (default: ``6500.0``)
    Second incoming particle’s momentum, in GeV

* ``PMOD`` (default: ``2`` = elastic)
    First incoming particle’s remnant mode or `structure functions modelling </structure-functions>`_

* ``EMOD`` (default: ``2`` = elastic)
    Second incoming particle’s remnant mode or `structure functions modelling </structure-functions>`_

* ``Q2MN``\ &\ ``Q2MX`` (default: ``0.0`` → ``1D5``)
    Q² range for the exchanged parton, in GeV²

Central system
~~~~~~~~~~~~~~

* ``PAIR`` (default: ``13`` = $\\mu$)
    Outgoing leptons’ PDG identifier

* ``PTCT`` (default: ``3.0``)
    Minimal transverse momentum for any single central particle, in  GeV/c

* ``MSCT`` (default: ``0.0``)
    Minimal central system mass, in GeV/c²

* ``ECUT`` (default: ``0.0``)
    Minimal energy for any single central particle, in GeV

* ``ETMN``\ &\ ``ETMX`` (default: ``-2.5`` → ``2.5``)
    Pseudo-rapidity range for central outgoing particles

* ``YMIN``\ &\ ``YMAX`` (default: ``-5.0`` → ``5.0``)
    Rapidity range for central outgoing particles

Outgoing protons / remnants
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``HADR``
    Hadronisation algorithm

* ``MXMN`` & ``MXMX`` (default: ``1.07`` :math:`=m_p+m _ {\pi^{0}}` → ``320.0``)
    Invariant mass range of proton remnants, in GeV/c²

* ``XIMN`` & ``XIMX`` (default: ``0.0`` → ``1.0``)
    .. versionadded:: 0.9.2

    Fractional longitudinal momentum loss :math:`\xi = \Delta p/p`

.. _configuration-card-example:

Configuration card example
--------------------------

The generation of 100k single-dissociative $\\gg{\mu^+\mu^-}$ events at 13 TeV with the `LPAIR matrix element </processes/lpair>`_ implementation with the following phase space cuts:

* :math:`\pt(\mu^\pm)>` 25 GeV, :math:`\lvert\eta(\mu^\pm)\rvert<` 2.5
* 1.07 $< M_X <$ 1000 GeV

can be steered using the following card:

.. code:: fortran

   PROC lpair
   MODE 3      ! inelastic-elastic
   PAIR 13     ! muons
   IN1P 6500.
   IN2P 6500.
   PMOD 11     ! Suri-Yennie
   PTCT 25.
   ETMN -2.5
   ETMX 2.5
   ECUT 0.
   MXMN 1.07
   MXMX 1000.
   NGEN 100000 ! generate 100k events

This configuration is equivalent to the *Python card* shown `here <cards-python#configuration-card-example>`_.
