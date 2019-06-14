:orphan:

===========================
Installation for developers
===========================

Active development branches
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* ``devel``
  Main branch where most of the new code is being pushed. The Pythia 8 hadronisation part handles by default the collinear emission of a valence quark from which the two-photon system arises.

* ``diffvm``
  Version including the `DiffVM </processes/diffvm>`_ diffractive vector meson production process.

* ``herwig-fragmentation``
  Implementation of the Herwig (aka Herwig 7) cluster fragmentation algorithm.

General recipe
--------------

To obtain a list of mandatory and optional dependencies, have a look at the `general installation procedure <install>`_.
Note that prior to any clone/... steps you need to:

* `create an account on HepForge <https://www.hepforge.org/register>`_,
* contact us to link your account to the list of registered developers, and
* register your public key on `HepForge’s phabricator <https://phab.hepforge.org/>`_ instance:

   <https://phab.hepforge.org/settings/user/[YOUR HEPFORGE USERNAME]/page/ssh/>

.. code:: sh

   git clone ssh://vcs@phab.hepforge.org/source/cepgen.git
   cd cepgen-dev
   mkdir build
   cd build

Then check a development branch out (see the list above for a detailed description of major branch features):

.. code:: sh

   git checkout <branch name>

Then the building procedure can be launched:

.. code:: sh

   cmake ..
   make
