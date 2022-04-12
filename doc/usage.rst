General usage
=============

The installation guide for the library and general examples may be found `here <install>`_.

Several tools are shipped to ease the user interaction with the library.
All of these can be compiled (after the ``CMake`` environment is set as described `here <install>`_) through a simple ``make <tool_name>``.
The associated **output executable** may then be found in your ``build/`` directory.
Conventionally, a ``make install`` command can be run, either as administrator to install these system-wise, or locally by overriding the ``CMAKE_INSTALL_PREFIX`` environment variable.

For instance, one may find the following:

* ``cepgen`` (see `here <https://github.com/cepgen/cepgen/blob/master/src/cepgen.cpp>`_) is the main CepGen executable.
  It allows computing the **process cross section** and launching an **events generation**.
  According to the steering card content, the events may or may not be stored on disk after production.
  For instance, in `the LPAIR process Python steering card <https://github.com/cepgen/cepgen/blob/master/Cards/lpair_cfg.py>`_, the event content is printed on screen by the ``dump`` module every 10,000 generations.

  A few extra arguments can be passed through the command line. As for other utilitaries, these are generally documented adding the ``-h``, or ``--help`` flag after the command name.
  Beside this "help" flag a few additional "standard" flags are handled, i.e.

  * ``-d`` or ``--debug`` turns on the verbose mode, and produces *a lot* of useful debugging information
  * ``-v`` or ``--version`` shows the compiled CepGen version tag and returns.

  For the ``cepgen`` executable, one also finds:

  * ``-i`` or ``--config`` allows to provide the input steering card path,
  * ``-n`` or ``--num-events`` lets the user override the number of events to be generated,
  * ``-o`` or ``--output`` allows to define a comma-separated list of output modules names to which events shall be fed (with standard parameters values)

  In addition, two useful flags are:

  * ``-s`` or ``--safe-mode``, to disable the loading in the runtime environment of all add-ons found, either in the ``build`` directory, or in the ``CEPGEN_PATH``.
    Useful in the case where a package under its experimental phase is corrupted.
  * ``-a`` or ``--add-ons``, to specify a comma-separated list of shared libraries to load at the runtime.
    For instance in combination to the earlier, it allows to spot precisely which package is inducing an error, or to load new processes generated and pre-compiled externally.

  Finally, the ``-l`` or ``--list-modules`` flag allows to dump a list of modules (output, event modification, drawing, functional parser, ...) accessible in the runtime environment after the loading of the various shared libraries.

* ``cepgenDescribeModules`` gives a user-readable description of one, or all modules defined in the runtime environment. The ``-a``, or ``--all`` flag produces a complete list of all modules, similar to the one given in `this raw list <raw-modules>`_.

We suggest you to use CepGen through the main executable, in combination with the various **steering modules** available.
These are, to explicitely name a few,

* Python, `structured configurations <cards-python>`_
* `LPAIR-like steering cards <cards-lpair>`_
* JSON, XML, INFO-formatted configurations, as handled by the `Boost.PropertyTree <https://www.boost.org/doc/libs/1_78_0/doc/html/property_tree.html>`_ serialisation library (requires ``CepGenBoost`` to be loaded).

We can however assist you in producing a packaged version, fitting to your event generation scheme.
Please contact us for more details.

