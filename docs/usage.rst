General usage
=============

The installation guide for the library and general examples may be found `here <install>`_.

Several tools are shipped to ease the user interaction with the library.
All of these can be compiled (after the ``CMake`` environment is set as described `here <install>`_) through a simple ``make <tool_name>``.
The associated **output executable** may then be found in the ``test/`` directory of the build environment (e.g. in ``test/tool_name`` for the example above).

These tools are:

* ``cepgen`` (``test/cepgen.cpp>``) for a simple run computing the **process cross section** and launching an **events generation**.
  According to the steering card content, the events may or may not be stored on disk after production.
  Additionally, the per-event callback function dumps events in the terminal at a frequency defined by the steering card.
  This tool may be used to validate a configuration file and modified to fit the user’s needs.
* ``cepgen-event`` (``test/cepgen-event.cpp>``) to **compute the cross section** and **store events** in one of the supported ASCII output format.
* ``cepgen-root`` (``test/cepgen-root.cxx>``), generating a **ROOT file** with **events** and **run information**.
  This latter may be used in parallel to the two ``ROOT::CepGenRun`` and ``ROOT::CepGenEvent`` helper reader objects for a compact analysis workflow.

CepGen may either be steered through the modification of its **internal parameters definition**, or using the various **steering modules** available.
For the latter, two types of cards are currently supported for the run parameterisation.

* `Structured configurations <cards-python>`_
* `LPAIR-like steering cards <cards-lpair>`_
