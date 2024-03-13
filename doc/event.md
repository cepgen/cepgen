# Event content definition

## Event content object

A couple of specific objects containing the full information of a physics event at a given snapshot of its processing are defined.
The {cpp:class}`cepgen::Event` object can be seen as an extended collection of {cpp:class}`cepgen::Particle` components with the attributes defined hereafter.

### Particle properties

The {cpp:class}`cepgen::Particle` object holds the minimal amount of information needed to define a particle: a momentum (as a {cpp:class}`cepgen::Momentum` wrapper of the 4-momentum kinematics), some status (as a {cpp:type}`cepgen::Particle::Status` enum) and role in event (the {cpp:type}`cepgen::Particle::Role` enum) flags, a {cpp:class}`cepgen::pdgid_t` PDG identifier, a mothers/daughters parentage, helicity, ...

The detailed description of the object's attributes can be shown below:

````{toggle}
```{doxygenclass} cepgen::Particle
:members:
```
````

### Event properties

In addition to a map between a particles' role in the event and a collection of all particles with this role, the event object {cpp:class}`cepgen::Event` also handles a few useful accessors, helpers, and properties.
To quote a few:

- {cpp:func}`cepgen::Event::cmEnergy` allows to compute the centre-of-mass energy of the incoming particles states,
- {cpp:func}`cepgen::Event::compress` returns a `compressed` event, where all intermediate "vertex-like" states are skipped, only to keep the minimal 2-to-N event content,
- {cpp:func}`cepgen::Event::missingEnergy` allows to compute the missing momentum from a "standard" detector's perspective (the sum of all invisible particles, such as neutrinos/neutralinos/...)

The detailed description of the object's methods and attributes can be displayed below:

````{toggle}
```{doxygenclass} cepgen::Event
:members:
```
````

## Event accessors

The interaction with event content was simplified from CepGen version `1.1.0` to introduce a base placeholder object with the full knowledge of the {cpp:class}`cepgen::Event` object.

````{toggle}
```{doxygenclass} cepgen::EventHandler
:members:
```
````

The implementations of event accessors can be classified into three categories: importers, modifiers, and exporters.

### Event importer

This category was introduced in CepGen [version `1.2.0`](/changelog.md#release-1-2-0), and allows to load {cpp:class}`cepgen::Event` objects from an external source.
It can be useful when, for instance, CepGen is used as an interface to several modification algorithms (as described below), or to convert a dataset between two event formats.

````{toggle}
```{doxygenclass} cepgen::EventImporter
:members:
```
````

### Event modifier

Event modifiers are basic buiding blocks for e.g. hadronisation algorithms interfaces, or treatment algorithms to match one event content definition to another.
They are fully described in the [Event modification algorithms](/event-modifiers.md) section of this documentation.

### Event exporter

Event exporters, or output modules are providing an interface to major output formats definitions and writers.
Again, a full documentation of this sub-class can be found in the [Output formats](/output-formats.md) section.
