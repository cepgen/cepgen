(processes-devel)=

# Processes development

We recommend the implementation of all new processes in C++, as described below.
However, we conveniently provide an interface to Fortran implementations for historical reasons.
A general description of this latter can be seen down this page.

In general, two possible tracks are given to the CepGen process developer: the implementation of the full matrix element, starting from the beams kinematics definition down to the event weight and event content, giving a full flexibility to the whole process definition, or using a parton-from-beam-factorised simplified implementation.

The earlier relies on a sub-classing of the {cpp:class}`cepgen::Proc::Process` object, already providing some useful kinematics variables of interest:

````{toggle}
```{doxygenclass} cepgen::proc::Process
:members:
:protected-members:
```
````

As noted above, to avoid code redundancy, the following helper accessors can be used to retrieve and set the various event content kinematics at the level of the process:

- {cpp:func}`cepgen::proc::Process::pA`, and {cpp:func}`cepgen::proc::Process::pX`, for the positive-\$z\$ incoming/outgoing beam particle's kinematics ;
- {cpp:func}`cepgen::proc::Process::pB`, and {cpp:func}`cepgen::proc::Process::pY`, for the negative-\$z\$ incoming/outgoing beam particle's kinematics ;
- {cpp:func}`cepgen::proc::Process::q1`, and {cpp:func}`cepgen::proc::Process::q2`, for the positive, and negative-\$z\$ incoming parton kinematics ;
- {cpp:func}`cepgen::proc::Process::pc`, with an index (starting at 0) accessor for the central system particles' 4-momentum.

However, in most cases, such an access to the full beam dynamics and kinematics content is not necessarily required.
Therefore, the so-called factorised processes definition can be useful for most of user cases.

## Factorised processes

```{versionadded} 0.9
```

CepGen attempts to give the user a relative freedom in its implementation of factorised, two-parton level processes.
The kinematics of the latter can indeed be generated collinearly to the incoming beams kinematics (as done in major modern central exclusive processes generator), or with a physical transverse momentum, according to the approach described in [the reference papers](/zz-bibliography.md).
The $\kt$-factorisation approach allows a direct factorisation of any hard process (e.g. photon- or gluon-induced productions) while accounting for transverse components of parton virtualities.
For instance, a $pp\to p^{(\ast)}(\ggx)p^{(\ast)}$ matrix element can be factorised through the following formalism:

```{math}
\mathrm d\sigma = \int \frac{\mathrm d^2{\mathbf q_{\mathrm T}^2}_1}{\pi {\mathbf q_{\mathrm T}^2}_1}
                    {\cal F}_{\gamma/p}^{\rm el/inel}(x_1,{\mathbf q_{\mathrm T}^2}_1)
                    \int \frac{\mathrm d^2{\mathbf q_{\mathrm T}^2}_2}{\pi {\mathbf q_{\mathrm T}^2}_2}
                    {\cal F}_{\gamma/p}^{\rm el/inel}(x_2,{\mathbf q_{\mathrm T}^2}_2) ~ \mathrm d\sigma^\ast,
```

where $\mathcal F_{\gamma/p}^{\rm el/inel}(x_i,\vecqt_i)$ are unintegrated parton densities, and $\mathrm d\sigma^\ast$ the hard process factorised out of the total matrix element.

### Unintegrated photon fluxes

Elastic unintegrated photon densities are expressed as functions of the proton electric and magnetic form factors $G_E$ and $G_M$:

```{math}
\mathcal F_{\gamma/p}^{\rm el}(\xi,\vecqt^2) = \frac{\alpha}{\pi}\left[(1-\xi)\left(\frac{\vecqt^2}{\vecqt^2+\xi^2 m_p^2}\right)^2 F_E(Q^2)+\frac{\xi^2}{4}\left(\frac{\vecqt^2}{\vecqt^2+\xi^2 m_p^2}\right) F_M(Q^2)\right].
```

The inelastic contribution further requires both the diffractive state four-momentum norm $M_X$ and a [proton structure functions parameterisation](/structure-functions.md) as an input:

```{math}
\mathcal F_{\gamma/p}^{\rm inel}(\xi,\vecqt^2) = \frac{\alpha}{\pi}\Bigg[(1-\xi)\left(\frac{\vecqt^2}{\vecqt^2+\xi(M_X^2-m_p^2)+\xi^2 m_p^2}\right)^2\frac{F_2(\xbj,Q^2)}{Q^2+M_X^2-m_p^2}+{}\\
  {}+\frac{\xi^2}{4}\frac{1}{\xbj^2} \left(\frac{\vecqt^2}{\vecqt^2+\xi(M_X^2-m_p^2)+\xi^2 m_p^2}\right) \frac{2\xbj F_1(\xbj,Q^2)}{Q^2+M_X^2-m_p^2}\Bigg],
```

with $\xbj = {Q^2}/({Q^2+M_X^2-m_p^2})$ the Bjorken scaling variable.

### C++ interface

The {class}`cepgen::proc::FactorisedProcess` helper derivated-class of the earlier is introduced to allow the parton emission part to be transparent to the process developper.

```{doxygenclass} cepgen::proc::FactorisedProcess
:outline:
```

````{toggle}
```{doxygenclass} cepgen::proc::FactorisedProcess
:members:
```
````

As this object is pure virtual, the two following members have to be defined for any factorised process:

- {cpp:func}`cepgen::proc::FactorisedProcess::prepareFactorisedPhaseSpace`, to initialise all the process-local kinematic variables (if needed) ;
- {cpp:func}`cepgen::proc::FactorisedProcess::computeFactorisedMatrixElement` to retrieve the central matrix element weight that will be convoluted to the incoming parton fluxes and the event Jacobian to form the total matrix element weight.

The two-parton (or more generally the full central kinematics) can be generated automatically by CepGen according to the two scenarios enumerated above: collinear, or transverse momentum-dependent parton emission.

An experimental interfacing to other phase space mappers (such as Rambo) was recently made technically possible, and is currently under development.

The basic object of interest for this mapping is:

````{toggle}
```{doxygenclass} cepgen::PhaseSpaceGenerator
:members:
:protected-members:
```
````

Which currently contains one instantiation, for 2-to-4 processes (such as `pptoff` and `pptoww`):

````{toggle}
```{doxygenclass} cepgen::PhaseSpaceGenerator2to4
:protected-members:
```
````

Again, the same accessors as for {cpp:class}`cepgen::proc::Process` can be used for the retrieval of the event kinematics used in computing the event weight.

### Fortran interface

For development or testing purposes (and increased flexibility in general), Fortran definitions of physics processes can be fed and used by CepGen for cross section computation and event generation.
In this page a summary and hands-on example of such a Fortran implementation and linkage is described.

Before other things, you will have to provide the definition of your process and its topology for CepGen to handle it properly.
This requires you to fill in a set of predefined common blocks shared between your process definition and the core CepGen instance.

All these new processes are then linked to CepGen through the construction and association of each matrix element definition and interface to the {cpp:class}`cepgen::proc::FortranFactorisedProcess` object.

````{toggle}
```{doxygenclass} cepgen::proc::FortranFactorisedProcess
:members:
```
````

In a later paragraph, this linking recipe will be described.

#### Output event kinematics

This common block lets CepGen access the whole event structure and kinematics information.
Beside the outgoing beam particles 4-momenta `px` and `py`, it contains the PDG id and 4-momentum of all central system particles in `pc`.

For this latter, a maximum multiplicity of 10 particles is handled by default in CepGen.

```fortran
common/evtkin/nout,ipdg,idum2,pc,px,py
integer nout,idum2,ipdg(10)
double precision pc(4,10),px(4),py(4)
```

Should this maximal value not be sufficient for your implementation purposes, please get in touch with us.

#### Process definition

This block specifies all useful run-level parameters, such as:

- `icontri`, the kinematics mode considered (`1` = elastic, `2-3` = single-dissociative, `4` = double-dissociative),
- `iflux1`, `iflux2`, the type of $\kt$-dependent flux used to parameterise the incoming partons’ kinematics,
- `sfmod`, the structure functions modelling used [^f1],
- `inp1`, `inp2`, the longitudinal beam momenta (in GeV/c).

Additionally, for heavy ions initial states, the following parameters are used:

- `a_nuc1`, `a_nuc2`, the mass numbers of positive- and negative-z incoming beam particles,
- `z_nuc1`, `z_nuc2`, the atomic numbers of positive- and negative-z incoming beam particles,

For increased modularity, a set of process parameters directly parsed from the input cards (see e.g. [the processParameters block](/cards-python.md#process-processparameters-parameters-block)) may also be introduced.
At the Fortran process definition level, all parameters may be accessed through a set of interfacing functions:

```{doxygenfunction} cepgen_param_int_
```

```{doxygenfunction} cepgen_param_real_
```

````{note}
These latter may be translated into the following Fortran signature:

```fortran
external CepGen_param_int,CepGen_param_real
integer CepGen_param_int
double precision CepGen_param_real
```
````

Those two can be called from Fortran using e.g.

```fortran
integer ival
real*8 rval
! [...]
ival = CepGen_param_int('int_parameter', 0)      ! first argument is parameter name
rval = CepGen_param_real('real_parameter', 0.d0) ! second argument is default value
```

#### $k_{\rm T}$-factorisation kinematics

This common block holds the following $k_{\rm T}$-factorisation processes kinematics variables (generated and filled by the core CepGen instance):

- `q1t` and `q2t`, the 2-norm of transverse incoming parton virtualities,
- `phiq1t` and `phiq2t`, the equivalent azimutal angle $\phi$ of these virtualities,
- `y1` and `y2` the outgoing system particles’ rapidities,
- `ptdiff` and `phiptdiff` the 2-norm and azimutal angle of the outgoing system’s transverse momentum difference,
- `am_x` and `am_y` the dissociative outgoing beams’ invariant masses (or $m_p$ for elastic parton emission from proton).

```fortran
   common/ktkin/q1t,q2t,phiq1t,phiq2t,y1,y2,ptdiff,phiptdiff,
&     am_x,am_y
   double precision q1t,q2t,phiq1t,phiq2t,y1,y2,ptdiff,phiptdiff,
&     am_x,am_y
```

#### Phase space cuts

A limited set of pre-registered phase space cuts (one boolean switch, a lower and an upper value to apply) may be found in this interfacing library, such as:

- single particles $\pt$, energy, pseudo-rapidity $\eta$,
- central system invariant mass, $\pt$,
- central system correlations: $\Delta y$.

```fortran
   common/constants/am_p,units,pi
   common/kincuts/ipt,iene,ieta,iinvm,iptsum,idely,
&     pt_min,pt_max,ene_min,ene_max,eta_min,eta_max,
&     invm_min,invm_max,ptsum_min,ptsum_max,
&     dely_min,dely_max
   logical ipt,iene,ieta,iinvm,iptsum,idely
   double precision pt_min,pt_max,ene_min,ene_max,eta_min,eta_max,
&     invm_min,invm_max,ptsum_min,ptsum_max,
&     dely_min,dely_max
```

Should this collection not be sufficient for your purposes, please contact us.

#### Physics constants

This block introduces all useful and standard physical constants in double precision to help the process definition: \$m_p\$, GeV\${}^2\$/barn, or \$pi\$.

```fortran
common/constants/am_p,units,pi
double precision am_p,units,pi
```

#### Matrix element definition

The core definition of your process is to be implemented within a new `.f` file you will store in the `External/Processes/` directory of your local CepGen install.

Currently, the only following Fortran process is defined: `nucl_to_ff`

For illustration in this documentation page, we will be using a dummy process name, such as `dummy_f77_process`.
For the sake of simplicity in your processes definition bookkeeping we suggest you to match the filename to the matrix element weighting function name.
In our example, this corresponds to `External/Processes/dummy_f77_process.f`

We expect your weighting function to follow the C-equivalent signature:

```c
double dummy_f77_process_(void);
```

This might get translated, for instance, into the following minimal working example (here using the free form Fortran coding convention):

```{literalinclude} ../CepGenProcesses/Examples/dummy_f77_process.f
:language: fortran
:linenos: true
```

With the kinematics common blocks defined through the include statement described above.

#### Helper tools

To ease the work of the process developer, we provide utilitaries for the evaluation of unintegrated parton fluxes and other physics quantities through the CepGen core C++ library methods.
The external methods exposed to the Fortran process through the include statement above are:

```{doxygenfunction} cepgen_kt_flux_
```

```{doxygenfunction} cepgen_kt_flux_hi_
```

```{doxygenfunction}} cepgen_particle_charge_
```

```{doxygenfunction} cepgen_particle_mass_
```

These can be translated in the following Fortran subroutines/functions definitions:

```fortran
external CepGen_kT_flux,CepGen_kT_flux_HI
double precision CepGen_kT_flux,CepGen_kT_flux_HI
external CepGen_particle_charge,CepGen_particle_mass
double precision CepGen_particle_charge,CepGen_particle_mass
```

#### Overall linking

To interface your process to CepGen, edit the `CepGenProcesses/ProcessesWrapper.cpp` file to link your function implementation to the core processes module.

The following macros are used to declare the function name and link it to CepGen.

```{doxygendefine} DECLARE_FORTRAN_FUNCTION
```

```{doxygendefine} REGISTER_FORTRAN_PROCESS
```

```{note}
No single, nor double trailing underscore (\_) is required for this name registration.
```

For the latter, the signature is the following:

```C
REGISTER_FORTRAN_PROCESS(my_dummy_process, dummy_f77_process,
                         "this process is only for illustration purposes")
```

Or, following this nomenclature:

- the first argument corresponds to the process name (i.e. may be used with the my_dummy_process string in steering cards),
- the second argument is the function name as defined above,
- the third one is a one-line, human-readable process description that will be shown at the run start.

For this version, the following process is already registered:

```{literalinclude} ../CepGenProcesses/ProcessesWrapper.cpp
:caption: List of registered external Fortran processes, as imported from ``CepGenProcesses/ProcessesWrapper.cpp``.
:language: cpp
```

[^f1]: See [this page](/structure-functions.md) for a complete list of integer-type structure functions definitions.
