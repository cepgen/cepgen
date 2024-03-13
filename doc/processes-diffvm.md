---
orphan: true
---

```{title} Diffractive production of vector mesons 𝛽
```

# Diffractive vector meson production

The $ep$-induced $\gamma^{(\ast)}\Pom\rightarrow$ vector meson production process may be selected from the `diffvm` label.
It was designed for HERA physics cases, hence the asymmetric electron-proton initial beam kinematics.

## Process-specific options

- `vmFlavour`: PDG id of the VM produced

- `vmMode`: vector meson production mode,
  `protonMode`: proton side emission type

  - `BeamMode.GluonFragmentation := -1`
  - `BeamMode.Elastic := 0`
  - `BeamMode.StandardFragmentation := 1`
  - `BeamMode.NucleonPionsDecay := 2`

- `photonMode`: photon emission type

  - `PhotonMode.Fixed := -1`
  - `PhotonMode.InvK := 0`
  - `PhotonMode.WWA := 1`
  - `PhotonMode.ABTSmith = 2`
  - `PhotonMode.AandS := 3`

Furthermore, the following parameters can be steered in this process:

### `slopeParameters`

- `wb0` is the $w _ {\gamma p}$ centre-of-mass energy (in GeV) at which `b0` $=b_0$ was measured
- `amxb0` is the mass $M_X$ of the diffractively dissociating hadronic system $X$ for which $b_0$ was measured
- `anexp` is a power-law exponent

### `pomeronParameters`

- `epsilonW` and `epsilonM` are controlling the intercept of the pomeron trajectory (minus 1).
  The first one steers the rise of $\sigma _ {\gamma p}$ with $W$, while the second controls the $M_X$ spectrum
- `alpha1` and `alpha1m` control the pomerons trajectory’s $\alpha'$ (therefore, expressed in GeV¯²)

### `vmParameters`

- `lambda` and `eprop` control the $Q^2$-dependence of the
  total production cross section, through

```{math}
\sigma(Q^2) = \sigma_0\left(1 + Q^2/\Lambda^2\right)^{-\epsilon _ {\rm prop}}
```

- `xi` and `chi` control the behaviour of the
  longitudinal-to-transverse cross section ratio through

```{math}
\frac{\sigma_L(Q^2)}{\sigma_T(Q^2)}=\frac{\xi Q^2/m^2}{1+\xi\chi Q^2/m^2}.
```

```{note}
In this scheme, $\sigma_L/\sigma_T\to\xi Q^2/m^2$ for low-$Q^2$, and $1/\chi$ for high-$Q^2$.
```

```{doxygenclass} cepgen::proc::DiffVM
```
