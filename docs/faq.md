# Frequently Asked Questions

- [Are excitations peak or RMS?](#are-excitations-peak-or-rms)
- [Does the peak/RMS choice affect the capacitance, inductance, or resistance matrices?](#does-the-peakrms-choice-affect-the-capacitance-inductance-or-resistance-matrices)
- [Why is my axisymmetric loss off by a factor of 2*pi*r?](#why-is-my-axisymmetric-loss-off-by-a-factor-of-2pir)
- [Why does the solver reject `simulation.frequency`?](#why-does-the-solver-reject-simulationfrequency)
- [Why must every domain attribute be claimed by a region?](#why-must-every-domain-attribute-be-claimed-by-a-region)
- [Do the reported losses include hysteresis?](#do-the-reported-losses-include-hysteresis)
- [Why is a Robin boundary condition rejected?](#why-is-a-robin-boundary-condition-rejected)

---

## Are excitations peak or RMS?

**Peak (amplitude).** In time-harmonic (magnetoquasistatic) runs, every
`excitations[].value` is interpreted as the amplitude of the phasor, not its
RMS value.

This matters because the time-averaged Joule loss is computed as

```
P = (1/2) * sigma * |E|^2
```

where the factor of one half is precisely the time average of a product of
*peak-amplitude* phasors. Supplying RMS phasors applies that averaging twice
and **halves every reported loss**.

The hazard is that a nameplate current is conventionally an RMS quantity, so
transcribing one directly into a config file is the natural mistake. There is
deliberately no rms/peak selector in the schema, and:

> The convention is **documented but not enforced.** No code validates or
> converts the excitation. An RMS value is silently accepted and produces a
> plausible, quietly wrong answer.

If your source data is RMS, multiply by `sqrt(2)` before putting it in the
config file.

The value is carried unscaled from `excitations[].value` through the
right-hand side into the solved potentials and port voltages, so the convention
you choose propagates to every downstream quantity. This matches the convention
used by common commercial AC/DC and eddy-current tools.

## Does the peak/RMS choice affect the capacitance, inductance, or resistance matrices?

**No.** Extracted coupling matrices are *ratios* of solved quantities to the
drive that produced them -- resistance is `Re(V/I)`, capacitance is `Q/V`, and
inductance is `flux/I`. A common scale factor on the excitation cancels between
numerator and denominator, so the matrices are convention-free.

The asymmetry is worth internalizing:

| Output | Sensitive to peak vs. RMS? | Why |
|--------|----------------------------|-----|
| Capacitance matrix | No | Ratio `Q/V` |
| Inductance matrix | No | Ratio `flux/I` |
| Resistance matrix | No | Ratio `Re(V/I)` |
| Loss density / region losses | **Yes** | Absolute quadratic quantity |
| Exported field values | **Yes** | Scale linearly with the drive |

In short: anything linear-and-ratiometric is safe; anything absolute and
quadratic is not. A `coupling_matrix` run is immune. A `field` run that reports
losses is not.

## Why is my axisymmetric loss off by a factor of 2*pi*r?

Loss *density* is a per-unit-volume quantity. Turning it into a total loss
requires the volume measure, and in an axisymmetric model that measure is
`2*pi*r`, supplied by the integrator rather than by the coefficient.

Integrating the loss density with a planar integrator on an axisymmetric
problem omits that measure and yields a **per-radian** loss, not a total. Note
that this `2*pi*r` volume measure is entirely distinct from the physical `1/r`
in the azimuthal drive field `E_phi = V/(2*pi*r)`, which the coefficient does
handle internally.

## Why does the solver reject `simulation.frequency`?

Frequency is a property of a *scenario*, not of a run, so that one file can
sweep frequency without duplicating the rest of the model. Putting it in
`simulation` would make the sweep forms impossible to express.

Move it onto each scenario. See
[`frequency`](config_reference.md#frequency-mqs-only) for the single, list, and
sweep forms.

## Why must every domain attribute be claimed by a region?

Because the failure mode of *not* checking is silent. An unclaimed domain
attribute would receive default vacuum properties and still solve, producing a
plausible but wrong answer. Requiring exactly one region per domain attribute
converts that into an error you see before the solve.

Doubly-claimed attributes are rejected for the same reason: the winner would
otherwise depend on ordering.

## Do the reported losses include hysteresis?

**No.** The loss density is resistive (eddy-current) dissipation only, assuming
linear, isotropic, non-hysteretic conduction. Hysteresis and excess (anomalous)
losses are not represented, so total loss in laminated steel is
**underestimated**.

The result is also single-frequency: it is meaningless for a transient solution
or a superposition of several frequencies.

## Why is a Robin boundary condition rejected?

`robin` is parsed and reserved in the schema but not implemented in any solver.
Rather than silently ignoring it -- which would apply a homogeneous Neumann
condition and quietly change the physics -- the validator rejects it outright.

For open-boundary problems, truncate the domain further out with a Dirichlet
condition and verify that the result is insensitive to the truncation radius.
