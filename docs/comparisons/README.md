# Per-notebook comparison data

One page per notebook, 128 in all, holding every number and every figure that was compared between the reference implementation and OpenEM. Across the set that is **6,003 paired numbers** (5,219 of them compared; the rest are near-zero) and **1,125 figures**.

The notebooks belong to Flexcompute and are **not** redistributed here. Each page links to the notebook on the official example library; what this directory contains is our own measurement of the difference between the two runs.

The summary table with one row per notebook, including the verdict and a one-line note on each, is [`../validation.md`](../validation.md).

## How a comparison is produced

Both sides run the same notebook top to bottom. The reference side is either the output archived with the official notebook, or a rerun against the reference implementation; each page says which. Cells are then matched in order, and within a matched cell:

- **Numbers.** Printed lines carrying numbers are aligned between the two sides. Within an aligned line the numbers are paired left to right, so a line printing three numbers yields three pairs.
- **Figures.** Images produced by the same cell are compared pixel by pixel, and also by eye. The pixel fraction below is only a screen; see the field table for what it does and does not mean.

## What the fields mean

| Field | Meaning |
|---|---|
| `numbers[].cell` | index of the notebook cell the value was printed from |
| `numbers[].label` | the printed line the value came from, truncated to 90 characters |
| `numbers[].reference` / `.openem` | the two values |
| `numbers[].rel_diff` | `abs(reference - openem) / denominator`, as a **fraction**, not a percentage. The denominator is `max(abs(reference), 0.01 * L)` where `L` is the largest absolute value on that same printed line. Using the line maximum keeps a tiny value on a line of large ones from producing a meaningless ratio. |
| `numbers[].near_zero` | true when `abs(reference) < 0.01 * L`, i.e. this value is at least 100 times smaller than the largest number on its line. These are **excluded** from `n_numbers_compared`, `max_rel_diff`, `n_over_5pct` and `n_over_20pct`. |
| `figures[].pixel_diff` | fraction of pixels whose largest per-channel RGB difference exceeds 30 out of 255. A **screening number, not a physical error**: an autoscaled colour bar, antialiasing, or a shifted tick label all register, so a figure can differ by tens of percent by this measure and still show the same physics. Every figure was also compared by eye, and that is what the verdict reflects. |
| `cells.mismatched` | cells whose printed lines could not be aligned one to one |

## Why a scraped number can differ wildly while the notebook still agrees

The numbers below are **scraped automatically** from printed output. That is a screen for finding places worth looking at, not the basis of any verdict. Three things routinely produce a large difference that is not a solver difference:

1. **The pairing lines up values that are not counterparts.** Numbers are paired by position within a matched printed line. When a notebook prints a table row, a dict dump or a bare unlabelled number, the two sides' lines can match while the individual values do not correspond. `Bandstructure` prints a resonance table whose two sides found a different number of resonances, so the columns pair across different rows; `ParameterScan` prints one dict per scan point and pairs values from different points.
2. **The value is not a physical quantity.** `GeneticAlgorithmReflector` prints the index and bit string of each individual in a genetic algorithm. All 97 of its over-5% entries are those listings, which are stochastic by construction and carry no physics.
3. **It is a per-iteration value of a non-convex optimization.** `grad_norm` and `objective` printed each iteration follow the optimization trajectory, and two trajectories separate even when the solver agrees step for step. The notebooks where this happens are judged on the first step: `Autograd5BoundaryGradients` on its initial objective, 0.722 against 0.721, and `Autograd21GaPLightExtractor` on its first-step objective, 9.46 against 9.463.

Concretely: of the 99 notebooks judged **agree**, 13 contain at least one scraped number over 5%, and every one of those traces to a cause above. The conclusion quantity of each, which is what the verdict is set on, is in its row of [`../validation.md`](../validation.md) together with a note on what was found.

## The verdicts

| Verdict | Notebooks | Meaning |
|---|---:|---|
| **agree** | 99 | largest relative difference at most 5% on every comparable number and conclusion quantity, with no real difference in any figure |
| **close** | 11 | between 5% and 20%, or the difference sits on near-zero quantities or on a metric that does not itself converge |
| **trajectory divergence** | 18 | the forward solve and the first-step gradient match, but after a full non-convex optimization the two sides reach different designs; this kind is inherently incomparable |

## Index

### Tutorial (34)

| Notebook | Verdict | Numbers | Largest rel. diff | Figures | Largest pixel diff |
|---|---|---:|---:|---:|---:|
| [AbsorbingBoundaryReflection](AbsorbingBoundaryReflection.md) | close | 0 | - | 5 | 7.52% |
| [AnimationTutorial](AnimationTutorial.md) | agree | 0 | - | 3 | 0.28% |
| [BeerLambert](BeerLambert.md) | agree | 4 | 1.01% | 3 | 6.00% |
| [BoundaryConditions](BoundaryConditions.md) | agree | 0 | - | 14 | 37.87% |
| [CavityFOM](CavityFOM.md) | agree | 8 | 0.04% | 3 | 0.18% |
| [CustomFieldSource](CustomFieldSource.md) | agree | 0 | - | 5 | 8.93% |
| [Dispersion](Dispersion.md) | agree | 0 | - | 7 | 1.29% |
| [FarFieldProjectionNonDecayingField](FarFieldProjectionNonDecayingField.md) | agree | 0 | - | 4 | 34.28% |
| [FieldProjections](FieldProjections.md) | agree | 0 | - | 18 | 18.02% |
| [FullyAnisotropic](FullyAnisotropic.md) | agree | 0 | - | 5 | 11.57% |
| [GDSExport](GDSExport.md) | agree | 0 | - | 3 | 4.73% |
| [GratingEfficiency](GratingEfficiency.md) | agree | 19 | 0.13% | 7 | 8.78% |
| [GroupDelayCalculation](GroupDelayCalculation.md) | agree | 0 | - | 4 | 1.26% |
| [Gyrotropic](Gyrotropic.md) | agree | 0 | - | 3 | 2.05% |
| [MMIMeepBenchmark](MMIMeepBenchmark.md) | agree | 0 | - | 2 | 0.46% |
| [ModalSourcesMonitors](ModalSourcesMonitors.md) | agree | 23 | 93.63% | 11 | 10.23% |
| [ModeSimulation](ModeSimulation.md) | agree | 246 | 0.33% | 20 | 9.75% |
| [ModeSolver](ModeSolver.md) | agree | 12 | 1.13% | 19 | 13.77% |
| [Near2FarSphereRCS](Near2FarSphereRCS.md) | agree | 1 | 33.91% | 3 | 2.55% |
| [PECSphereRCS](PECSphereRCS.md) | agree | 0 | - | 2 | 4.84% |
| [ParameterScan](ParameterScan.md) | agree | 79 | 92.28% | 5 | 6.27% |
| [ParameterScanWebRun](ParameterScanWebRun.md) | agree | 30 | 680% | 5 | 9.32% |
| [PlasmonicNanoparticle](PlasmonicNanoparticle.md) | agree | 0 | - | 6 | 17.65% |
| [Primer](Primer.md) | agree | 90 | 0 | 6 | 0 |
| [ResonanceFinder](ResonanceFinder.md) | agree | 8 | 2.96% | 4 | 41.44% |
| [STLImport](STLImport.md) | agree | 0 | - | 11 | 25.13% |
| [SourceNormalization](SourceNormalization.md) | agree | 3 | 0.75% | 12 | 6.73% |
| [StartHere](StartHere.md) | agree | 4 | 0 | 1 | 4.23% |
| [Symmetry](Symmetry.md) | agree | 0 | - | 23 | 1.02% |
| [TidyFab0GC](TidyFab0GC.md) | trajectory divergence | 157 | 1,091% | 9 | 12.73% |
| [TimeModulationTutorial](TimeModulationTutorial.md) | close | 2 | 15.91% | 3 | 8.72% |
| [VizData](VizData.md) | agree | 11 | 0.21% | 8 | 1.11% |
| [WaveguidePluginDemonstration](WaveguidePluginDemonstration.md) | agree | 32 | 3.92e-05% | 18 | 7.20% |
| [WebAPI](WebAPI.md) | agree | 0 | - | 0 | - |

### Example Library (68)

| Notebook | Verdict | Numbers | Largest rel. diff | Figures | Largest pixel diff |
|---|---|---:|---:|---:|---:|
| [8ChannelDemultiplexer](8ChannelDemultiplexer.md) | agree | 0 | - | 17 | 7.66% |
| [90BendPolarizationSplitterRotator](90BendPolarizationSplitterRotator.md) | agree | 0 | - | 10 | 13.30% |
| [90OpticalHybrid](90OpticalHybrid.md) | agree | 0 | - | 10 | 6.24% |
| [AllDielectricStructuralColor](AllDielectricStructuralColor.md) | close | 114 | 94.16% | 2 | 2.73% |
| [AnisotropicMetamaterialBroadbandPBS](AnisotropicMetamaterialBroadbandPBS.md) | agree | 2 | 0.07% | 8 | 11.21% |
| [Bandstructure](Bandstructure.md) | agree | 2 | 318% | 5 | 1.15% |
| [BilayerSiNSiGC](BilayerSiNSiGC.md) | close | 271 | 180% | 10 | 7.85% |
| [BilevelPSR](BilevelPSR.md) | agree | 0 | - | 10 | 16.40% |
| [BiosensorGrating](BiosensorGrating.md) | agree | 0 | - | 3 | 0.69% |
| [BraggGratings](BraggGratings.md) | agree | 0 | - | 3 | 2.93% |
| [BroadbandDirectionalCoupler](BroadbandDirectionalCoupler.md) | agree | 4 | 1.92e-09% | 13 | 11.67% |
| [BullseyeCavityPSO](BullseyeCavityPSO.md) | agree | 12 | 0 | 7 | 62.26% |
| [CMOSRGBSensor](CMOSRGBSensor.md) | close | 1 | 0 | 6 | 13.55% |
| [DielectricMetasurfaceAbsorber](DielectricMetasurfaceAbsorber.md) | agree | 0 | - | 7 | 0.68% |
| [DirectionalCoupler](DirectionalCoupler.md) | agree | 2 | 0.11% | 5 | 5.26% |
| [DirectionalScatteringNanodisks](DirectionalScatteringNanodisks.md) | agree | 0 | - | 3 | 6.18% |
| [DisorderedPlasmonicColor](DisorderedPlasmonicColor.md) | agree | 1 | 0 | 4 | 10.28% |
| [DistributedBraggReflectorCavity](DistributedBraggReflectorCavity.md) | agree | 1 | 0 | 4 | 12.34% |
| [EdgeCoupler](EdgeCoupler.md) | agree | 0 | - | 8 | 7.76% |
| [EffectiveIndexApproximation](EffectiveIndexApproximation.md) | agree | 12 | 3.32% | 10 | 8.34% |
| [EulerWaveguideBend](EulerWaveguideBend.md) | close | 0 | - | 7 | 5.94% |
| [FreeFormCoupler](FreeFormCoupler.md) | agree | 0 | - | 3 | 7.84% |
| [GeneticAlgorithmReflector](GeneticAlgorithmReflector.md) | agree | 511 | 425% | 5 | 9.87% |
| [GradientMetasurfaceReflector](GradientMetasurfaceReflector.md) | agree | 0 | - | 5 | 7.36% |
| [GrapheneMetamaterial](GrapheneMetamaterial.md) | agree | 0 | - | 5 | 6.05% |
| [GratingCoupler](GratingCoupler.md) | agree | 3 | 0 | 4 | 7.84% |
| [HexagonalLatticeBands](HexagonalLatticeBands.md) | agree | 20 | 0 | 5 | 8.06% |
| [HighQGe](HighQGe.md) | agree | 0 | - | 2 | 6.40% |
| [HighQSi](HighQSi.md) | agree | 0 | - | 2 | 5.64% |
| [MIMResonator](MIMResonator.md) | agree | 0 | - | 3 | 6.20% |
| [MMI1x4](MMI1x4.md) | agree | 0 | - | 9 | 16.42% |
| [MMIPowerSplitter2x2](MMIPowerSplitter2x2.md) | agree | 5 | 0 | 4 | 0.34% |
| [MaxwellStressTensor](MaxwellStressTensor.md) | agree | 5 | 0 | 4 | 4.94% |
| [MicrowaveFrequencySelectiveSurface](MicrowaveFrequencySelectiveSurface.md) | agree | 0 | - | 3 | 8.22% |
| [MoS2Waveguide](MoS2Waveguide.md) | agree | 0 | - | 6 | 1.22% |
| [MultipoleExpansion](MultipoleExpansion.md) | agree | 9 | 99.27% | 5 | 5.78% |
| [NanobeamCavity](NanobeamCavity.md) | close | 17 | 730% | 7 | 10.76% |
| [NanostructuredBoronNitride](NanostructuredBoronNitride.md) | agree | 0 | - | 5 | 0.07% |
| [NonHermitianMetagratings](NonHermitianMetagratings.md) | agree | 1 | 0.03% | 6 | 14.22% |
| [OpticalLuneburgLens](OpticalLuneburgLens.md) | agree | 0 | - | 9 | 5.69% |
| [OpticalSwitchDBS](OpticalSwitchDBS.md) | agree | 25 | 0.12% | 5 | 17.74% |
| [OptimizedL3](OptimizedL3.md) | agree | 1 | 5.74e-04% | 5 | 11.36% |
| [ParticleSwarmOptimizedPBS](ParticleSwarmOptimizedPBS.md) | trajectory divergence | 19 | 33.13% | 6 | 14.89% |
| [PhaseChangeAntennas](PhaseChangeAntennas.md) | agree | 14 | 0.72% | 4 | 25.01% |
| [PhotonicCrystalWaveguidePolarizationFilter](PhotonicCrystalWaveguidePolarizationFilter.md) | agree | 0 | - | 5 | 11.10% |
| [PhotonicSpinSelector](PhotonicSpinSelector.md) | trajectory divergence | 54 | 1,265% | 4 | 34.84% |
| [PlasmonicNanorodArray](PlasmonicNanorodArray.md) | agree | 1 | 0 | 4 | 26.83% |
| [PlasmonicYagiUdaNanoantenna](PlasmonicYagiUdaNanoantenna.md) | agree | 0 | - | 7 | 0.94% |
| [PolarizationSplitterRotator](PolarizationSplitterRotator.md) | agree | 0 | - | 7 | 7.49% |
| [QMRSMetasurface](QMRSMetasurface.md) | agree | 0 | - | 7 | 52.45% |
| [RadiativeCoolingGlass](RadiativeCoolingGlass.md) | agree | 2 | 0 | 7 | 12.37% |
| [RingResonator](RingResonator.md) | agree | 0 | - | 3 | 6.90% |
| [SWGBroadbandPolarizer](SWGBroadbandPolarizer.md) | agree | 0 | - | 7 | 16.57% |
| [SWGWaveguideCrossing](SWGWaveguideCrossing.md) | close | 8 | 6.81% | 3 | 9.35% |
| [SbendCMAES](SbendCMAES.md) | agree | 18 | 1.22% | 5 | 8.69% |
| [ScaleInvariantWaveguide](ScaleInvariantWaveguide.md) | agree | 0 | - | 7 | 23.51% |
| [StripToSlotConverters](StripToSlotConverters.md) | agree | 0 | - | 12 | 0.51% |
| [THzDemultiplexerFilter](THzDemultiplexerFilter.md) | agree | 0 | - | 5 | 9.09% |
| [TopoQuantumPhC](TopoQuantumPhC.md) | agree | 0 | - | 6 | 10.00% |
| [TunableChiralMetasurface](TunableChiralMetasurface.md) | agree | 0 | - | 7 | 9.62% |
| [VerticalGratingCoupler](VerticalGratingCoupler.md) | agree | 0 | - | 4 | 14.51% |
| [VortexMetasurface](VortexMetasurface.md) | agree | 1 | 99.96% | 5 | 9.83% |
| [WaveguideCrossing](WaveguideCrossing.md) | agree | 0 | - | 3 | 4.09% |
| [WaveguideGratingAntenna](WaveguideGratingAntenna.md) | agree | 1 | 0 | 3 | 4.76% |
| [WaveguideSizeConverter](WaveguideSizeConverter.md) | agree | 0 | - | 12 | 7.86% |
| [WaveguideToRingCoupling](WaveguideToRingCoupling.md) | close | 5 | 378% | 6 | 19.73% |
| [YJunction](YJunction.md) | agree | 0 | - | 6 | 5.45% |
| [ZonePlateFieldProjection](ZonePlateFieldProjection.md) | agree | 0 | - | 4 | 12.54% |

### Inverse Design (26)

| Notebook | Verdict | Numbers | Largest rel. diff | Figures | Largest pixel diff |
|---|---|---:|---:|---:|---:|
| [Autograd0Quickstart](Autograd0Quickstart.md) | agree | 21 | 1.59% | 1 | 4.60% |
| [Autograd0QuickstartII](Autograd0QuickstartII.md) | agree | 14 | 0 | 1 | 3.52% |
| [Autograd10YBranchLevelSet](Autograd10YBranchLevelSet.md) | trajectory divergence | 702 | 1,014% | 10 | 22.16% |
| [Autograd12LightExtractor](Autograd12LightExtractor.md) | trajectory divergence | 5 | 94.27% | 11 | 53.01% |
| [Autograd13Metasurface](Autograd13Metasurface.md) | agree | 10 | 374% | 7 | 22.10% |
| [Autograd15Antenna](Autograd15Antenna.md) | trajectory divergence | 128 | 574% | 28 | 8.59% |
| [Autograd16BilayerCoupler](Autograd16BilayerCoupler.md) | trajectory divergence | 96 | 34.10% | 37 | 20.42% |
| [Autograd17BandPassFilter](Autograd17BandPassFilter.md) | trajectory divergence | 149 | 66.38% | 43 | 28.05% |
| [Autograd18TopologyBend](Autograd18TopologyBend.md) | trajectory divergence | 100 | 71.65% | 32 | 10.19% |
| [Autograd1Intro](Autograd1Intro.md) | agree | 7 | 2.93% | 4 | 3.98% |
| [Autograd21GaPLightExtractor](Autograd21GaPLightExtractor.md) | agree | 210 | 41.18% | 73 | 11.43% |
| [Autograd22PhotonicCrystal](Autograd22PhotonicCrystal.md) | trajectory divergence | 31 | 68.12% | 8 | 12.38% |
| [Autograd24DigitalSplitter](Autograd24DigitalSplitter.md) | trajectory divergence | 125 | 31.15% | 19 | 12.98% |
| [Autograd25WaveguideCrossing](Autograd25WaveguideCrossing.md) | close | 79 | 77.05% | 33 | 14.59% |
| [Autograd26DiffractiveBeamSplitter](Autograd26DiffractiveBeamSplitter.md) | trajectory divergence | 11 | 5,179% | 5 | 29.56% |
| [Autograd27Smatrix](Autograd27Smatrix.md) | trajectory divergence | 108 | 1,825% | 33 | 9.59% |
| [Autograd29SourceGradients](Autograd29SourceGradients.md) | trajectory divergence | 91 | 53.68% | 6 | 8.54% |
| [Autograd2GradientChecking](Autograd2GradientChecking.md) | agree | 32 | 0.14% | 1 | 5.00% |
| [Autograd30ParallelAdjoint](Autograd30ParallelAdjoint.md) | agree | 6 | 53.42% | 1 | 4.03% |
| [Autograd31GratingCouplerWithBeamOptimization](Autograd31GratingCouplerWithBeamOptimization.md) | trajectory divergence | 126 | 138% | 7 | 8.92% |
| [Autograd3InverseDesign](Autograd3InverseDesign.md) | trajectory divergence | 118 | 70.54% | 26 | 9.97% |
| [Autograd5BoundaryGradients](Autograd5BoundaryGradients.md) | agree | 252 | 80.51% | 8 | 8.48% |
| [Autograd6GratingCoupler](Autograd6GratingCoupler.md) | trajectory divergence | 377 | 4,019% | 4 | 14.19% |
| [Autograd7Metalens](Autograd7Metalens.md) | close | 186 | 2,363% | 6 | 10.79% |
| [Autograd8WaveguideBend](Autograd8WaveguideBend.md) | agree | 151 | 317% | 8 | 6.18% |
| [Autograd9WDM](Autograd9WDM.md) | trajectory divergence | 213 | 95.55% | 58 | 28.19% |
