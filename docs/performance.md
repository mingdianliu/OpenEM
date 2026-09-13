# 性能：与参考实现的逐例对照

单张 NVIDIA H800，fp32。对照口径是参考实现自己日志里的分段计时（准备 / 步进 / 读出），
我们按同样三段拆。**步进段可比性最强**——两边做的是同一件事：在同一套网格上推同样多的时间步。

## 汇总（118 个 simulation）

**小结**：118 个 simulation 的三段计时。以步进段（两边做同一件事，可比性最强）看，我们比云端快的有 **89 个**、慢的 29 个，比值中位 **0.66**（逐例中位数看我们更快，而按总时长加权的合计比是 1.54，说明慢的都是大例）。慢得最多的是 MetalOxideSunscreen_sim0（30）、BullseyeCavityPSO_sim0（14）、PlasmonicNanoparticle_sim0（3.9）、VortexMetasurface_sim0（3.6），快得最多的是 THzDemultiplexerFilter_sim0（0.06）、ResonanceFinder_sim0（0.11）、MetasurfaceBIC_th8（0.13）、8ChannelDemultiplexer_sim0（0.18）（数值为我们÷云端）。准备段两边构成不同（我们含容器冷启动与内核编译），读出段取决于监视器数量，都不作判定依据。

**总量对比（118 个仿真加总）**：步进段云端 **6,278 s** 对我们 **9,647 s**，总步进比 **1.54×**；端到端云端 **7,424 s** 对我们 **9,934 s**，总端到端比 **1.34×**。分类看（步进 / 端到端）：Tutorial 1.39× / 1.12×；Example Library 1.55× / 1.37×；InverseDesign 0.77× / 0.45×。端到端比小于步进比，是因为准备段云端占 13%、我们只占 2%（我们那 228 s 里还含容器冷启动与内核编译）；我们的时间 97% 花在步进上，云端是 85%。**一句话：整体上我们比云端慢，步进慢 54%、端到端慢 34%；但差距集中在少数大例上，逐例中位数反而是我们快。**

## 成本

**小结**：同样 118 个 simulation，按整卡占用秒数折算人民币。我们更便宜的 **99 个**、云端更便宜的 19 个，成本比中位 **2.0×**。省得最多的是 BeerLambert_sim0（32）、BoundaryConditions_sim0（27）、THzDemultiplexerFilter_sim0（23）、CustomFieldSource_sim_3（16）（倍数=云端÷我们）；云端更划算的集中在大例上，因为 credit 不与墙钟成正比、有下限与档位，而我们按 GPU 秒实打实算。

## 停机步与对称性缩减

**小结**：118 个 simulation 里判定一致 **78** 个（其中 66 个停机步数完全相同），步数有差异 35 个、折叠有差异 3 个。我们停得晚的 20 个（最多的是 MetalOxideSunscreen_sim0（5.4）、CustomFieldSource_sim0（4）、CustomFieldSource_sim_3（4）），停得早的 18 个（THzDemultiplexerFilter_sim0（0.03）、QMRSMetasurface_sim0（0.17）、MMIMeepBenchmark_sim0（0.22））。停得晚多是相量收敛判据比场衰减保守，停得早多是能量衰减先触发；两者都不改频域结果，只影响耗时。

## 怎么读这些数字

**总量比我们慢、逐例中位比我们快**，两个都是真的，不矛盾：总量被几个大例子主导，
而那几个正是我们最吃亏的。所以看单个例子的加速比要看它自己那一行，别用总量推。

参考实现是运行在云端的商业服务，它的「准备」段包含我们没有的环节（任务排队、上传、网格与材料生成——
那部分我们直接取用它客户端的结果，见 [validation.md](validation.md) 开头）。
所以**端到端的比值不代表求解器快慢**，只有步进段是干净的对照。

下表只列 [validation.md](validation.md) 里那 128 本对应的 simulation，按计算量从小到大。


## Tutorial（入门教程）（36 个）

| 例子 | 计算量（格×步） | 参考 合计 | OpenEM 合计 | 时间比 |
|---|---:|---:|---:|---:|
| [SourceNormalization_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/SourceNormalization/) | 3.06e+06 | 0.9 | 0.6 | 1.77 |
| [Gyrotropic_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/Gyrotropic/) | 2.69e+07 | 1.0 | 0.5 | 0.64 |
| [Primer_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/Primer/) | 4.62e+07 | 1.1 | 1.4 | 2.43 |
| [CustomFieldSource_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/CustomFieldSource/) | 2.82e+08 | 1.2 | 1.0 | 0.96 |
| [TimeModulationTutorial_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/TimeModulationTutorial/) | 8.17e+08 | 2.2 | 0.5 | 0.33 |
| [BeerLambert_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BeerLambert/) | 8.45e+08 | 1.4 | 0.2 | 0.20 |
| [WebAPI_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/WebAPI/) | 1.18e+09 | 1.0 | 0.7 | 0.57 |
| [BoundaryConditions_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BoundaryConditions/) | 1.55e+09 | 1.2 | 0.3 | 0.30 |
| [ModeSimulation_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ModeSimulation/) | 1.69e+09 | 4.8 | 0.6 | 0.27 |
| [ModeSolver_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ModeSolver/) | 1.69e+09 | 2.6 | 0.5 | 0.32 |
| [VizData_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/VizData/) | 2.96e+09 | 5.1 | 2.3 | 0.53 |
| [FieldProjections_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/FieldProjections/) | 5.37e+09 | 3.0 | 1.1 | 0.37 |
| [FieldProjections_sim4](https://www.flexcompute.com/tidy3d/examples/notebooks/FieldProjections/) | 5.37e+09 | 3.3 | 1.4 | 0.33 |
| [ModalSourcesMonitors_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ModalSourcesMonitors/) | 7.29e+09 | 4.6 | 0.9 | 0.53 |
| [GratingEfficiency_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/GratingEfficiency/) | 8.84e+09 | 4.2 | 0.8 | 0.21 |
| [StartHere_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/StartHere/) | 9.64e+09 | 3.1 | 0.8 | 0.35 |
| [STLImport_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/STLImport/) | 1.14e+10 | 5.5 | 1.4 | 0.25 |
| [TidyFab0GC_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/TidyFab0GC/) | 1.16e+10 | 6.1 | 4.4 | 0.83 |
| [FieldProjections_sim3](https://www.flexcompute.com/tidy3d/examples/notebooks/FieldProjections/) | 1.28e+10 | 4.3 | 1.6 | 0.44 |
| [Symmetry_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/Symmetry/) | 2.54e+10 | 4.9 | 5.3 | 1.21 |
| [FieldProjections_sim5](https://www.flexcompute.com/tidy3d/examples/notebooks/FieldProjections/) | 2.95e+10 | 6.1 | 3.0 | 0.54 |
| [WaveguidePluginDemonstration_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguidePluginDemonstration/) | 3.81e+10 | 22.4 | 7.6 | 0.89 |
| [ResonanceFinder_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ResonanceFinder/) | 4.28e+10 | 28.8 | 3.1 | 0.11 |
| [AnimationTutorial_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/AnimationTutorial/) | 5.01e+10 | 4.9 | 3.0 | 0.71 |
| [AbsorbingBoundaryReflection_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/AbsorbingBoundaryReflection/) | 6.07e+10 | 8.6 | 1.7 | 0.25 |
| [Dispersion_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/Dispersion/) | 7.70e+10 | 8.9 | 14.7 | 1.68 |
| [MMIMeepBenchmark_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MMIMeepBenchmark/) | 1.96e+11 | 25.7 | 10.7 | 0.57 |
| [CavityFOM_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/CavityFOM/) | 2.50e+11 | 42.7 | 30.7 | 0.74 |
| [FarFieldProjectionNonDecayingField_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/FarFieldProjectionNonDecayingField/) | 3.38e+11 | 25.3 | 22.7 | 1.06 |
| [FullyAnisotropic_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/FullyAnisotropic/) | 3.91e+11 | 37.4 | 39.5 | 0.87 |
| [GroupDelayCalculation_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/GroupDelayCalculation/) | 5.75e+11 | 43.6 | 44.2 | 1.20 |
| [ParameterScanWebRun_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ParameterScanWebRun/) | 7.29e+11 | 55.5 | 43.0 | 0.90 |
| [ParameterScan_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ParameterScan/) | 9.22e+11 | 75.8 | 52.8 | 0.84 |
| [Near2FarSphereRCS_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/Near2FarSphereRCS/) | 1.11e+12 | 77.9 | 67.3 | 0.87 |
| [PECSphereRCS_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PECSphereRCS/) | 1.18e+12 | 43.3 | 57.4 | 1.84 |
| [PlasmonicNanoparticle_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicNanoparticle/) | 1.39e+12 | 122.6 | 447.8 | 3.85 |

## Example Library（案例库）（63 个）

| 例子 | 计算量（格×步） | 参考 合计 | OpenEM 合计 | 时间比 |
|---|---:|---:|---:|---:|
| [DistributedBraggReflectorCavity_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/DistributedBraggReflectorCavity/) | 3.86e+05 | 0.9 | 0.9 | 1.49 |
| [HighQGe_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/HighQGe/) | 8.85e+07 | 0.8 | 0.5 | 1.07 |
| [HighQSi_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/HighQSi/) | 1.14e+08 | 1.2 | 0.8 | 0.66 |
| [EffectiveIndexApproximation_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/EffectiveIndexApproximation/) | 1.39e+08 | 1.4 | 0.9 | 0.96 |
| [GradientMetasurfaceReflector_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/GradientMetasurfaceReflector/) | 1.80e+08 | 1.4 | 1.0 | 0.72 |
| [BilayerSiNSiGC_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BilayerSiNSiGC/) | 2.38e+08 | 6.2 | 0.7 | 0.97 |
| [DirectionalScatteringNanodisks_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/DirectionalScatteringNanodisks/) | 5.40e+08 | 1.7 | 2.9 | 0.35 |
| [MicrowaveFrequencySelectiveSurface_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MicrowaveFrequencySelectiveSurface/) | 1.12e+09 | 1.9 | 2.9 | 2.62 |
| [GrapheneMetamaterial_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/GrapheneMetamaterial/) | 1.57e+09 | 1.4 | 1.1 | 0.87 |
| [VerticalGratingCoupler_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/VerticalGratingCoupler/) | 1.81e+09 | 7.2 | 1.7 | 0.95 |
| [GratingCoupler_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/GratingCoupler/) | 2.18e+09 | 7.4 | 1.6 | 0.73 |
| [BiosensorGrating_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BiosensorGrating/) | 3.37e+09 | 5.1 | 5.0 | 1.17 |
| [GeneticAlgorithmReflector_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/GeneticAlgorithmReflector/) | 6.49e+09 | 4.2 | 2.4 | 0.46 |
| [HexagonalLatticeBands_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/HexagonalLatticeBands/) | 6.58e+09 | 10.1 | 2.3 | 0.19 |
| [AllDielectricStructuralColor_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/AllDielectricStructuralColor/) | 7.21e+09 | 4.3 | 0.9 | 0.22 |
| [MaxwellStressTensor_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MaxwellStressTensor/) | 7.93e+09 | 6.2 | 1.1 | 0.20 |
| [ParticleSwarmOptimizedPBS_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ParticleSwarmOptimizedPBS/) | 8.70e+09 | 6.5 | 1.0 | 0.26 |
| [OpticalSwitchDBS_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/OpticalSwitchDBS/) | 1.13e+10 | 5.0 | 2.7 | 0.53 |
| [DisorderedPlasmonicColor_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/DisorderedPlasmonicColor/) | 1.32e+10 | 14.4 | 10.8 | 0.73 |
| [BroadbandDirectionalCoupler_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BroadbandDirectionalCoupler/) | 1.48e+10 | 18.3 | 5.6 | 0.44 |
| [SbendCMAES_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/SbendCMAES/) | 1.86e+10 | 9.3 | 3.8 | 0.26 |
| [PhotonicSpinSelector_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PhotonicSpinSelector/) | 3.32e+10 | 13.6 | 3.0 | 0.33 |
| [MIMResonator_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MIMResonator/) | 3.62e+10 | 7.3 | 7.9 | 1.28 |
| [StripToSlotConverters_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/StripToSlotConverters/) | 4.07e+10 | 11.2 | 3.0 | 0.62 |
| [EulerWaveguideBend_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/EulerWaveguideBend/) | 4.71e+10 | 22.7 | 5.0 | 0.31 |
| [NonHermitianMetagratings_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/NonHermitianMetagratings/) | 5.89e+10 | 9.6 | 8.2 | 0.99 |
| [MMIPowerSplitter2x2_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MMIPowerSplitter2x2/) | 7.06e+10 | 18.3 | 6.2 | 0.59 |
| [WaveguideGratingAntenna_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideGratingAntenna/) | 7.55e+10 | 22.3 | 8.3 | 0.25 |
| [OpticalLuneburgLens_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/OpticalLuneburgLens/) | 7.85e+10 | 7.3 | 2.0 | 0.36 |
| [PlasmonicYagiUdaNanoantenna_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicYagiUdaNanoantenna/) | 8.06e+10 | 13.8 | 7.2 | 0.56 |
| [DirectionalCoupler_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/DirectionalCoupler/) | 8.08e+10 | 18.4 | 7.5 | 0.66 |
| [YJunction_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/YJunction/) | 9.91e+10 | 38.2 | 22.1 | 0.65 |
| [MultipoleExpansion_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MultipoleExpansion/) | 1.10e+11 | 56.3 | 27.0 | 0.45 |
| [TunableChiralMetasurface_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/TunableChiralMetasurface/) | 1.40e+11 | 325.0 | 86.2 | 0.25 |
| [PhaseChangeAntennas_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PhaseChangeAntennas/) | 1.46e+11 | 23.9 | 8.7 | 0.34 |
| [EdgeCoupler_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/EdgeCoupler/) | 1.70e+11 | 21.7 | 5.0 | 0.36 |
| [TopoQuantumPhC_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/TopoQuantumPhC/) | 1.96e+11 | 29.0 | 42.2 | 1.64 |
| [DielectricMetasurfaceAbsorber_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/DielectricMetasurfaceAbsorber/) | 2.30e+11 | 67.9 | 24.2 | 0.36 |
| [WaveguideSizeConverter_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideSizeConverter/) | 2.44e+11 | 20.5 | 11.6 | 0.63 |
| [90OpticalHybrid_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/90OpticalHybrid/) | 3.28e+11 | 40.9 | 22.0 | 0.87 |
| [NanobeamCavity_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/NanobeamCavity/) | 3.35e+11 | 94.5 | 34.5 | 0.41 |
| [SWGBroadbandPolarizer_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/SWGBroadbandPolarizer/) | 3.58e+11 | 84.6 | 77.9 | 1.48 |
| [AnisotropicMetamaterialBroadbandPBS_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/AnisotropicMetamaterialBroadbandPBS/) | 3.64e+11 | 84.5 | 44.1 | 0.63 |
| [8ChannelDemultiplexer_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/8ChannelDemultiplexer/) | 3.72e+11 | 61.2 | 10.4 | 0.18 |
| [WaveguideToRingCoupling_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideToRingCoupling/) | 5.70e+11 | 131.0 | 19.5 | 0.18 |
| [WaveguideCrossing_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideCrossing/) | 7.15e+11 | 69.3 | 41.5 | 0.72 |
| [QMRSMetasurface_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/QMRSMetasurface/) | 7.18e+11 | 222.2 | 53.0 | 0.23 |
| [BullseyeCavityPSO_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BullseyeCavityPSO/) | 8.92e+11 | 79.8 | 247.2 | 13.56 |
| [PolarizationSplitterRotator_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PolarizationSplitterRotator/) | 1.10e+12 | 112.1 | 76.9 | 0.89 |
| [FreeFormCoupler_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/FreeFormCoupler/) | 1.23e+12 | 159.7 | 33.7 | 0.36 |
| [MMI1x4_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MMI1x4/) | 1.36e+12 | 38.7 | 8.8 | 0.32 |
| [ZonePlateFieldProjection_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/ZonePlateFieldProjection/) | 1.39e+12 | 160.0 | 70.2 | 0.70 |
| [PhotonicCrystalWaveguidePolarizationFilter_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/PhotonicCrystalWaveguidePolarizationFilter/) | 1.40e+12 | 69.4 | 78.4 | 1.17 |
| [SWGWaveguideCrossing_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/SWGWaveguideCrossing/) | 1.57e+12 | 181.4 | 367.0 | 2.33 |
| [90BendPolarizationSplitterRotator_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/90BendPolarizationSplitterRotator/) | 1.66e+12 | 120.5 | 233.9 | 2.18 |
| [OptimizedL3_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/OptimizedL3/) | 1.79e+12 | 124.2 | 108.4 | 0.89 |
| [MoS2Waveguide_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/MoS2Waveguide/) | 1.97e+12 | 292.2 | 182.7 | 1.02 |
| [BraggGratings_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BraggGratings/) | 2.55e+12 | 154.9 | 199.6 | 1.35 |
| [VortexMetasurface_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/VortexMetasurface/) | 4.11e+12 | 100.4 | 289.1 | 3.61 |
| [NanostructuredBoronNitride_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/NanostructuredBoronNitride/) | 4.60e+12 | 234.3 | 381.4 | 1.65 |
| [RingResonator_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/RingResonator/) | 5.47e+12 | 472.6 | 819.4 | 1.94 |
| [THzDemultiplexerFilter_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/THzDemultiplexerFilter/) | 1.73e+13 | 1069.6 | 66.0 | 0.06 |
| [BilevelPSR_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/BilevelPSR/) | 4.14e+13 | 1234.6 | 2508.0 | 2.11 |

## Inverse Design（逆设计）（1 个）

| 例子 | 计算量（格×步） | 参考 合计 | OpenEM 合计 | 时间比 |
|---|---:|---:|---:|---:|
| [Autograd26DiffractiveBeamSplitter_sim0](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd26DiffractiveBeamSplitter/) | 1.60e+09 | 2.4 | 0.7 | 0.22 |

## 优化开关

速度相关的开关（内核融合、内域瘦内核、批量清算等）都在 `openem/knobs.py` 里登记，
默认值是在 H800 上实测定标的。**所有开关都不改数值**：切换后 22 个金标准仿真的输出必须按位相同。
换硬件可能需要重新定标，但那只影响速度。

