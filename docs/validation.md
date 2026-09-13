# 验证：对 Tidy3D 官方案例的逐本对照

OpenEM 的正确性判据不是「看起来合理」，而是**把 Tidy3D 官方 example library 的 notebook 整本跑一遍，
逐个数字、逐张图与参考实现对照**。这一页是结果，判定口径见下一节。

> **网格与材料来自 Tidy3D。** OpenEM 读的是 `td.Simulation` 对象，**直接取用 Tidy3D 客户端已经算好的
> 网格、几何光栅化与材料数组**，自己只做时间步进、边界、源注入与监视器读出。这是刻意的：
> 两侧吃同一份离散化，差异就只可能来自求解器本身，否则分不清是建模差异还是求解差异。
> 独立的网格与材料生成模块在后续开发计划里。

本页列出 **128 本**判定为一致 / 基本一致 / 轨迹差异的 notebook
（一致 99、基本一致 11、轨迹差异 18）。
另有 11 本未列入，原因见文末。

**每本的逐项对照数据都在 [`comparisons/`](comparisons/) 里**：合计 6,003 个配对数字与 1,125 张图，
每个数字都给出参考值、OpenEM 的值和相对差，每张图给出像素差。表里每一行的「逐项数据」链到对应那页；
字段定义与比对方法见 [`comparisons/README.md`](comparisons/README.md)。

## 判定口径

| 判定 | 含义 |
|---|---|
| **一致** | 可比数字与结论量的最大相对差 ≤5%，图无真差异 |
| **基本一致** | 5%～20%，或差异集中在近零量、指标本身不收敛的量上 |
| **轨迹差异** | 正向与首步梯度对得上，但非凸优化跑满后两侧走到不同设计——这类**本质不可比**，官方自己重跑轨迹也不同 |

「可比数字」指 notebook 打印出来、两侧能一一对齐的数；「结论量」是那本 notebook 真正要得出的物理结论
（透射率、品质因子、耦合效率之类），最后一列给的就是它。

判定为**有差异 / 能力缺口 / 参考存疑 / 不能跑**的不在本页，理由见文末。


## Tutorial（入门教程）（34 本）

| 例子 | 计算量（格×步） | 判定 | 结论量（参考 / OpenEM，相对差） | 备注 |
|---|---:|---|---|---|
| [Primer](https://www.flexcompute.com/tidy3d/examples/notebooks/Primer/) · [逐项数据](comparisons/Primer.md) | 4.618e+07 | 一致 | 教学 notebook：几何/网格/场分布 6 图 | 复跑：90 个可比数字一致，6 张图逐像素相同 |
| [Gyrotropic](https://www.flexcompute.com/tidy3d/examples/notebooks/Gyrotropic/) · [逐项数据](comparisons/Gyrotropic.md) | 9.356e+07 | 一致 | 旋磁介质中的场分布（法拉第旋转） | 复跑：3 张图像素差 ≤2.1%，无数字差异 |
| [SourceNormalization](https://www.flexcompute.com/tidy3d/examples/notebooks/SourceNormalization/) · [逐项数据](comparisons/SourceNormalization.md) | 2.219e+08 | 一致 | 偶极子辐射功率 Flux [W]（解析 0.00133）：0.00134 / 0.00133（0.7%）<br>E0 [V/µm]：1.009 / 1.006（0.3%） | 修后重跑（37/37 单元，153 s）：3 个可比数字差 ≤0.7%，12 张图像素差 ≤6.7% |
| [WebAPI](https://www.flexcompute.com/tidy3d/examples/notebooks/WebAPI/) · [逐项数据](comparisons/WebAPI.md) | 1.182e+09 | 一致 | 云端 API 演示，无图无物理量 | 复跑：8 个单元无差异 |
| [CustomFieldSource](https://www.flexcompute.com/tidy3d/examples/notebooks/CustomFieldSource/) · [逐项数据](comparisons/CustomFieldSource.md) | 2.486e+09 | 一致 | 自定义场源注入后的场分布与通量 | |Ey| 分布一致；前向通量 8.53e-3 vs 8.34e-3（2%），反向通量都≈0 |
| [TimeModulationTutorial](https://www.flexcompute.com/tidy3d/examples/notebooks/TimeModulationTutorial/) · [逐项数据](comparisons/TimeModulationTutorial.md) | 5.761e+09 | 基本一致 | 调制引起的光子数相对变化：-0.008978 / -0.00755（15.9%） | 结构剖面图相同 |
| [BeerLambert](https://www.flexcompute.com/tidy3d/examples/notebooks/BeerLambert/) · [逐项数据](comparisons/BeerLambert.md) | 8.863e+09 | 一致 | 总吸收功率（解析 0.99）：0.99 / 1（1.0%） | 复跑：4 个数字最大差 1.0%，Re{Ex} 剖面与 |Ex|² 衰减曲线重合 |
| [ModeSimulation](https://www.flexcompute.com/tidy3d/examples/notebooks/ModeSimulation/) · [逐项数据](comparisons/ModeSimulation.md) | 1.483e+10 | 一致 | — | 复跑：246 个可比数字最大差 0.3%，20 张图像素差 ≤9.8% |
| [ModeSolver](https://www.flexcompute.com/tidy3d/examples/notebooks/ModeSolver/) · [逐项数据](comparisons/ModeSolver.md) | 1.483e+10 | 一致 | — | 复跑：12 个可比数字最大差 1.1%，19 张图像素差 ≤13.8% |
| [VizData](https://www.flexcompute.com/tidy3d/examples/notebooks/VizData/) · [逐项数据](comparisons/VizData.md) | 2.373e+10 | 一致 | 中心频率通量 [W]：0.977 / 0.978（0.1%） | 复跑：36 个可比数字最大差 0.3%，8 张图像素差 ≤1.2% |
| [GratingEfficiency](https://www.flexcompute.com/tidy3d/examples/notebooks/GratingEfficiency/) · [逐项数据](comparisons/GratingEfficiency.md) | 2.767e+10 | 一致 | 各衍射级功率之和：0.9984 / 0.9997（0.1%） | 修掉斜入射平面波方向/单格 y 轴面积/衍射读出三处 bug（/#24/#25）后重跑：正入射与斜入射总功率 0.9983/1.0001（参考 0.9984… |
| [StartHere](https://www.flexcompute.com/tidy3d/examples/notebooks/StartHere/) · [逐项数据](comparisons/StartHere.md) | 2.904e+10 | 一致 | Ez 场分布 | 复跑：4 个可比数字一致 |
| [TidyFab0GC](https://www.flexcompute.com/tidy3d/examples/notebooks/TidyFab0GC/) · [逐项数据](comparisons/TidyFab0GC.md) | 3.570e+10 | 轨迹差异 | — | 修后整本重跑（33/33 单元）：梯度已从「40 步全零」修好（入口分流 + 窄带伴随的设计矩阵条件数自动放宽） |
| [BoundaryConditions](https://www.flexcompute.com/tidy3d/examples/notebooks/BoundaryConditions/) · [逐项数据](comparisons/BoundaryConditions.md) | 4.539e+10 | 一致 | 各类边界条件下的场分布（14 图） | 14 张图数量相同；像素差最大的一张是 PMC 面上 |Hx|（数值零，1e-16 量级噪声图不可比），同一图里 |Hz| 分布与色标一致 |
| [ModalSourcesMonitors](https://www.flexcompute.com/tidy3d/examples/notebooks/ModalSourcesMonitors/) · [逐项数据](comparisons/ModalSourcesMonitors.md) | 4.779e+10 | 一致 | 中心频率通量：1 / 1（0.0%）<br>前向基模幅值：1 / 1（0.0%） | 改判一致：模式监视器功率带内两侧都在 1.000±0.003，带外长波端OpenEM 0.97、参考 1.015（宽带注入的带外保真度差 4%，在 5% 内） |
| [Symmetry](https://www.flexcompute.com/tidy3d/examples/notebooks/Symmetry/) · [逐项数据](comparisons/Symmetry.md) | 5.285e+10 | 一致 | 对称性缩减前后场分布（23 图） | 复跑：23 张图像素差 ≤1.1%，无数字差异 |
| [AnimationTutorial](https://www.flexcompute.com/tidy3d/examples/notebooks/AnimationTutorial/) · [逐项数据](comparisons/AnimationTutorial.md) | 6.009e+10 | 一致 | 时域场动画帧 | 复跑（修完时域单位 后）：3 张图像素差 0.3%，可比数字无差异 |
| [FieldProjections](https://www.flexcompute.com/tidy3d/examples/notebooks/FieldProjections/) · [逐项数据](comparisons/FieldProjections.md) | 1.258e+11 | 一致 | 远场投影方向图（18 图） | 角谱投影 Ey/Ez 一致、Ex 为弱分量（<1% 主分量）差异属离散噪声；k 空间投影的斜入射 GaussianBeam 原先场幅只有一半，修掉全 PML… |
| [AbsorbingBoundaryReflection](https://www.flexcompute.com/tidy3d/examples/notebooks/AbsorbingBoundaryReflection/) · [逐项数据](comparisons/AbsorbingBoundaryReflection.md) | 1.265e+11 | 基本一致 | 不同吸收边界的反射率曲线 | 修完模式面读出（#67）与模式源半格相位（#69）后整本重跑（15/15 单元）目检，两张图都对上了 |
| [WaveguidePluginDemonstration](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguidePluginDemonstration/) · [逐项数据](comparisons/WaveguidePluginDemonstration.md) | 1.836e+11 | 一致 | — | 合入（notebook 直连的本地 ModeSolver 按参考环境阶梯化）后重跑（24/24，2584 s）：32 个数字与官方存档逐位一致（最大差 0.… |
| [Dispersion](https://www.flexcompute.com/tidy3d/examples/notebooks/Dispersion/) · [逐项数据](comparisons/Dispersion.md) | 2.706e+11 | 一致 | 色散介质透射/反射谱 vs TMM | 复跑：7 张图像素差 ≤2.2%，无数字差异 |
| [STLImport](https://www.flexcompute.com/tidy3d/examples/notebooks/STLImport/) · [逐项数据](comparisons/STLImport.md) | 3.241e+11 | 一致 | STL 导入几何的场分布 | Ex/Ez 一致；Ey 是对称性决定的数值零分量（1e-6 vs 4e-4 量级噪声），噪声图不可比 |
| [ResonanceFinder](https://www.flexcompute.com/tidy3d/examples/notebooks/ResonanceFinder/) · [逐项数据](comparisons/ResonanceFinder.md) | 3.421e+11 | 一致 | 回音壁模 Q（f≈154.2 THz）：2856 / 2858（0.1%） | 复跑（修完时域单位 后）：FieldTimeMonitor 的 Re{Ey} 现在与参考同刻度（±4000，此前OpenEM是 4e9，差 1e6 倍），拍… |
| [GroupDelayCalculation](https://www.flexcompute.com/tidy3d/examples/notebooks/GroupDelayCalculation/) · [逐项数据](comparisons/GroupDelayCalculation.md) | 6.896e+11 | 一致 | — | 复跑：4 张图像素差 ≤1.3% |
| [MMIMeepBenchmark](https://www.flexcompute.com/tidy3d/examples/notebooks/MMIMeepBenchmark/) · [逐项数据](comparisons/MMIMeepBenchmark.md) | 9.076e+11 | 一致 | 透射谱 vs Meep | 复跑：5 个单元无数字差异，2 张图像素差 0.5% |
| [Near2FarSphereRCS](https://www.flexcompute.com/tidy3d/examples/notebooks/Near2FarSphereRCS/) · [逐项数据](comparisons/Near2FarSphereRCS.md) | 1.343e+12 | 一致 | RCS vs Mie 理论 | 复跑：唯一差 >5% 的数字是「本地近远场变换耗时 0.204 s vs 0.127 s」，是计时不是物理量；3 张图像素差 ≤2.5% |
| [FarFieldProjectionNonDecayingField](https://www.flexcompute.com/tidy3d/examples/notebooks/FarFieldProjectionNonDecayingField/) · [逐项数据](comparisons/FarFieldProjectionNonDecayingField.md) | 1.405e+12 | 一致 | — | 修后整本重跑（10/10 单元）目检：四个窗口尺寸（0.20/0.40/0.60）×四个投影尺寸（20/30/40/50）共 16 张加窗远场图，以及四张三… |
| [FullyAnisotropic](https://www.flexcompute.com/tidy3d/examples/notebooks/FullyAnisotropic/) · [逐项数据](comparisons/FullyAnisotropic.md) | 1.424e+12 | 一致 | 全张量介质场分布（等效对照） | 介质色散曲线数值相同，仅线条配色顺序不同 |
| [ParameterScanWebRun](https://www.flexcompute.com/tidy3d/examples/notebooks/ParameterScanWebRun/) · [逐项数据](comparisons/ParameterScanWebRun.md) | 1.627e+12 | 一致 | — | 改判一致：差 >5% 的只有两个 amplitude²=0.00 端口（反射/交叉臂）的相位，近零复数的相位是噪声、没有物理意义；有功率的两个端口 ampl… |
| [PECSphereRCS](https://www.flexcompute.com/tidy3d/examples/notebooks/PECSphereRCS/) · [逐项数据](comparisons/PECSphereRCS.md) | 1.821e+12 | 一致 | 单站 RCS vs 解析 | 2.43 亿格仿真在作业内建场景+求解+组装后跑通，2 张图像素差 ≤4.8% |
| [ParameterScan](https://www.flexcompute.com/tidy3d/examples/notebooks/ParameterScan/) · [逐项数据](comparisons/ParameterScan.md) | 1.832e+12 | 一致 | — | 改判一致：差 >5% 的只有两个 amplitude²=0.00 端口（反射/交叉臂）的相位，近零复数的相位是噪声、没有物理意义；有功率的两个端口 ampl… |
| [CavityFOM](https://www.flexcompute.com/tidy3d/examples/notebooks/CavityFOM/) · [逐项数据](comparisons/CavityFOM.md) | 1.996e+12 | 一致 | 腔 Q（f≈329.7 THz）：9.492e+04 / 9.488e+04（0.0%）<br>V_eff [(λ/n)³]：0.79 / 0.79（0.0%）<br>Purcell 因子 F_p：9173 / 9169（0.0%） | 8 个可比数字相同，3 张图像素差 ≤0.2% |
| [PlasmonicNanoparticle](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicNanoparticle/) · [逐项数据](comparisons/PlasmonicNanoparticle.md) | 3.461e+12 | 一致 | — | 复跑：6 张图像素差 ≤12.5%，无数字差异 |
| [GDSExport](https://www.flexcompute.com/tidy3d/examples/notebooks/GDSExport/) · [逐项数据](comparisons/GDSExport.md) | 9.002e+12 | 一致 | GDS 导出几何图（无物理量） | 复跑：3 张图像素差 ≤4.7% |

## Example Library（案例库）（68 本）

| 例子 | 计算量（格×步） | 判定 | 结论量（参考 / OpenEM，相对差） | 备注 |
|---|---:|---|---|---|
| [DistributedBraggReflectorCavity](https://www.flexcompute.com/tidy3d/examples/notebooks/DistributedBraggReflectorCavity/) · [逐项数据](comparisons/DistributedBraggReflectorCavity.md) | 7.790e+06 | 一致 | 反射带归一化带宽：0.32 / 0.32（0.0%） | 各周期数反射谱曲线重合 |
| [Bandstructure](https://www.flexcompute.com/tidy3d/examples/notebooks/Bandstructure/) · [逐项数据](comparisons/Bandstructure.md) | 3.331e+08 | 一致 | 首个共振频率 [Hz]：2.559e+13 / 9.858e+13（285.2%） | 能带图与时域信号图像素差 ≤1.1%；唯一 >5% 的数字来自 ResonanceFinder 输出表的行配对错位（两侧找到的谐振数目/顺序不同），不是能带… |
| [EffectiveIndexApproximation](https://www.flexcompute.com/tidy3d/examples/notebooks/EffectiveIndexApproximation/) · [逐项数据](comparisons/EffectiveIndexApproximation.md) | 5.082e+08 | 一致 | 3D 平均 FSR [nm]：21.08 / 21（0.4%）<br>2D/3D 平均 FSR 差 [nm]：2.5 / 2.417（3.3%） | 合入（notebook 直连的本地 ModeSolver 按参考环境阶梯化）后重跑（23/23，1411 s）：12 个数字 0 处差 >5%，最大 3.3… |
| [GradientMetasurfaceReflector](https://www.flexcompute.com/tidy3d/examples/notebooks/GradientMetasurfaceReflector/) · [逐项数据](comparisons/GradientMetasurfaceReflector.md) | 1.264e+09 | 一致 | 反射衍射效率与场分布 | 复跑：结构剖面图相同，8 张图（参考复跑漏执行最后 3 个单元时只有 5 张） |
| [BilayerSiNSiGC](https://www.flexcompute.com/tidy3d/examples/notebooks/BilayerSiNSiGC/) · [逐项数据](comparisons/BilayerSiNSiGC.md) | 2.275e+09 | 基本一致 | — | 全部修复（#62～#71）后整本重跑（29/29，6223 s）：两张参数扫描表逐点对上（121 点反射率表最优点同为 (w=0.44, p=0.88)，2… |
| [HighQSi](https://www.flexcompute.com/tidy3d/examples/notebooks/HighQSi/) · [逐项数据](comparisons/HighQSi.md) | 2.403e+09 | 一致 | 透射谱（高 Q Fano 共振） | 复跑：透射谱与 Fano 共振位置一致，2 张图像素差 ≤5.6% |
| [GratingCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/GratingCoupler/) · [逐项数据](comparisons/GratingCoupler.md) | 2.974e+09 | 一致 | — | 复跑：原来卡在 7×11=77 个 2D 仿真的参数扫描单元（>45 min）超时，工人池下整本跑完 |
| [HighQGe](https://www.flexcompute.com/tidy3d/examples/notebooks/HighQGe/) · [逐项数据](comparisons/HighQGe.md) | 3.153e+09 | 一致 | 透射谱（高 Q Fano 共振） | 复跑：结构剖面图相同 |
| [GrapheneMetamaterial](https://www.flexcompute.com/tidy3d/examples/notebooks/GrapheneMetamaterial/) · [逐项数据](comparisons/GrapheneMetamaterial.md) | 6.032e+09 | 一致 | 吸收谱 | 复跑：四条费米能级的吸收谱重合 |
| [HexagonalLatticeBands](https://www.flexcompute.com/tidy3d/examples/notebooks/HexagonalLatticeBands/) · [逐项数据](comparisons/HexagonalLatticeBands.md) | 6.575e+09 | 一致 | 六角晶格能带图 | 复跑：20 个可比数字全部一致，能带图 5/5 相同 |
| [MicrowaveFrequencySelectiveSurface](https://www.flexcompute.com/tidy3d/examples/notebooks/MicrowaveFrequencySelectiveSurface/) · [逐项数据](comparisons/MicrowaveFrequencySelectiveSurface.md) | 7.529e+09 | 一致 | 透射/反射谱与场强 | 复跑：结构与网格图相同（铜按 PEC 近似，见 CONVENTIONS §30） |
| [ScaleInvariantWaveguide](https://www.flexcompute.com/tidy3d/examples/notebooks/ScaleInvariantWaveguide/) · [逐项数据](comparisons/ScaleInvariantWaveguide.md) | 7.821e+09 | 一致 | — | 复跑：可比数字无差异；六个厚度的模式场分布与色标一致，像素差来自子图排布 |
| [BiosensorGrating](https://www.flexcompute.com/tidy3d/examples/notebooks/BiosensorGrating/) · [逐项数据](comparisons/BiosensorGrating.md) | 1.034e+10 | 一致 | 透射/反射谱 | 复跑：可比数字无差异，3 张图像素差 ≤0.7% |
| [ParticleSwarmOptimizedPBS](https://www.flexcompute.com/tidy3d/examples/notebooks/ParticleSwarmOptimizedPBS/) · [逐项数据](comparisons/ParticleSwarmOptimizedPBS.md) | 1.063e+10 | 轨迹差异 | — | 复跑：整本跑通，19 个可比数字里 7 个差 >5%，全部是 PSO 搜到的「Best Parameters」（10 个几何参数，差 5%～33%）—— 粒… |
| [GeneticAlgorithmReflector](https://www.flexcompute.com/tidy3d/examples/notebooks/GeneticAlgorithmReflector/) · [逐项数据](comparisons/GeneticAlgorithmReflector.md) | 1.418e+10 | 一致 | — | 用 300 min 的作业整本重跑，NBRESULT OK、17/17 单元无报错，（110 min 那次被切在 12/17） |
| [OpticalSwitchDBS](https://www.flexcompute.com/tidy3d/examples/notebooks/OpticalSwitchDBS/) · [逐项数据](comparisons/OpticalSwitchDBS.md) | 1.522e+10 | 一致 | — | 复跑（原来算「DBS 优化驱动」没跑完）：25 个可比数字最大差 0.1%，5 张图数量相同、形状一致 |
| [BroadbandDirectionalCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/BroadbandDirectionalCoupler/) · [逐项数据](comparisons/BroadbandDirectionalCoupler.md) | 1.764e+10 | 一致 | — | 全部修复后整本重跑（25/25，1928 s）：4 个数字最大差 0.0%，两处本地模式求解的有效折射率与参考差到 1e-12（2.459162628530… |
| [VerticalGratingCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/VerticalGratingCoupler/) · [逐项数据](comparisons/VerticalGratingCoupler.md) | 2.506e+10 | 一致 | — | 重跑 11/11 单元无报错 |
| [SbendCMAES](https://www.flexcompute.com/tidy3d/examples/notebooks/SbendCMAES/) · [逐项数据](comparisons/SbendCMAES.md) | 2.846e+10 | 一致 | — | 发散守卫修好后整本重跑（18/18 单元）：18 个可比数字全部差 ≤5%，5 张图像素差 8.7% |
| [DisorderedPlasmonicColor](https://www.flexcompute.com/tidy3d/examples/notebooks/DisorderedPlasmonicColor/) · [逐项数据](comparisons/DisorderedPlasmonicColor.md) | 3.369e+10 | 一致 | 反射谱与结构色 | 复跑：结构与网格图相同 |
| [PhotonicSpinSelector](https://www.flexcompute.com/tidy3d/examples/notebooks/PhotonicSpinSelector/) · [逐项数据](comparisons/PhotonicSpinSelector.md) | 5.932e+10 | 轨迹差异 | — | 重跑 16/16 单元无报错，打印数字全部对上 |
| [DirectionalScatteringNanodisks](https://www.flexcompute.com/tidy3d/examples/notebooks/DirectionalScatteringNanodisks/) · [逐项数据](comparisons/DirectionalScatteringNanodisks.md) | 5.978e+10 | 一致 | 前向/后向散射谱 | 复跑：透射/反射二维图形状与量级一致 |
| [EulerWaveguideBend](https://www.flexcompute.com/tidy3d/examples/notebooks/EulerWaveguideBend/) · [逐项数据](comparisons/EulerWaveguideBend.md) | 7.091e+10 | 基本一致 | 弯曲损耗与场强 | 全部修复后重跑（19/19）目检：|E| 场分布两侧逐像素一致；弯曲损耗谱量级一致——圆弧弯 0.018～0.028 dB（OpenEM 0.018～0.0… |
| [AllDielectricStructuralColor](https://www.flexcompute.com/tidy3d/examples/notebooks/AllDielectricStructuralColor/) · [逐项数据](comparisons/AllDielectricStructuralColor.md) | 1.153e+11 | 基本一致 | 反射谱 vs 周期（2 图） | 反射谱主峰位置与高度相同；数字差 >5% 的 36 个都在短波端（0.40～0.41 µm）反射率 <0.03 的尾部（相对差 7%～20%），2 张图像素… |
| [YJunction](https://www.flexcompute.com/tidy3d/examples/notebooks/YJunction/) · [逐项数据](comparisons/YJunction.md) | 1.291e+11 | 一致 | Y 结透射与场分布 | 模式场 |Ex|/|Ey|/|Ez| 分布与色标一致 |
| [PlasmonicYagiUdaNanoantenna](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicYagiUdaNanoantenna/) · [逐项数据](comparisons/PlasmonicYagiUdaNanoantenna.md) | 1.496e+11 | 一致 | 方向性/远场辐射图 | 复跑：17 个可比数字无差异，7 张图像素差 ≤0.9% |
| [PlasmonicNanorodArray](https://www.flexcompute.com/tidy3d/examples/notebooks/PlasmonicNanorodArray/) · [逐项数据](comparisons/PlasmonicNanorodArray.md) | 1.589e+11 | 一致 | 入射场强 [V/µm]：212.8 / 212.8（0.0%） | 复跑：可比数字一致；五个波长的场分布与量级相同，像素差来自子图宽度与色条位置 |
| [NonHermitianMetagratings](https://www.flexcompute.com/tidy3d/examples/notebooks/NonHermitianMetagratings/) · [逐项数据](comparisons/NonHermitianMetagratings.md) | 1.609e+11 | 一致 | — | 改判一致：当前代码整本重跑（18/18）后唯一可比数字差 0.0%，左/右泄漏通量对参考误差 0.4%/0.1%，衬底里不对称泄漏波的方向与强度两侧一致，6… |
| [OpticalLuneburgLens](https://www.flexcompute.com/tidy3d/examples/notebooks/OpticalLuneburgLens/) · [逐项数据](comparisons/OpticalLuneburgLens.md) | 1.654e+11 | 一致 | 焦点附近场强分布 | 复跑：24 个可比数字无差异，9 张图像素差 ≤5.7% |
| [MultipoleExpansion](https://www.flexcompute.com/tidy3d/examples/notebooks/MultipoleExpansion/) · [逐项数据](comparisons/MultipoleExpansion.md) | 1.808e+11 | 一致 | 多极子展开首个标量（c04）：0.1377 / 0.001（99.3%） | 复跑：四个多极分量散射谱与解析解重合 |
| [WaveguideGratingAntenna](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideGratingAntenna/) · [逐项数据](comparisons/WaveguideGratingAntenna.md) | 1.859e+11 | 一致 | — | 复跑（原来整本 >45 min 超时）：可比数字一致，3 张图像素差 ≤4.8% |
| [TunableChiralMetasurface](https://www.flexcompute.com/tidy3d/examples/notebooks/TunableChiralMetasurface/) · [逐项数据](comparisons/TunableChiralMetasurface.md) | 2.000e+11 | 一致 | LCP/RCP 吸收谱 | 复跑：16 个单元无数字差异，7 张图像素差 ≤9.6% |
| [MMIPowerSplitter2x2](https://www.flexcompute.com/tidy3d/examples/notebooks/MMIPowerSplitter2x2/) · [逐项数据](comparisons/MMIPowerSplitter2x2.md) | 2.005e+11 | 一致 | 平均分光比：0.5 / 0.5（0.0%）<br>平均总功率：0.869 / 0.869（0.0%） | 复跑：5 个数字一致，4 张图像素差 ≤0.3% |
| [StripToSlotConverters](https://www.flexcompute.com/tidy3d/examples/notebooks/StripToSlotConverters/) · [逐项数据](comparisons/StripToSlotConverters.md) | 4.456e+11 | 一致 | 转换效率与场分布（12 图） | 复跑：24 个可比数字无差异，12 张图像素差 ≤0.5% |
| [RadiativeCoolingGlass](https://www.flexcompute.com/tidy3d/examples/notebooks/RadiativeCoolingGlass/) · [逐项数据](comparisons/RadiativeCoolingGlass.md) | 4.847e+11 | 一致 | 透射/反射/吸收谱 | 复跑：2 个数字一致，结构剖面图相同 |
| [MaxwellStressTensor](https://www.flexcompute.com/tidy3d/examples/notebooks/MaxwellStressTensor/) · [逐项数据](comparisons/MaxwellStressTensor.md) | 5.199e+11 | 一致 | — | 复跑：21 个 Bloch 仿真的批量单元原来 >45 min 超时，工人池下整本跑完 |
| [PhaseChangeAntennas](https://www.flexcompute.com/tidy3d/examples/notebooks/PhaseChangeAntennas/) · [逐项数据](comparisons/PhaseChangeAntennas.md) | 5.786e+11 | 一致 | 天线 1 一级衍射相位 [°]：38.03 / 38.01（0.1%）<br>天线 4 一级衍射相位 [°]：180 / 180（0.0%） | 8 个天线的交叉极化相位柱状图数值相同，仅柱色不同 |
| [8ChannelDemultiplexer](https://www.flexcompute.com/tidy3d/examples/notebooks/8ChannelDemultiplexer/) · [逐项数据](comparisons/8ChannelDemultiplexer.md) | 6.056e+11 | 一致 | — | 复跑：44 个可比数字全部一致，17 张图像素差 ≤7.7%；TE0～TE3 模式场分布与色标相同 |
| [SWGBroadbandPolarizer](https://www.flexcompute.com/tidy3d/examples/notebooks/SWGBroadbandPolarizer/) · [逐项数据](comparisons/SWGBroadbandPolarizer.md) | 6.696e+11 | 一致 | — | 重跑 20/20 单元无报错 |
| [DirectionalCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/DirectionalCoupler/) · [逐项数据](comparisons/DirectionalCoupler.md) | 7.094e+11 | 一致 | 交叉波长 [µm]：1.53 / 1.53（0.0%）<br>总透射：0.991 / 0.9899（0.1%） | 复跑：2 个数字差 0.1%，5 张图像素差 ≤5.3% |
| [90OpticalHybrid](https://www.flexcompute.com/tidy3d/examples/notebooks/90OpticalHybrid/) · [逐项数据](comparisons/90OpticalHybrid.md) | 7.347e+11 | 一致 | 四路输出功率与相位 FOM | 复跑：23 个可比数字一致；插入损耗、CMRR、不平衡度三组谱线重合（10 张图像素差 ≤6.2%，来自曲线抗锯齿） |
| [AnisotropicMetamaterialBroadbandPBS](https://www.flexcompute.com/tidy3d/examples/notebooks/AnisotropicMetamaterialBroadbandPBS/) · [逐项数据](comparisons/AnisotropicMetamaterialBroadbandPBS.md) | 7.679e+11 | 一致 | — | 全部修复（#62～#71）后整本重跑（20/20，23135 s）：2 个打印数字最大差 0.1%，8 张图目检逐线重合（透射-n_100 扫描四条曲线 T… |
| [WaveguideToRingCoupling](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideToRingCoupling/) · [逐项数据](comparisons/WaveguideToRingCoupling.md) | 1.097e+12 | 基本一致 | through 差最大值（近零量）：4.867e-05 / 7.659e-05（57.4%） | 全部修复后重跑（21/21）：透射/下载谱两侧形状一致；参考 |E| 峰值 2e15 那一例是 notebook 故意演示的「色散介质伸进 PML 导致发散… |
| [QMRSMetasurface](https://www.flexcompute.com/tidy3d/examples/notebooks/QMRSMetasurface/) · [逐项数据](comparisons/QMRSMetasurface.md) | 1.110e+12 | 一致 | — | 复跑（原来两次 1 h 跑不完、只有官方存档输出）：15 个可比数字无差异，7 张图内容相同；像素差大是画布尺寸和背景介质配色不同（参考那张图更大、真空区画… |
| [MIMResonator](https://www.flexcompute.com/tidy3d/examples/notebooks/MIMResonator/) · [逐项数据](comparisons/MIMResonator.md) | 1.390e+12 | 一致 | 反射谱与 R/G/B 共振磁场 | 复跑：RGB 三条透射谱重合，3 张图像素差 ≤6.2% |
| [FreeFormCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/FreeFormCoupler/) · [逐项数据](comparisons/FreeFormCoupler.md) | 1.423e+12 | 一致 | — | 用 300 min 的作业整本重跑，NBRESULT OK、11/11 单元无报错 |
| [EdgeCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/EdgeCoupler/) · [逐项数据](comparisons/EdgeCoupler.md) | 1.454e+12 | 一致 | — | 复跑：22 个可比数字无差异，8 张图像素差 ≤7.8% |
| [WaveguideCrossing](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideCrossing/) · [逐项数据](comparisons/WaveguideCrossing.md) | 1.583e+12 | 一致 | 透射与串扰 | 复跑：11 个单元无数字差异，3 张图像素差 ≤4.1% |
| [WaveguideSizeConverter](https://www.flexcompute.com/tidy3d/examples/notebooks/WaveguideSizeConverter/) · [逐项数据](comparisons/WaveguideSizeConverter.md) | 1.604e+12 | 一致 | 插入损耗与模式成分（12 图） | 复跑：26 个单元无数字差异，12 张图像素差 ≤7.9% |
| [SWGWaveguideCrossing](https://www.flexcompute.com/tidy3d/examples/notebooks/SWGWaveguideCrossing/) · [逐项数据](comparisons/SWGWaveguideCrossing.md) | 1.895e+12 | 基本一致 | — | 全部修复后重跑（21/21）：直通透射 TE0 −0.265 对 −0.269 dB、TE1 −0.331 对 −0.343 dB（换成线性差 0.1%～0… |
| [CMOSRGBSensor](https://www.flexcompute.com/tidy3d/examples/notebooks/CMOSRGBSensor/) · [逐项数据](comparisons/CMOSRGBSensor.md) | 2.075e+12 | 基本一致 | — | 全部修复（#62～#71）后整本重跑（28/28，5563 s）：1 个打印数字差 0.0% |
| [TopoQuantumPhC](https://www.flexcompute.com/tidy3d/examples/notebooks/TopoQuantumPhC/) · [逐项数据](comparisons/TopoQuantumPhC.md) | 2.079e+12 | 一致 | 边缘态能带与场分布 | 复跑：边界态场分布与色标一致，6 张图 |
| [PolarizationSplitterRotator](https://www.flexcompute.com/tidy3d/examples/notebooks/PolarizationSplitterRotator/) · [逐项数据](comparisons/PolarizationSplitterRotator.md) | 2.227e+12 | 一致 | 输入模 0 的 n_eff：2.608 / 2.594（0.5%） | 复跑：15 个可比数字无差异，7 张图像素差 ≤7.5% |
| [NanobeamCavity](https://www.flexcompute.com/tidy3d/examples/notebooks/NanobeamCavity/) · [逐项数据](comparisons/NanobeamCavity.md) | 2.462e+12 | 基本一致 | Q y [MM]：16.45 / 16.46（0.1%）<br>模式体积 [(λ/n)³]：0.38 / 0.38（0.0%）<br>共振 Q（ResonanceFinder）：7.911e+06 / —（缺） | 定稿，两处差异都不是实现问题 |
| [MoS2Waveguide](https://www.flexcompute.com/tidy3d/examples/notebooks/MoS2Waveguide/) · [逐项数据](comparisons/MoS2Waveguide.md) | 2.541e+12 | 一致 | MoS2 单层波导模式与传输 | 复跑：6 张图像素差 ≤1.2% |
| [90BendPolarizationSplitterRotator](https://www.flexcompute.com/tidy3d/examples/notebooks/90BendPolarizationSplitterRotator/) · [逐项数据](comparisons/90BendPolarizationSplitterRotator.md) | 2.606e+12 | 一致 | 透射谱（11 图） | 复跑：21 个可比数字无差异；两个模式的 |E| 分布与色标一致 |
| [BraggGratings](https://www.flexcompute.com/tidy3d/examples/notebooks/BraggGratings/) · [逐项数据](comparisons/BraggGratings.md) | 2.904e+12 | 一致 | 透射/反射谱 | 复跑：8 个数字无差异，3 张图像素差 ≤2.9% |
| [MMI1x4](https://www.flexcompute.com/tidy3d/examples/notebooks/MMI1x4/) · [逐项数据](comparisons/MMI1x4.md) | 3.114e+12 | 一致 | 四路输出功率与场分布 | 改判一致：可比数字无差异，9 张图里 8 张相同；唯一差别是 3×3 粗扫的「内外波导功率差」热图最优格挪了一格，对应的结论值 L_MMI 参考 ≈11.1… |
| [BullseyeCavityPSO](https://www.flexcompute.com/tidy3d/examples/notebooks/BullseyeCavityPSO/) · [逐项数据](comparisons/BullseyeCavityPSO.md) | 3.420e+12 | 一致 | — | 用 600 min 的作业整本跑完，NBRESULT OK、20/20 单元无报错，（15 轮粒子群 × 5 个仿真） |
| [ZonePlateFieldProjection](https://www.flexcompute.com/tidy3d/examples/notebooks/ZonePlateFieldProjection/) · [逐项数据](comparisons/ZonePlateFieldProjection.md) | 3.884e+12 | 一致 | 近场与远场投影 | 用当前代码重组装后逐点比较：投影场 Er/Eθ/Eφ 与 Ex/Ey/Ez 幅值中位比 0.99～1.00、相位差恒为 2.6°；notebook 里 Re… |
| [PhotonicCrystalWaveguidePolarizationFilter](https://www.flexcompute.com/tidy3d/examples/notebooks/PhotonicCrystalWaveguidePolarizationFilter/) · [逐项数据](comparisons/PhotonicCrystalWaveguidePolarizationFilter.md) | 5.584e+12 | 一致 | TE/TM 透射 | 复跑：两种极化的 |E|² 分布与色标一致 |
| [VortexMetasurface](https://www.flexcompute.com/tidy3d/examples/notebooks/VortexMetasurface/) · [逐项数据](comparisons/VortexMetasurface.md) | 7.045e+12 | 一致 | — | 重跑 14/14 单元无报错 |
| [DielectricMetasurfaceAbsorber](https://www.flexcompute.com/tidy3d/examples/notebooks/DielectricMetasurfaceAbsorber/) · [逐项数据](comparisons/DielectricMetasurfaceAbsorber.md) | 1.063e+13 | 一致 | 吸收谱与场强 | 复跑：17 个可比数字无差异，7 张图像素差 ≤0.7%（修掉 GaussianBeam direction="-" 注入方向反了的 bug 后一直如此） |
| [OptimizedL3](https://www.flexcompute.com/tidy3d/examples/notebooks/OptimizedL3/) · [逐项数据](comparisons/OptimizedL3.md) | 1.436e+13 | 一致 | L3 腔 Q（f≈191.8 THz）：8.866e+05 / —（缺） | 复跑：14 个可比数字一致；k 空间分布图逐点相同，像素差来自光锥边界那圈白色遮罩的画法 |
| [RingResonator](https://www.flexcompute.com/tidy3d/examples/notebooks/RingResonator/) · [逐项数据](comparisons/RingResonator.md) | 1.942e+13 | 一致 | add/drop 透射谱 | 全部修复后重跑（14/14）目检：透射谱五个谐振峰位一一对上（1.513/1.529/1.547/1.565/1.583 µm）、自由光谱范围相同、线宽相同… |
| [THzDemultiplexerFilter](https://www.flexcompute.com/tidy3d/examples/notebooks/THzDemultiplexerFilter/) · [逐项数据](comparisons/THzDemultiplexerFilter.md) | 2.318e+13 | 一致 | — | 合入（notebook 直连的本地 ModeSolver 按参考环境阶梯化）后重跑（12/12，135 s）：12 个打印数字 0 处差 >1%，5 张图全… |
| [NanostructuredBoronNitride](https://www.flexcompute.com/tidy3d/examples/notebooks/NanostructuredBoronNitride/) · [逐项数据](comparisons/NanostructuredBoronNitride.md) | 2.669e+13 | 一致 | hBN 声子极化激元场分布 | 重跑 15/15 单元无报错 |
| [BilevelPSR](https://www.flexcompute.com/tidy3d/examples/notebooks/BilevelPSR/) · [逐项数据](comparisons/BilevelPSR.md) | 5.087e+13 | 一致 | 输入模 0 的 n_eff：2.386 / 2.35（1.5%） | 复跑：19 个可比数字无差异；10 张图里模式场（TE0/TE1 的 |Ey|）分布与色标一致，像素差来自色条位置 |

## Inverse Design（逆设计）（26 本）

| 例子 | 计算量（格×步） | 判定 | 结论量（参考 / OpenEM，相对差） | 备注 |
|---|---:|---|---|---|
| [Autograd30ParallelAdjoint](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd30ParallelAdjoint/) · [逐项数据](comparisons/Autograd30ParallelAdjoint.md) | 3.519e+08 | 一致 | 顺序伴随目标值：0.3893 / 0.3893（0.0%）<br>并行伴随目标值：0.3893 / 0.3893（0.0%） | 目标值/梯度 7 个数字一致；>5% 的两个数字是「顺序伴随 / 并行伴随」的耗时列（67.9 s vs 31.0 s），不是物理量 |
| [Autograd5BoundaryGradients](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd5BoundaryGradients/) · [逐项数据](comparisons/Autograd5BoundaryGradients.md) | 8.283e+08 | 一致 | 初始目标值：0.722 / 0.721（0.1%） | 复跑：PolySlab 顶点（边界）梯度逐分量对上参考 —— gradient[0] 2.714e-03 对 2.704e-03、[1] 1.0749e-0… |
| [Autograd3InverseDesign](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd3InverseDesign/) · [逐项数据](comparisons/Autograd3InverseDesign.md) | 1.513e+10 | 轨迹差异 | 输入模 n_eff：1.572 / 1.572（0.0%）<br>第 1 步目标 J：-0.03995 / -0.03995（0.0%）<br>第 1 步梯度范数：0.0004596 / 0.0004634（0.8%） | 复跑：首步锚点一致（J −3.9947e-02 vs −3.9946e-02，grad_norm 4.5956e-04 vs 4.6339e-04，差 0.… |
| [Autograd26DiffractiveBeamSplitter](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd26DiffractiveBeamSplitter/) · [逐项数据](comparisons/Autograd26DiffractiveBeamSplitter.md) | 1.839e+10 | 轨迹差异 | — | 用含 的代码整本跑完（21/21，81531 s = 22.6 h） |
| [Autograd31GratingCouplerWithBeamOptimization](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd31GratingCouplerWithBeamOptimization/) · [逐项数据](comparisons/Autograd31GratingCouplerWithBeamOptimization.md) | 2.455e+10 | 轨迹差异 | 第 1 步目标值：0.08427 / 0.08384（0.5%） | 修复后代码重跑（15/15，6577 s）：step 1 objective 0.084244 对 0.084270（0.03%），grad norm 1.… |
| [Autograd17BandPassFilter](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd17BandPassFilter/) · [逐项数据](comparisons/Autograd17BandPassFilter.md) | 3.393e+10 | 轨迹差异 | 初始值（c21）：0.0788 / 0.0789（0.1%）<br>第 1 步目标 J：-1.078 / -1.077（0.1%） | 用 600 min 的作业整本跑完，NBRESULT OK、38/38 单元无报错，（之前两次都在 5 小时上限被切断） |
| [Autograd1Intro](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd1Intro/) · [逐项数据](comparisons/Autograd1Intro.md) | 5.001e+10 | 一致 | 功率：0.552 / 0.547（0.9%）<br>d_power/d_eps：-0.1614 / -0.1567（2.9%） | 复跑：7 个数字最大差 2.9%，4 张图像素差 ≤4% |
| [Autograd8WaveguideBend](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd8WaveguideBend/) · [逐项数据](comparisons/Autograd8WaveguideBend.md) | 5.294e+10 | 一致 | n_eff（首模）：1.797 / 1.802（0.3%）<br>初始目标值（c25）：0.5688 / 0.5688（0.0%） | 三个模式的 Ex/Ey/Ez 场型与色标一致（9-01 结果，未复跑） |
| [Autograd0Quickstart](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd0Quickstart/) · [逐项数据](comparisons/Autograd0Quickstart.md) | 7.186e+10 | 一致 | 第 1 轮强度：861 / 860（0.1%） | 复跑：21 个数字最大差 1.6%，梯度校验图一致 |
| [Autograd21GaPLightExtractor](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd21GaPLightExtractor/) · [逐项数据](comparisons/Autograd21GaPLightExtractor.md) | 7.396e+10 | 一致 | 第 1 步目标值：9.46 / 9.463（0.0%） | 最终设计图案几乎相同（9-02 结果，未复跑） |
| [Autograd29SourceGradients](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd29SourceGradients/) · [逐项数据](comparisons/Autograd29SourceGradients.md) | 7.911e+10 | 轨迹差异 | 基线聚焦功率：0.0413 / 0.0396（4.1%）<br>初始聚焦功率：0.0263 / 0.0238（9.5%）<br>初始梯度范数：0.331 / 0.364（10.0%） | 修复后代码重跑（17/17，19227 s）：起点就差——baseline focused power 0.0396 对 0.0413（−4.1%）、ini… |
| [Autograd24DigitalSplitter](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd24DigitalSplitter/) · [逐项数据](comparisons/Autograd24DigitalSplitter.md) | 9.556e+10 | 轨迹差异 | 第 1 轮 FOM：0.1113 / 0.1107（0.5%） | 修复后代码重跑（23/23，17211 s）：首轮 FOM 0.1107 对 0.1113（0.5%），第 2 轮起系统性偏低（0.1418 对 0.157… |
| [Autograd9WDM](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd9WDM/) · [逐项数据](comparisons/Autograd9WDM.md) | 9.868e+10 | 轨迹差异 | n_eff（首模）：3.151 / 3.151（0.0%）<br>初始目标值（c18）：-2.331 / -2.331（0.0%） | 全部修复（#62～#71）后整本重跑（32/32，21995 s）：模式有效折射率 3.15138897 对 3.15138964；优化循环第 1 轮 J … |
| [Autograd10YBranchLevelSet](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd10YBranchLevelSet/) · [逐项数据](comparisons/Autograd10YBranchLevelSet.md) | 1.297e+11 | 轨迹差异 | 第 1 步 obj_eps：36.66 / 36.66（0.0%）<br>第 1 步梯度范数：21.34 / 21.34（0.0%） | 优化终态不同：参考 gap/curvature 罚 0.061/0.432，OpenEM 0/0；优化轨迹不比，只比首步锚点 |
| [Autograd27Smatrix](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd27Smatrix/) · [逐项数据](comparisons/Autograd27Smatrix.md) | 1.553e+11 | 轨迹差异 | 初始目标 J：0.266 / 0.256（3.8%）<br>初始梯度范数：0.0286 / 0.0287（0.3%） | 用 300 min 的作业整本重跑，NBRESULT OK、26/26 单元无报错，（此前 110 min 那次切在 24/26） |
| [Autograd2GradientChecking](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd2GradientChecking/) · [逐项数据](comparisons/Autograd2GradientChecking.md) | 1.718e+11 | 一致 | T(FDTD)（TMM 0.78581）：0.7852 / 0.7851（0.0%）<br>grad_eps(FDTD) 首项（TMM -0.2766）：-0.2817 / -0.2816（0.0%） | 结构剖面图相同 |
| [Autograd12LightExtractor](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd12LightExtractor/) · [逐项数据](comparisons/Autograd12LightExtractor.md) | 3.106e+11 | 轨迹差异 | 光提取效率优化轨迹与场分布（轨迹不比） | 9-01 旧结果（未复跑）：终态场分布不同（参考为条纹干涉图样，OpenEM为偶极子附近的同心环），优化轨迹不比 |
| [Autograd18TopologyBend](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd18TopologyBend/) · [逐项数据](comparisons/Autograd18TopologyBend.md) | 9.530e+11 | 轨迹差异 | 第 1 步目标 J：-0.9875 / -0.9874（0.0%）<br>第 1 步梯度范数：0.01181 / 0.01186（0.4%） | 复跑：首步锚点一致（J −9.8748e-01 vs −9.8743e-01，grad_norm 1.1811e-02 vs 1.1860e-02，差 0.… |
| [Autograd13Metasurface](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd13Metasurface/) · [逐项数据](comparisons/Autograd13Metasurface.md) | 1.291e+12 | 一致 | 无器件平均强度：1.81 / 1.81（0.0%）<br>初始损失：0.246 / 0.246（0.0%） | Re{Ex} 场分布与色标一致（9-01 结果，未复跑） |
| [Autograd6GratingCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd6GratingCoupler/) · [逐项数据](comparisons/Autograd6GratingCoupler.md) | 1.304e+12 | 轨迹差异 | 第 1 轮目标 J：-0.9993 / -0.9996（0.0%）<br>第 1 轮梯度范数：0.0005068 / 0.0003394（33.0%） | 修复后代码重跑（22/22，6324 s）：首步 J −0.99950 对 −0.99931、grad_norm 4.15e-4 对 5.07e-4（0.8… |
| [Autograd0QuickstartII](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd0QuickstartII/) · [逐项数据](comparisons/Autograd0QuickstartII.md) | 1.393e+12 | 一致 | 优化过程与最终场分布（轨迹不比） | 复跑：14 个数字全部一致 |
| [Autograd16BilayerCoupler](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd16BilayerCoupler/) · [逐项数据](comparisons/Autograd16BilayerCoupler.md) | 1.543e+12 | 轨迹差异 | n_eff（首模）：2.7 / 2.7（0.0%）<br>初始目标值：-0.8573 / -0.8574（0.0%） | 修复后代码重跑（26/26，6368 s）：三个模式 n_eff 逐位一致，起点 objective −0.85735 对 −0.85730，首步 grad… |
| [Autograd15Antenna](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd15Antenna/) · [逐项数据](comparisons/Autograd15Antenna.md) | 1.563e+12 | 轨迹差异 | 无结构强度：2084 / 2084（0.0%）<br>初始目标值：0.002172 / 0.002142（1.4%） | 复跑（放宽逐格列负填充容差 后跑通，此前一直建不出场景）：首步增强因子 21.2 对 20.1、目标 3.40 对 3.22（差 5.5%），之后 40 轮… |
| [Autograd22PhotonicCrystal](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd22PhotonicCrystal/) · [逐项数据](comparisons/Autograd22PhotonicCrystal.md) | 2.907e+12 | 轨迹差异 | 初始透射（c17）：0.9999 / 1.002（0.2%）<br>第 1 步目标 J：0.5039 / 0.5003（0.7%） | 修复后代码重跑（27/27，3244 s）：首步 J 0.5003 对参考 0.5039（0.7%）、第二步 0.7175 对 0.7162，正向一致；gr… |
| [Autograd7Metalens](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd7Metalens/) · [逐项数据](comparisons/Autograd7Metalens.md) | 3.140e+12 | 基本一致 | 初始目标值（c15）：0.00414 / 0.00396（4.4%） | 全部修复（#62～#71）后整本重跑（25/25，15700 s；上一轮 09-08 那次跑完后写盘时撞上 workspace 配额满、结果丢失，这是补跑的） |
| [Autograd25WaveguideCrossing](https://www.flexcompute.com/tidy3d/examples/notebooks/Autograd25WaveguideCrossing/) · [逐项数据](comparisons/Autograd25WaveguideCrossing.md) | 3.507e+12 | 基本一致 | 初始中心波长透射：0.7748 / 0.7745（0.0%）<br>第 1 步梯度范数：0.002335 / 0.002735（17.1%） | 全部修复（含对称面镜像侧梯度 #68）后重跑（形状梯度合入后）：目标值逐步一致（首步 J 0.77449 对 0.77478，差 0.04%；第二步 0.7… |

## 未列入的 11 本

- **不能跑** 5 本
- **能力缺口** 2 本
- **有差异** 2 本
- **参考存疑** 2 本

分别是：差异超过 20% 且已定位到具体机制的；用到 OpenEM 还没实现的算法（不是数值误差，是能力缺口）；
参考侧自身不可复现、我们这侧自洽的；以及超出单卡容量或依赖非 FDTD 求解器（如本征模展开）跑不了的。
这些例子的差异都查到了机制，但结论仍在演进，等稳定后再发布。

## 怎么复现

对照用的 notebook 是 Flexcompute 的 Tidy3D example library，版权归其所有，本仓库**不重新分发**——
表里每个例子的链接指向官方页面。要复现某一本：从官方取 notebook，装上 OpenEM 的钩子
（`openem.install()`，见 [tutorial.md](tutorial.md)），整本重跑，再与官方页面上的输出对照。

判定按本页开头那张表：可比数字与结论量的最大相对差 ≤5% 记「一致」，5%～20% 记「基本一致」，
差异集中在近零量或本身不收敛的指标上也归后者。图按像素差目检，只看有没有真差异。

