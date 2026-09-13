# OpenEM

[English](README.md) · [简体中文](README.zh-CN.md)

一个 CUDA FDTD 求解器：输入是 Tidy3D 的 `Simulation` 对象，在本地单张 GPU 上解，fp32、Yee 网格。
装上钩子（`openem.install()`）之后，`tidy3d.web.run(sim)` 和 autograd 的 `grad` 都落到这张卡上，
原来的 notebook 一行都不用改。

> **网格和材料来自 Tidy3D。** OpenEM 读进一个 `td.Simulation`，**网格、几何光栅化、材料数组直接取用
> Tidy3D 客户端已经算好的结果**，自己只做时间步进、边界、源注入和监视器读出。这是刻意的：两侧吃同一份
> 离散化，结果上的差异就只可能来自求解器本身，否则分不清是建模差异还是求解差异。我们自己的网格与材料
> 生成模块在后续开发计划里，不在这一版。

## 安装

```bash
git clone <仓库地址>
cd OpenEM
pip install -e ".[gpu]"      # 已经装好 cupy 的话，pip install -e . 就够
```

需要 Python ≥ 3.10、`tidy3d` ≥ 2.8、NumPy、SciPy，以及一张 CUDA 12 及以上的 NVIDIA GPU 和
`cupy-cuda12x` ≥ 13（`gpu` 这个 extra 装的就是它）。CUDA 内核由 cupy 在首次用到时编译，
装的时候不需要编译步骤，也不需要 nvcc。没有 GPU 的机器上包仍然能 import，不依赖 CUDA 的测试也能跑
（CI 就是这么配的），只是解不了仿真。

## 30 秒上手

平面波正入射打一块 0.3 µm 厚、ε=4 的介质板，用通量监视器读透射：

```python
import tidy3d as td
import openem

openem.install()                       # 之后 td.web.run 就在本地解

f0 = 2.998e14                          # 1 µm
periodic = td.Boundary(plus=td.Periodic(), minus=td.Periodic())
sim = td.Simulation(
    size=(0.6, 0.6, 3.0),
    grid_spec=td.GridSpec.uniform(dl=0.02),
    structures=[td.Structure(
        geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, 0.3)),
        medium=td.Medium(permittivity=4.0))],
    sources=[td.PlaneWave(
        center=(0, 0, -1.0), size=(td.inf, td.inf, 0), direction="+",
        source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
    monitors=[td.FluxMonitor(center=(0, 0, 1.0), size=(td.inf, td.inf, 0),
                             freqs=[f0], name="T")],
    run_time=4e-13,
    boundary_spec=td.BoundarySpec(x=periodic, y=periodic, z=td.Boundary.pml()))

sim_data = td.web.run(sim, task_name="slab")
print("transmission =", float(sim_data["T"].flux.values[0]))
```

```
transmission = 0.8310044973694679
```

同一块板，解析的单层薄膜公式给 0.8373。`examples/01_slab_planewave.py` 是这个仿真的 11 频点版本，
还会打印反射和能量守恒；[`docs/tutorial.md`](docs/tutorial.md) 把它、模式源和求梯度各讲了一遍。

`openem.install()` 同样接管 `web.Batch` 和 `web.Job`，用批量提交的 notebook 照跑不误；成本查询接口
返回一个占位数，不会去碰云端。不想打补丁的话，`openem.run(sim)` 直接解一个 `Simulation`
（也接受列表或 `{名字: sim}` 字典），返回真正的 `td.SimulationData`；这条路会把中间文件留在磁盘上
并且能复用，教程里有说明。

## 支持范围

| | |
|---|---|
| **介质** | `Medium`（可带 `conductivity`）、所有能归约成 `PoleResidue` 的色散介质（`Lorentz`、`Drude`、`Sellmeier`、`Debye`，以及它们逐格变化的 `Custom*` 形式）、`CustomMedium`、`AnisotropicMedium`（对角）、`FullyAnisotropicMedium`（全张量，阶梯化）、结构上的 `PECMedium`、经客户端折算成等效体的 `Medium2D`、时变介质（`modulation_spec`） |
| **边界** | `PML` / `StablePML`（CPML）、`Absorber`、`Periodic`、`BlochBoundary`、`PECBoundary`、`PMCBoundary`，以及镜像对称面（`symmetry`） |
| **源** | `PlaneWave`（正入射与斜入射）、`ModeSource`、`PointDipole`、`UniformCurrentSource`、`GaussianBeam`、`TFSF`、`CustomFieldSource`、`CustomCurrentSource` |
| **监视器** | `FluxMonitor`（面或盒）、`FluxTimeMonitor`、`FieldMonitor`、`FieldTimeMonitor`、`ModeMonitor`、`ModeSolverMonitor`、`PermittivityMonitor`、`DiffractionMonitor`，以及远场投影监视器（角度域与 k 空间） |
| **梯度** | tidy3d autograd 接口背后的离散伴随：介电常数、`CustomMedium` 的逐格像素，以及 `Box` / `Cylinder` / `PolySlab` / `GeometryGroup` 的几何参数 |

这张表之外的东西一律 fail closed——场景构建器直接抛 `NotImplementedError`，不会悄悄换个别的做法。

## 验证

- **128 本** Tidy3D 官方 example library 的 notebook，整本跑完、逐个数字逐张图与参考输出对照：
  **99 本一致**、**11 本基本一致**、**18 本只有轨迹差异**（非凸优化两侧走到不同设计，但正向结果和首步
  梯度对得上）。逐本的表和具体数字在 [`docs/validation.md`](docs/validation.md)，另有 11 本暂未发布，
  原因也在那页。
- **22 个金标准仿真**（`tools/golden.py`）与基线**按位比对**。只有通量目标的梯度带容差，理由见下面的限制。
- 测试套件逐块覆盖求解器（`pytest tests`）；没有 GPU 的机器上 CUDA 相关的用例自动跳过。
- 与参考实现的耗时对照（118 个仿真，按准备 / 步进 / 读出三段拆）：[`docs/performance.md`](docs/performance.md)。
- 「一致」是怎么判的，以及每类判定对应多大的偏差：见 [`docs/validation.md`](docs/validation.md)
  开头的判定口径表。

## 限制

- **单卡**。不支持多卡，也不支持多机。一个仿真必须放得进一张卡；`install()` 默认拒绝 8000 万格以上的
  场景，要放宽用 `install(max_cells=...)`。
- **不做网格与材料生成**，两者都来自 Tidy3D 客户端，见开头那段。subpixel 平均需要 `tidy3d-extras`；
  拿不到的时候场景构建器退回阶梯化，并且会把这件事打印出来。
- **fail closed**。不支持的物理直接抛 `NotImplementedError`，不静默降级——静默降级会让后面的偏差
  归因不到任何地方。
- **已知的一处非确定性**。通量目标的伴随源是 `CustomCurrentSource`，会展开成几千个点偶极子；批量注入的
  内核用 `atomicAdd` 累加，float 加法不满足结合律，到达次序一变末位就变。这类梯度只能按容差比
  （5 次运行实测相对差 4e-7～3.7e-6）。模式目标的梯度和全部正向数据都按位可复现。
- 全程 fp32。

## 文档

| | |
|---|---|
| [`docs/tutorial.md`](docs/tutorial.md) | 三步上手：平板透射对解析公式、波导模式源与模式监视器、用 autograd 求梯度 |
| [`docs/architecture.md`](docs/architecture.md) | 目录树、每个模块一行职责，以及测试锁住的两条 import 边界 |
| [`docs/validation.md`](docs/validation.md) | 128 本 notebook 的逐本对照 |
| [`docs/comparisons/`](docs/comparisons/) | 同一批对照的完整数据：6,003 个配对数字、1,125 张图，一本一页外加一份 JSON |
| [`docs/performance.md`](docs/performance.md) | 与参考实现的逐例耗时对照 |
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | 改求解器代码的规矩：金标准必须保持按位相同 |

## 许可证与引用

GPL-3.0 或更新版本，见 [LICENSE](LICENSE)。基于本项目的衍生作品同样要以 GPL 发布。如果 OpenEM 对你的工作有帮助，引用信息在
[CITATION.cff](CITATION.cff) 里。
