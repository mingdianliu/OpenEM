# 上手教程

三步：平面波打介质板读通量并与解析公式对照、波导的模式源与模式监视器、用 autograd 对一个几何参数
求梯度。每一步的代码都能直接存成 `.py` 跑，输出是实跑贴回来的。

装好之后，唯一要记住的就是这一行：

```python
import openem
openem.install()      # 之后 td.web.run / web.Batch / autograd 的 grad 都落到本地 GPU
```

原来的 Tidy3D 脚本或 notebook 一行都不用改。不想打补丁的话，用 `openem.run(sim)` 直接解一个
`Simulation`，返回真正的 `td.SimulationData`。

**下面的数字是在什么环境上跑的**：单张 NVIDIA H800、tidy3d 2.12.0、cupy 13.3.0，并且装了
`tidy3d-extras`（客户端能做 subpixel 平均）。换硬件末位会有差异。**没装 `tidy3d-extras` 时，
OpenEM 拿不到平均后的 ε 张量，会打印一行提示并退回阶梯化**，第二步的有效折射率会明显不同——
具体差多少那一节里给了。

---

## 第一步：平面波打介质板

厚 0.3 µm、ε=4（n=2）的无损平板泡在真空里；x/y 周期边界、z 两侧 PML；平面波从 −z 侧正入射，
平板前后各一个通量监视器，中间一个 xz 面的频域场监视器。完整文件是
[`examples/01_slab_planewave.py`](../examples/01_slab_planewave.py)。

```python
import numpy as np
import tidy3d as td

import openem

C0 = 2.99792458e8          # 真空光速 m/s
LAM0 = 1.0                 # 中心波长 µm
N_SLAB, D_SLAB = 2.0, 0.3  # 平板折射率与厚度 µm


def analytic_transmission(freqs_hz):
    """空气/介质/空气单层薄膜、正入射的解析透射率（无损，故 T = 1 − |r|²）。"""
    lam_um = C0 / np.asarray(freqs_hz) * 1e6
    r1 = (1.0 - N_SLAB) / (1.0 + N_SLAB)          # 前界面；后界面是 −r1
    phase = np.exp(2j * (2 * np.pi * N_SLAB * D_SLAB / lam_um))
    r = (r1 - r1 * phase) / (1.0 - r1 * r1 * phase)
    return 1.0 - np.abs(r) ** 2


def main():
    openem.install()                       # 之后 td.web.run 落到本地 GPU

    f0 = C0 / (LAM0 * 1e-6)
    freqs = np.linspace(0.85 * f0, 1.15 * f0, 11)
    periodic = td.Boundary(plus=td.Periodic(), minus=td.Periodic())

    sim = td.Simulation(
        size=(0.6, 0.6, 3.0),
        grid_spec=td.GridSpec.uniform(dl=0.02),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, td.inf, D_SLAB)),
            medium=td.Medium(permittivity=N_SLAB ** 2))],
        sources=[td.PlaneWave(
            center=(0, 0, -1.0), size=(td.inf, td.inf, 0), direction="+",
            source_time=td.GaussianPulse(freq0=f0, fwidth=0.2 * f0))],
        monitors=[
            td.FluxMonitor(center=(0, 0, 1.0), size=(td.inf, td.inf, 0),
                           freqs=list(freqs), name="T"),
            td.FluxMonitor(center=(0, 0, -1.2), size=(td.inf, td.inf, 0),
                           freqs=list(freqs), name="R"),
            td.FieldMonitor(center=(0, 0, 0), size=(td.inf, 0, td.inf),
                            freqs=[f0], name="fxz"),
        ],
        run_time=4e-13,
        boundary_spec=td.BoundarySpec(x=periodic, y=periodic, z=td.Boundary.pml()))

    sim_data = td.web.run(sim, task_name="slab_planewave")

    T = np.asarray(sim_data["T"].flux.values, dtype=float)
    R = -np.asarray(sim_data["R"].flux.values, dtype=float)   # 反射波朝 −z，通量为负
    T_ref = analytic_transmission(freqs)
    i0 = len(freqs) // 2                                      # 正中间就是 f0

    print()
    print(f"格数           : {np.prod(sim.grid.num_cells)}  {tuple(sim.grid.num_cells)}")
    print(f"T(λ=1µm)       : {T[i0]:.4f}   解析 {T_ref[i0]:.4f}")
    print(f"R(λ=1µm)       : {R[i0]:.4f}   解析 {1 - T_ref[i0]:.4f}")
    print(f"T+R 全频点     : min {(T + R).min():.4f}  max {(T + R).max():.4f}")
    print(f"|T−解析|/解析  : max {np.abs(T - T_ref).max() / T_ref.max():.2e}")
    Ex = np.asarray(sim_data["fxz"].Ex.values)
    print(f"xz 面 max|Ex|  : {np.abs(Ex).max():.4f}   数组形状 {Ex.shape}")


if __name__ == "__main__":
    main()
```

跑出来：

```
格数           : 156600  (30, 30, 174)
T(λ=1µm)       : 0.8310   解析 0.8373
R(λ=1µm)       : 0.1689   解析 0.1627
T+R 全频点     : min 1.0000  max 1.0001
|T−解析|/解析  : max 8.21e-03
xz 面 max|Ex|  : 64.5887   数组形状 (31, 1, 175, 1)
```

怎么读这几个数：

- **先看 T+R。** 无损介质，能量必须守恒，11 个频点上 T+R 都是 1.0000～1.0001。这一条不过，
  后面所有数字都不用看了——它同时检查了源注入、两个通量面的口径和符号。
- **再看 T 对解析值。** 0.8310 对 0.8373，相对差 0.75%（全频点最大 0.82%）。这个差来自网格：
  `dl=0.02` µm，板内波长是 λ/n = 0.5 µm，一个波长只有 25 格。把 `dl` 改成 `0.01` 重跑，
  T = **0.8372**，解析 0.8373——离散误差按预期收敛掉了。
- **频域场的绝对值不用对着 1 看。** `max|Ex| = 64.6` 是按源谱归一化之后的量，大小取决于脉冲谱本身。

---

## 第二步：波导的模式源与模式监视器

ε=12.1 的硅条波导（宽 0.5 µm、厚 0.22 µm）在真空里，六面 PML；x=−0.6 µm 处一个 `ModeSource`
注入基模，x=+0.5 µm 处一个 `ModeMonitor`（2 个模式、3 个频点）和一个 `FluxMonitor`，
中间一个 `ModeSolverMonitor` 读模式本身。完整文件是
[`examples/02_waveguide_mode.py`](../examples/02_waveguide_mode.py)。

```python
import numpy as np
import tidy3d as td

import openem

F0 = 1.934e14                 # ≈ 1.55 µm


def main():
    openem.install()

    mode_spec = td.ModeSpec(num_modes=2)
    plane = dict(size=(0, 1.8, 1.4))
    sim = td.Simulation(
        size=(2.0, 2.4, 2.0),
        grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[td.Structure(
            geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
            medium=td.Medium(permittivity=12.1))],
        sources=[td.ModeSource(
            center=(-0.6, 0, 0), **plane, direction="+", mode_index=0,
            mode_spec=td.ModeSpec(num_modes=1),
            source_time=td.GaussianPulse(freq0=F0, fwidth=0.1 * F0))],
        monitors=[
            td.ModeMonitor(center=(0.5, 0, 0), **plane, mode_spec=mode_spec,
                           freqs=[0.98 * F0, F0, 1.02 * F0], name="m"),
            td.FluxMonitor(center=(0.5, 0, 0), **plane, freqs=[F0], name="fl"),
            td.ModeSolverMonitor(center=(0.0, 0, 0), **plane, mode_spec=mode_spec,
                                 freqs=[F0], name="ms"),
        ],
        run_time=2e-13,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))

    sim_data = td.web.run(sim, task_name="waveguide_mode")

    amps = sim_data["m"].amps.sel(direction="+")
    p_fwd = np.abs(np.asarray(amps.values)) ** 2          # (f, mode_index)
    flux = float(np.real(sim_data["fl"].flux.values[0]))
    n_eff = np.real(np.asarray(sim_data["ms"].n_eff.values)).ravel()

    print()
    print(f"格数              : {np.prod(sim.grid.num_cells)}  {tuple(sim.grid.num_cells)}")
    print(f"n_eff (f0)        : mode0 {n_eff[0]:.4f}   mode1 {n_eff[1]:.4f}")
    print(f"|amp+|² mode0     : " + "  ".join(f"{v:.4f}" for v in p_fwd[:, 0]))
    print(f"|amp+|² mode1     : " + "  ".join(f"{v:.3e}" for v in p_fwd[:, 1]))
    print(f"总通量 (f0)       : {flux:.4f}")
    print(f"|amp0|²/通量 (f0) : {p_fwd[1, 0] / flux:.4f}")


if __name__ == "__main__":
    main()
```

跑出来：

```
格数              : 200704  (56, 64, 56)
n_eff (f0)        : mode0 2.3900   mode1 1.3080
|amp+|² mode0     : 0.9997  0.9994  0.9988
|amp+|² mode1     : 3.743e-19  1.693e-19  1.018e-19
总通量 (f0)       : 0.9994
|amp0|²/通量 (f0) : 1.0000
```

怎么读：

- **n_eff 落在包层和芯层之间。** 芯层 √12.1 = 3.478，包层 1，基模 2.3900、一阶模 1.3080，
  基模更靠近芯层折射率——导波模式本该如此。
- **`|amp+|²` mode0 ≈ 1。** 直波导无损、无散射，注入多少功率就到达多少。三个频点从 0.9997 掉到
  0.9988，是带边的正常衰减。这个数同时说明模式源的注入幅值和模式分解的口径是一致的。
- **一阶模是 1e-19。** 源只注入 mode0，直波导不会把功率耦合到高阶模，剩下的就是数值零。
- **模式分解与坡印廷通量对得上。** 总通量 0.9994，`|amp0|²/通量 = 1.0000`。
- **没装 `tidy3d-extras` 的话数字会变。** 同一段代码在退回阶梯化的环境里给 n_eff mode0 =
  **2.5775**、`|amp+|²` = 0.9996。波导厚 0.22 µm 在 `dl=0.05` 的网格上只有 4.4 格，阶梯化会把
  芯层算胖，n_eff 因此偏高。这不是求解器的自由度：ε 数组是 Tidy3D 客户端给的，OpenEM 只是照单用。

---

## 第三步：用 autograd 对一个几何参数求梯度

场景：一根宽 0.5 µm 的硅波导，中间插一段长 0.6 µm、宽 `w` 的展宽段。`w` 偏离 0.5 µm 会引起模式
失配，透射跟着变。目标 J 是模式监视器上前向基模的 `|amp|²`，要的是 dJ/dw。

装了钩子之后，写法就是 tidy3d 官方 autograd 的写法——`autograd.value_and_grad` 包住一个
「造仿真 → `td.web.run` → 取标量」的函数，正向和伴随都在本地 GPU 上跑：

```python
import autograd as ag
import autograd.numpy as anp
import tidy3d as td

import openem

F0 = 1.934e14                                   # ≈ 1.55 µm
WG = td.Structure(geometry=td.Box(center=(0, 0, 0), size=(td.inf, 0.5, 0.22)),
                  medium=td.Medium(permittivity=12.1))


def make_sim(w):
    """中间一段宽度为 w 的展宽段；其余部分是宽 0.5 µm 的直波导。"""
    seg = td.Structure(geometry=td.Box(center=(0, 0, 0), size=(0.6, w, 0.22)),
                       medium=td.Medium(permittivity=12.1))
    return td.Simulation(
        size=(2.4, 2.4, 1.4), grid_spec=td.GridSpec.uniform(dl=0.05),
        structures=[WG, seg],
        sources=[td.ModeSource(center=(-0.8, 0, 0), size=(0, 1.6, 1.0), direction="+",
                               source_time=td.GaussianPulse(freq0=F0, fwidth=0.1 * F0),
                               mode_spec=td.ModeSpec(num_modes=1), mode_index=0)],
        monitors=[td.ModeMonitor(center=(0.8, 0, 0), size=(0, 1.6, 1.0), freqs=[F0],
                                 mode_spec=td.ModeSpec(num_modes=1), name="m")],
        run_time=2.0e-13, shutoff=0.0,
        boundary_spec=td.BoundarySpec.all_sides(td.PML(num_layers=8)))


def objective(w):
    sim_data = td.web.run(make_sim(w), task_name="widen", verbose=False)
    amp = sim_data["m"].amps.sel(direction="+", f=F0, mode_index=0)
    return anp.sum(anp.abs(amp.values) ** 2)


openem.install()

w0 = 0.90
J0, g = ag.value_and_grad(objective)(w0)
w1 = w0 + 0.01 * float(g)          # 沿梯度走一步；g 是负的，所以 w 变小
J1 = objective(w1)
print()
print(f"J(w={w0:.4f})       : {float(J0):.6f}")
print(f"dJ/dw             : {float(g):+.5f}")
print(f"w <- w + 0.01*g   : {w1:.4f}")
print(f"J(w={w1:.4f})       : {float(J1):.6f}")
print(f"两点连线的斜率    : {(float(J1) - float(J0)) / (w1 - w0):+.2f}")
```

跑出来：

```
J(w=0.9000)       : 0.684555
dJ/dw             : -3.93767
w <- w + 0.01*g   : 0.8606
J(w=0.8606)       : 0.851527
两点连线的斜率    : -4.24
```

怎么读：

- **一次 `value_and_grad` 里跑了两个仿真**：正向一个，伴随一个（伴随源由 tidy3d 客户端造，
  步数照抄正向）。参数再多，代价也还是这两个仿真——这正是伴随法的意义。
- **梯度是负的**，说明在 w=0.9 µm 附近把这段变窄能提高透射。沿梯度走一步，J 确实升上去了。
- **两点连线的斜率和梯度同号、同量级**（差不到 10%）。这是这类校验能给到的程度，别指望更准，
  原因见下一节。

---

## 常见问题

### 抛了 `NotImplementedError`，是坏了吗

不是，是**故意的**。OpenEM 对不支持的物理一律 fail closed：场景构建器直接报错，不会悄悄换个
近似做法接着算。报错信息里会写清楚这个位置支持哪些类型，例如：

```
NotImplementedError: 只支持 PlaneWave / PointDipole / UniformCurrentSource /
CustomCurrentSource / CustomFieldSource / ModeSource / GaussianBeam / TFSF，收到 ...
```

为什么这么设计：静默降级会让后面所有的偏差都归不到因——你分不清一个 8% 的差异是求解器的问题，
还是某处悄悄换了个近似。宁可停下。README 的「支持范围」表列了全部支持项，表外的都会这样停。

有两种情况会**降级但打印出来**，不是报错：拿不到 `tidy3d-extras` 时退回阶梯化
（`[openem] subpixel ... 退回 subpixel=False 重建`），以及 `LossyMetalMedium` 按 PEC 阶梯化处理。
两条都会在 stdout 上说明。

`install()` 默认还有一道闸门：超过 8000 万格的仿真直接拒绝，防止误提巨型场景。真要跑大的，
用 `openem.install(max_cells=...)` 抬高上限。

### 显存不够

现象是 cupy 抛 `OutOfMemoryError`。能做的事，按性价比排：

1. **减格数**——加大 `dl` 或缩小仿真区。场数组随格数线性增长，这是最直接的一项。
2. **用对称性**。`td.Simulation(symmetry=...)` 声明的镜像对称会被折成半域 / 四分之一域 /
   八分之一域，显存和时间一起省。
3. **减频域监视器**。每个频点、每个分量、每个格子都要一份复数缓冲，宽频带的三维场监视器很吃显存。
   频点数和监视器体积都可以砍。
4. **`OPENEM_SUBPROCESS=1`**（只对 `openem.run(...)` 那条路有效，见下一条）。一个脚本里连着解几十个
   仿真时，进程内的显存不一定及时还回去；设了这个开关，每次求解起一个子进程，结束就还干净。
5. `install(max_cells=...)` 只是拒收的闸门，调大它不省显存。

### 工作目录和缓存：结果会被复用吗

要看你走的是哪条路，两条路行为不一样：

- **装了钩子之后的 `td.web.run`** 每次都真算。它在进程里直接建场景、跑求解器，不落中间文件，
  所以同一个仿真调两次就是解两次。
- **`openem.run(sim)`** 走的是另一条路（`nb/backend.py`），会在工作目录里留下 `scene.npz`（输入）
  和 `ours.npz`（原始输出），而且 **`ours.npz` 比 `scene.npz` 新就直接复用、不重算**。
  工作目录根由 `OPENEM_WORK_DIR` 指定，默认是当前目录下的 `openem_runs/`。
  反过来说：**改了求解器代码但场景没变，拿到的还是旧结果**——要强制重算就把对应的目录
  （或整个 `openem_runs/`）删掉。这条路还认 `OPENEM_SUBPROCESS=1`：每次求解起一个子进程，
  长脚本里能及时把显存还回去。

两条路都会缓存模式基底（日志里的「模式缓存 命中 / miss」就是它）。`OPENEM_MODE_CACHE=0` 关掉，
`OPENEM_MODE_CACHE_DIR` 指定存放位置。

包里所有 `OPENEM_*` 开关都在 `openem/knobs.py` 登记了默认值和用途；性能类的开关切换后结果
逐位不变。

### 几何参数的梯度，能用有限差分校验吗

只能粗校验，别当判据。原因有两层：

- 网格是不动的，几何一变只改变 ε 数组。**在没有 subpixel 平均的环境里，ε 图对几何参数是阶梯函数**
  ——把宽度从 1.000 µm 挪到 1.005 µm，在 `dl=0.05` 的网格上光栅化出来是同一份 ε，目标函数一个数都不变，
  有限差分会给 0。
- 有了 subpixel 平均，ε 随几何连续变化了，但它在格内**不是线性的**，所以有限差分的值还依赖你取的
  步长和界面落在格子里的位置。

OpenEM 的几何梯度算的是**离散目标函数的斜率**（ε 图扰动的离散链式法则），这正是优化器需要的量。
上面第三步用一个 0.04 µm 的步长做两点斜率，和伴随梯度差不到 10%，这个量级的一致就说明梯度是对的。
要更严格的梯度校验，看官方的 `Autograd2GradientChecking`——它是对介电常数（连续参数）做的，
那种校验才能对到几位。

### 还想看什么

- 128 本官方 notebook 的逐本对照结果：[validation.md](validation.md)
- 代码结构与模块职责：[architecture.md](architecture.md)
