# 代码结构

一句话：`scene/` 把 Tidy3D 的 `Simulation` 翻译成一份纯 numpy 的场景数据，`solver.py` 带着
`kernels/*.cu` 把它在 GPU 上推完，`nb/` 负责接住 notebook 的调用并把结果组装回 `td.SimulationData`。

一次求解的数据流：

```
td.Simulation
   │  scene.from_simulation()          取用客户端算好的网格 / 光栅化 / 材料，单位 µm → 米
   ▼
 Scene（纯 numpy 的数据模型，model.py）
   │  serialize.save()                 → scene.npz（只有 openem.run 那条路会落盘）
   ▼
 solver.run()                          建表 → 显存分配上传 → 主循环 → 监视器读出
   │                                   （device.py 是唯一发射内核的地方）
   ▼
 原始相量 / 时域序列
   │  normalize + td_readout           零拟合换算到 Tidy3D 的口径
   ▼
 td.SimulationData                     nb/backend.py 组装，下游 API 原样可用
```

## 目录树

```
openem/
├── __init__.py            包入口：install() 装钩子，run() 直接解一个 Simulation；模块边界写在它的 docstring 里
│
├── scene/                 td.Simulation → Scene。允许 import tidy3d 的两处之一，也是单位换算的地方
│   ├── build.py           把 media / boundaries / sources / monitors 拼成一个 Scene
│   ├── media.py           ε 与电导率，以及材料类型的 fail closed
│   ├── poles.py           材料 → 极点表（只和材料有关，不涉及网格）
│   ├── weights.py         逐格反解几何权重（只和几何有关，不涉及极点）
│   ├── dispersion.py      极点表 + 几何权重 → 求解器要的 CSR 表
│   ├── boundaries.py      边界类型、CPML 系数、Absorber 电导率、边界拓扑的 fail closed
│   ├── sources.py         平面波单向注入、点偶极子、对称面的镜像源
│   ├── modes.py           模式源与 ModeMonitor：把平面上的场投影到导波模式
│   ├── monitors.py        通量 / 频域场 / 时域场三类监视器
│   ├── projection.py      远场投影监视器拆成六个面的近场监视器，变换本身交给客户端
│   ├── modulation.py      时变介质（modulation_spec）的提取与 fail closed
│   ├── subpixel.py        local subpixel 的全局开关
│   ├── normalize_td.py    「原始 DFT → 监视器数据」的源谱除数
│   ├── cache.py           磁盘缓存的公共件：数组哈希、缓存目录、原子写
│   └── _util.py           Yee 坐标、分量键、符号、平面索引这些小工具
│
├── nb/                    notebook 接入层。允许 import tidy3d 的另一处
│   ├── backend.py         web.run 的平替：建场景 → 求解 → 归一化 → 组装 td.SimulationData
│   ├── solve_worker.py    解一个已序列化的场景（scene.npz → ours.npz），也可当子进程用
│   ├── autograd_hook.py   钩子本体：接管 tidy3d 的正向 / 伴随入口，以及 Job / Batch 壳
│   └── shapegrad.py       几何参数的伴随梯度，改走 ε 图扰动的离散链式法则
│
├── solver.py              主循环：每步按 update_h → 注入 → update_e → 注入 → 色散 → 吸收 → 采样发射内核
├── device.py              CUDA 侧的唯一入口：内核编译、发射、内核级计时
├── device_tables.py       显存分配与上传：host 侧建好的表 pad 到 pitch、搬上设备
│
├── sources_setup.py       建表：平面波 / TFSF 盒 / 模式源 / 偶极子 / 时变介质
├── monitors_setup.py      建表：四类监视器的缓冲、DFT 相位表、发射几何、共享快照与批量采样
├── dispersion_setup.py    建表：色散条目重排（去重 / 分桶 / 推迟 / 内域分治）与吸收体表
├── setup_tables.py        上面三个加 pitch.py 的再导出 shim，老调用点继续可用
│
├── model.py               求解器的数据模型（Scene 及其各部分），不 import tidy3d
├── grid.py                非均匀 Yee 网格
├── coeffs.py              E/H 更新系数、PML 轴表、PEC 掩膜、cacb 查表、psi 开关
├── cpml.py                CPML 剖面与递推系数
├── fusion.py              各条融合与快路径的门槛和发射参数
├── pitch.py               数组 pitch / pad / 上传工具（z 维 pitch 的唯一出处）
├── readout.py             监视器的步内采样、批清算与读出
├── flux.py                相量 → 坡印廷矢量 → 平面通量
├── colocate.py            Yee 原位 → 监视器坐标的共位
├── fold.py                对称性折叠：按 symmetry 折成半域 / 四分之一域 / 八分之一域
├── normalize.py           换算到 Tidy3D 的归一化口径，零拟合
├── td_readout.py          场相量转成 Tidy3D 的输出口径（不 import tidy3d，也不 import cupy）
├── shutoff.py             早停判据，三个量任一达标即停
├── apodization.py         切趾：频域监视器 DFT 累加的时间窗
├── tfsf1d.py              TFSF 盒的 1D 辅助网格
├── tfsf_oblique.py        斜入射 TFSF：辅助线沿 k̂，盒壳按投影位置插值
├── waveform.py            源的时间波形与频谱，转调 Tidy3D 不自己实现
├── serialize.py           Scene 的 npz 存取
├── results.py             求解结果的数据类
├── knobs.py               56 个 OPENEM_* 环境变量的总表：默认值、类别、用途
├── units.py               长度单位
├── autograd_hook.py       转发 shim → nb/autograd_hook.py
├── shapegrad.py           转发 shim → nb/shapegrad.py
│
└── kernels/               17 个 .cu 文件，一共 59 个内核，全是 extern "C" __global__
    ├── yee.cu             通用非均匀 Yee 更新 + 三轴 CPML（每轴独立地是周期或吸收，介质可有损）
    ├── yee_lean.cu        内域瘦内核 + 边界壳分治
    ├── yee_habs.cu        update_h 与 H 族吸收衰减融成一个内核
    ├── yee_fused.cu       H 与 E 一趟做完（默认关，实测更慢）
    ├── yee_edge.cu        上一条的收尾内核
    ├── yee_ade.cu         极点更新内联进 update_e（默认关，实测更慢）
    ├── dispersion.cu      色散材料的辅助微分方程：每个复极点一个一阶方程，梯形法离散
    ├── dispersion_mix.cu  混合平均格：斜界面上算术平均与调和平均的加权和
    ├── tensor.cu          全张量 ε/σ 格的 E 更新，走「D 形式 + 逆张量」
    ├── source_dft.cu      单向平面波注入、伴随的批量点偶极子注入，以及频域监视器的运行时 DFT
    ├── tfsf.cu            TFSF 盒的六面入射修正，通用面内核
    ├── bloch.cu           Bloch（k≠0）场景的复数场内核变体
    ├── absorber.cu        绝热吸收体：按层渐变的匹配有损介质
    ├── absorb_batch.cu    吸收层衰减批量版，一次发射覆盖所有 slab
    ├── modulation.cu      时变介质：每步把被调制格子的 ca/cb 换成当步的值
    ├── dft_sub.cu         远场投影的面监视器只累加切向分量
    └── tmon_batch.cu      时域监视器批量采样，一次发射处理全部监视器
```

其余目录：`tests/`（单元测试与哨兵）、`tools/`（金标准、发布体检、两个文档生成器）、
`examples/`（可直接跑的小例子）、`docs/`。

## 两条 import 边界

这两条不是风格偏好，是测试锁着的硬约束，改动时不要绕过：

- **只有 `scene/` 与 `nb/` 能 import tidy3d**（`tests/test_no_tidy3d_leak.py` 锁住求解链路）。
  唯一的例外是 `normalize.py` 里一处函数内的延迟 import，只有「云端没做源归一化」那几例会走到。
  这样求解器本身不依赖 tidy3d 的版本，`td_readout` 之类也能在没有 `scene/` 的环境里单独跑。
- **只有求解链路的少数模块能 import cupy**：`solver`、`device`、`device_tables`、`sources_setup`、
  `monitors_setup`、`dispersion_setup`、`pitch`、`coeffs`、`fusion`、`readout`
  （`setup_tables` 只是它们的再导出 shim），`tfsf1d` 只在可选的 GPU 路径里局部 import。
  其余模块（`model` `grid` `cpml` `waveform` `serialize` `flux` `fold` `normalize` `colocate`
  `td_readout` `shutoff` `apodization` `tfsf_oblique` `results` `units` `knobs`）只见 numpy，
  在没有 GPU 的机器上照样能 import 和测试——CI 就是这么跑的。`tests/test_no_cupy_leak.py` 锁住。

再加一条：**发射内核只经 `device.py`**。以后要换掉 cupy、或者把主循环搬进 C++，只改这一个文件；
`kernels/*.cu` 是真 CUDA C，可以直接被独立的 C++ 程序编译。

## 接入层怎么工作

进来的路有两条，落在不同的模块上：

- **`openem.install()` 之后的 `td.web.run(sim)`** → `nb/autograd_hook.py`。它在进程内建场景、
  调 `solver.run`、把结果组装成 `td.SimulationData`，**不落中间文件、不复用上次结果**。
  `web.Job` / `web.Batch` 的 `run` / `load` / `start` / `monitor` 同样换成本地替身，
  `estimate_cost` / `real_cost` 返回占位值——本地跑没有云端账单。
- **`openem.run(sim)`** → `nb/backend.py`。它把场景存成 `scene.npz`、交给 `solve_worker` 解出
  `ours.npz`，再归一化并组装。工作目录默认是当前目录下的 `openem_runs/`（`OPENEM_WORK_DIR` 可改），
  `ours.npz` 比 `scene.npz` 新就直接复用不重算，`OPENEM_SUBPROCESS=1` 时每次求解起一个子进程。

两条路出来的都是真正的 `td.SimulationData`，`sim_data.plot_field`、`sim_data["mon"].Ex`
这些下游 API 原样可用。

autograd 的逆设计走 `nb/autograd_hook.py`：tidy3d 的 `local_gradient` 模式自己加监视器、
自己造伴随源、自己组装梯度，执行只从「正向」和「伴随」两个口下去；钩子把这两个口换成 OpenEM，
梯度组装仍然是 tidy3d 那套。几何参数的梯度由 `nb/shapegrad.py` 改走 ε 图扰动的离散链式法则
——它算的是**离散目标函数**的斜率，而不是连续形状导数的面积分。

## 开关

56 个 `OPENEM_*` 环境变量都在 `openem/knobs.py` 登记（默认值、类别、用途），分三类：

- `perf` 性能开关：切换后结果**逐位不变**，只影响速度和显存，每条都有哨兵测试守着；
- `policy` 口径开关：会改结果，只在对照实验里用；
- `debug` 诊断与计时打印，不改结果。

求解器侧和 `scene/`、`nb/` 侧的开关都经
`knobs.env` 现读，不是 import 时读一次。
