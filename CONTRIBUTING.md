# 参与开发

## 先跑通测试

```bash
pip install -e ".[dev]"
pytest tests -q -p no:randomly
```

没有 GPU 的机器上，依赖 CUDA 的用例会自动跳过 —— 纯 CPU 跑下来是 316 项里 168 通过、148 跳过，
**一条都不该失败**。跳过这么多是正常的：求解器主体要 CUDA，另有一批用例要 tidy3d-extras
（Flexcompute 的附加包，提供亚格平均后的 ε）。带 GPU 时全量约 375 项。

## 改动求解器代码的规矩

这个求解器的价值在于**它的输出与 Tidy3D 逐位对得上**（见 `docs/validation.md`）。所以任何改动都要回答一个问题：数值有没有变？

1. **只改结构、不改数值的改动**（重构、改名、抽函数）：跑 `tools/golden.py`，22 个小仿真必须与基线**按位相同**。

   ```bash
   python tools/golden.py run base      # 改之前生成基线（或用已有的）
   # …改代码…
   python tools/golden.py run mine
   python tools/golden.py compare base mine    # 末行必须是 GOLDEN_SAME
   ```

   浮点运算的顺序和结合方式不能动：`a*b + a*c` 不能写成 `a*(b+c)`，两次乘法不能合成一次，kernel 发射的顺序与实参不能变。这些都会改舍入。

2. **有意改数值的改动**（修 bug、改物理约定）：在提交信息和相关模块的 docstring 里写清楚改了什么约定、依据是哪次实测、影响哪些例子。判据以 Tidy3D 的输出为参考。

3. **不支持的功能一律 fail closed**：抛 `NotImplementedError`，不要静默降级。静默降级会让后面的偏差归因不到任何地方。

## 已知的非确定性

通量目标的伴随源会展开成几千个点偶极子，批量注入用 `atomicAdd`，float 求和次序不定，所以那类梯度**只能按容差比**（`golden.py` 里对应 case 已标 `~`）。模式目标的伴随源是 `ModeSource`，按位可复现。

## 代码风格

- Python ≥ 3.10，`ruff check openem tools --select F,E9` 必须干净。
- 代码里的注释和 docstring 一律用英文；README 中英双语，`docs/` 目前是中文。
- 模块边界见 `openem/__init__.py` 的 docstring：`scene/` 与 `nb/` 是唯一允许 import tidy3d 的地方，`device.py` 及少数几个模块才能 import cupy——`tests/test_no_tidy3d_leak.py` 与 `tests/test_no_cupy_leak.py` 锁着这两条。

## 提交

一个 commit 只做一件事，信息写清楚改了什么、为什么、怎么验的。
