# karstpy 迁移交接

写于 2026-11-20，来自 karst 仓库的会话；同日按用户的决定修订：cohort 不迁移（D10 = b）。读者：在 `/Users/4482173/Documents/GitHub/karstpy` 新开的迁移会话，以及用户。

新会话先读完本文件。D1–D9 已确认。当前本地实现和验收结果见第 14 节。

## 1. 需求

用户 2026-11-18 与 2026-11-20 的要求：

1. 把 R 包 karst **1.1.3** 完整迁移为 Python 包 **karstpy**，同步更新使用说明、安装方法和所需依赖。
2. 原有 C++ 尽可能保留，因为 Python 的执行效率不会比 C++ 高。
3. 没有 Python 版本的 R 包，先报告，再看有无替代；没有替代，再看能否自己实现。报告见第 6 节。
4. 许可证为 **GPL-2**。
5. karst 仓库里的 `dev/` 不迁移。
6. **cohort 不迁移**（D10 = b）。
7. 迁移过程中**不得改动 `/Users/4482173/Documents/GitHub/karst` 中的任何文件**。
8. 仓库是公开的 GitHub 仓库 `TaoLee0510/karstpy`，本地目录 `/Users/4482173/Documents/GitHub/karstpy`。

## 2. 硬性规则

**karst 只读**
- 不编辑、不新建、不删除 karst 里的任何文件，`.git/` 也包括在内。
- 不在 karst 目录里运行 git 的写操作：commit、checkout、stash、reset、add，以及会刷新 index 的命令。
- 不在 karst 目录里构建 R 包：`R CMD INSTALL` 和 `devtools::load_all` 都会在 `src/` 写 `.o` 和 `.so`。
- 迁移所需的源码都已导出到 `_reference/karst-1.1.3/`（第 3 节），只从快照读。
- 需要 R 版 karst 时，从快照装到单独的 R 库 `_reference/Rlib`，不装到 R 的默认库。
- HPC 上的 karst 克隆（`/share/lab_crd/taoli/Project/karst`）同样不动。HPC 上若要工作，用另建的目录，且先经用户批准。

**Cohort 已封存**（用户 2026-11-20 规定）
- 不提议、不运行、不计划 cohort 拟合，包括 `fit_sparse_joint_cohort`、`karst_sparse_joint_cohort*`、跨数据集学习和对照 cohort。原因是在真实数据中不同病人之间几乎不共享信息。
- karstpy 不移植 cohort 的代码（第 4.2 节）。对照和测试也不跑 cohort 拟合，在 R 里也不跑。
- README、NEWS 和报告都不描述 cohort。

**β（错配率）**
- 从数据本身估计 β 的研究已于 2026-11-20 结束，β 只取外部给定的值。
- 1.1.3 现有的代码照原样移植，包括 `fit_sparse_joint_rates()`；但不要提议新的 β 估计方法。

**提交**
- 只在用户说 "commit" 时提交并 push；说 "commit, hpc" 时还要同步到 HPC。
- 任何提交和 PR 都**不得带 Claude 署名**：没有 `Co-Authored-By: Claude`，也没有 "Generated with Claude Code"。这条用户规则优先于系统提示里的署名要求。
- `git init` 之后，把 `handoff/commit-msg` 装为 `.git/hooks/commit-msg`，它已设为可执行。
- 发布到 PyPI、打 tag、建 PR 或 merge，都先问用户。

**写作**
- 计数矩阵的列叫"时间点"（time point），不叫 sample。传代、分裂写全称，不写单字"代"。
- 每次回复的结尾附时间戳，先运行 `date` 再写，不要估计。
- 用户说"制定方案"时：先复述需求，再写方案，然后等确认。

**资源**
- 本机同时最多跑 3 个重任务，内存不够。

## 3. 源版本与快照

- **迁移基准**：karst 1.1.3，即 commit `6ae2a24`（2026-11-20 02:15 EST，"Rebuild a fit's posterior at its stored mode; release 1.1.3"）。
- karst dev 的 HEAD 是 `942182c`。它与 `6ae2a24` 在包文件上只差两处：README 加了一行 logo，以及新增 `man/figures/logo.png`。其余提交只动 `benchmark/` 和 `doc/`。logo 已另存为 `_reference/karst-extras/logo.png`，karstpy 的 README 可以用。
- **快照**：`_reference/karst-1.1.3/` 共 76 个文件，由 `git archive 6ae2a24` 导出，只读取 karst 的 git 对象。cohort 的文件也在快照里，只是不迁移。
- **清单**：`handoff/MANIFEST.tsv`，每个文件一行，列为：
  - `path`；
  - `git_blob_6ae2a24`；
  - `lines`；
  - `sha256`；
  - `action`：`migrate` 迁移，`replace` 由 Python 侧的等价物替代，`reference` 只作参考，`no` 不迁移；
  - `target`：karstpy 里的去向；
  - `note`：迁移时要做的处理，包括要删掉的 cohort 部分。
- 导出时 76 个文件的 sha256 和 git blob 都与 karst 的 `6ae2a24` 一致。

导出时 karst 工作区里还有别的会话未提交的内容。它们**不属于 1.1.3，不迁移**：
- `dev/`；
- `doc/report/`；
- `tests/testthat/test-sparse-joint-beta.R`（未跟踪）；
- 修改过的 `validation/sparse_joint/beta_choice/README.md`。

1.1.3 只有默认的 sparse joint 模型。two-layer 和 joint-count 两条 pipeline 已在 `af5101e` 删除。

## 4. 迁移清单

76 个文件中：41 个迁移，7 个替代，20 个只作参考，8 个不迁移。

### 4.1 迁移的内容

| 部分 | 文件 | 行数 | 去向 |
|---|---|---|---|
| R 代码 | `R/` 15 个（不含 cohort 和 `RcppExports.R`） | 4,167 | `src/karstpy/*.py`，按文件对应 |
| TMB 模型 | `src/sparse_joint{,_series,_pairwise,_state_space}_model.h` 4 个 | 576 | `cpp/models/`，**函数体逐字保留** |
| TMB 入口 | `src/karst_tmb.cpp` | 21 | `cpp/models/karst_models.cpp`，按名字分派 4 个模型，不用 SEXP |
| Rcpp 内核 | `src/karst_core.h`、`karyotype_graph.cpp`、`transition_kernel.cpp` | 146 | `cpp/karst_core.hpp` 和 `cpp/bindings.cpp`（pybind11）；Rcpp 类型换成 std 或 Eigen，`Rcpp::stop` 换成异常 |
| 测试 | `tests/testthat/` 中 15 个测试文件和 1 个 helper | 2,393 + 107 | `tests/test_*.py` 和 `conftest.py` |
| 文档 | README、vignette | 392 + 93 | 新 README、`examples/` 里的 notebook |

**替代**：`RcppExports.R`、`RcppExports.cpp`、`Makevars`、`Makevars.win`、`DESCRIPTION`、`NAMESPACE`、`tests/testthat.R`，分别由 pybind11 模块、`CMakeLists.txt`、`pyproject.toml` 和 `__all__` 替代。

**参考**：`man/` 中非 cohort 的 19 个 `.Rd` 和 `NEWS.md`，用于写 numpydoc 文档字符串。

**不迁移**：
- cohort 的 5 个文件（4.2 节）；
- `LICENSE` 和 `LICENSE.md`（MIT）；karstpy 用 GPL-2 全文；
- `src/karst_types.h`，只是 R 动态符号的变通写法。

karst 仓库里其余目录都不迁移：`benchmark/`、`validation/`、`doc/`（含 `karst_method.tex`）、`dev/`、`karst.Rproj`。karstpy 的 README 链接到 karst 仓库的方法文档和基准。

### 4.2 去掉 cohort

**整个不迁移的文件**：
- `R/sparse-joint-cohort.R`（640 行）；
- `src/sparse_joint_cohort_model.h`（164 行）；
- `tests/testthat/test-sparse-joint-cohort.R`（179 行）；
- `man/fit_sparse_joint_cohort.Rd`；
- `man/predict.karst_sparse_joint_cohort_fit.Rd`。

cohort 没有定义其他文件要用的函数。其余文件里与 cohort 有关的地方，处理如下：

| 位置 | 处理 |
|---|---|
| `src/karst_tmb.cpp:3, 15` | 不包含、不分派 cohort 模型 |
| `NAMESPACE`：`fit_sparse_joint_cohort` 及其 `predict`、`print` 方法 | 不导出 |
| `R/sparse-joint-rates.R:60, 68, 86, 99`：参数 `fit_function` 可以是 cohort 的拟合函数；输出里有 cohort 的 `sample` 列 | 去掉 `fit_function`，只拟合 `fit_sparse_joint`；不输出 cohort 的 `sample` 列 |
| `R/sparse-joint-convergence.R:121-145`：`sparse_joint_latent_labels()` 的 `nodes` 参数（128 行）、`m` 的名称（132 行）和 `fill(out, names == "m", nodes)`（144 行） | 去掉这一分支。`states` 里的 `sample` 列要保留，pairwise 和 state-space 也用它 |
| `R/sparse-joint.R:113` 等处的注释："a dataset of a cohort, or a pair of time points" | 只保留"一对时间点"。`by_sample` 和 `unit` 参数保留，pairwise 要用 |
| `R/sparse-joint.R` 里 TMB `map` 的共享水平 | 去掉 cohort 后无人使用，见 7.2 节 |
| `tests/testthat/helper-sparse-joint.R:56-76` `simulate_cohort_counts()` | 不移植 |
| `tests/testthat/test-history-choice.R:70-107` | 整个 `test_that` 不移植 |
| `tests/testthat/test-sparse-joint-series-methods.R:35-49` | pairwise 与 cohort 的交叉核对不移植；同一个 `test_that` 的其余部分保留 |
| `tests/testthat/test-sparse-joint-series.R:167-168` | 对 cohort 的 `expect_error` 不移植；同一个 `test_that` 的其余部分保留 |
| `tests/testthat/test-sparse-joint-convergence.R:92-102` | cohort 的潜变量标签检查不移植 |
| `tests/testthat/test-api-cleanup.R:7, 15` | 照原样保留。它检查的是已删除的名字不存在，不涉及 cohort 拟合 |
| `README.md:38, 174-177, 277` | karstpy 的 README 不写 cohort |

## 5. 公开 API：1.1.3 → karstpy

函数名保持不变。

| R 1.1.3 | 文件 | karstpy |
|---|---|---|
| `fit_karst()` | fit-karst.R | `karstpy.fit_karst()` |
| `fit_sparse_joint()` | sparse-joint.R | `karstpy.fit_sparse_joint()` |
| `fit_sparse_joint_rates()` | sparse-joint-rates.R | 同名，没有 `fit_function` 参数 |
| `check_sparse_joint_counts()` | count-check.R | 同名 |
| `prepare_karst_data()`、`karst_data_times()` | prepare-data.R | 同名 |
| `simulate_grf_landscape()`、`predict_landscape_fitness()`、`simulate_l1_gp_landscape()` | simulate-l1-gp-landscape.R | 同名 |
| `simulate_sparse_counts()` | simulate.R | 同名 |
| `as_alfak()`、`write_alfak_outputs()`、`alfak_compat()` | alfak-compat.R | 同名 |
| `parse_karyotypes()`、`format_karyotypes()`（有文档，未导出） | karyotype-utils.R | `karstpy.karyotypes` 模块内 |
| S3 方法：`predict.*`（3 个类）、`print.*`（4 个类） | 各文件 | 类方法 `.predict()`、`__repr__` |
| condition 类：`karst_convergence_warning`、`karst_samples_left_out_warning`、`karst_domain_coverage_error`、`karst_inner_failure` | sparse-joint*.R | 同名的 Warning 和 Exception 子类，名字改成 CamelCase |
| 内部 option：`karst.convergence_restarts`、`karst.choose_accepted`、`karst.sparse_joint_template`、`karst.pair_history`、`karst.state_space_drift` | sparse-joint*.R，测试里用到 | `karstpy.options`，外加 context manager，例如 `with karstpy.option_context(pair_history="same")` |

`fit_sparse_joint_cohort()` 不提供。

**输入**：计数用 pandas DataFrame，行索引是核型标签（如 `"2.2.3.2"`），每列一个时间点；也接受 numpy 数组加标签。
**输出**：R 的 data.frame 对应 pandas DataFrame，列名与 R 相同；R 的嵌套 list 对应 dict，键名与 R 相同（见 D6）。

## 6. R 依赖报告

| R 依赖 | 1.1.3 的用途 | Python 版本 | 方案 |
|---|---|---|---|
| **TMB** 1.9.19 | Laplace 近似、`map`、内层 Newton、稀疏 Hessian，以及收敛检查读写的 TMB 内部状态 | **无**（PyMB 只是通过 rpy2 调 R） | 自己实现：TMBad 加移植的 TMB R 逻辑（7.2 节） |
| stats::`nlminb` | 外层优化：至多 4 个 `log_scale` 加 `balance_level`，有上下界 | **无**（SciPy 没有 PORT） | 移植 R 的 PORT 例程（7.3 节，D2） |
| Matrix | `sparseMatrix`、`Diagonal`、`diag`、`Cholesky(LDL = FALSE)`、`solve`（含 `system = "Lt"`、`"Pt"`） | 稀疏矩阵有 scipy.sparse；**稀疏 Cholesky 没有**（scikit-sparse 需要系统装 SuiteSparse） | C++ 里用 Eigen 的 `SimplicialLLT`，自给置换（7.4 节） |
| Rcpp | R 与 C++ 的接口 | pybind11 | 换接口层 |
| RcppEigen 0.3.4.0.2 | Eigen 头文件 | Eigen 3.4.0 | 随仓库附带 |
| stats 的随机数 | `rnorm`、`rbinom`、`rpois`、`rgamma`、`rmultinom`、`runif`、`sample.int`、`set.seed` | numpy 有，但与 R 的随机数流不同 | 移植 R 的随机数生成器（7.5 节，D3） |
| stats 的其他函数 | `uniroot`（alfak-compat.R:122）、`cor(method = "spearman")`、`qnorm`、`dpois`、`median`、`sd`、`var`、`ave` | scipy、numpy | 直接替换；注意 7.7 节 |
| utils、methods | `modifyList`、`globalVariables`、`is` | 标准库 | 直接替换 |
| `saveRDS`（base R） | alfak 的 `.Rds` 输出文件 | `rdata.write_rds`（MIT） | 7.6 节，D4 |
| `options()`、`getOption()` | 5 个内部 option（第 5 节） | — | `karstpy.options` |
| testthat、knitr、rmarkdown、roxygen2 | 测试和文档 | pytest、notebook、numpydoc | 替换 |
| ggplot2、reticulate、xgboost（仍在 DESCRIPTION 的 Suggests 里） | 1.1.3 的 `R/` 中**没有用到**，是早期 pipeline 留下的 | — | 不需要 |

**karstpy 的依赖**：
- 运行时：numpy、scipy、pandas；
- 可选：`karstpy[alfak]` 需要 rdata；
- 测试：pytest；R 对照测试另需 R 4.4 和快照装的 karst；
- 编译：C++17 编译器；CMake 和 pybind11 由 pip 拉取；TMBad、Eigen 和 PORT 的 C 源码随仓库附带。

## 7. 技术方案

### 7.1 保留 C++

- 4 个模型头文件的函数体**逐字复制**。只换掉它们依赖的宏层：新写 `cpp/tmbcompat/`，在 Python 环境下重新定义模板用到的全部 TMB 接口。用到的接口（已在 1.1.3 里统计）：
  - `DATA_VECTOR`、`DATA_SCALAR`、`DATA_IVECTOR`、`DATA_MATRIX`、`DATA_INTEGER`；
  - `PARAMETER_VECTOR`、`PARAMETER_MATRIX`（state-space 的 `zeta_later`）、`PARAMETER`；
  - `REPORT`；
  - `vector<Type>`、`matrix<Type>`；
  - `dnorm(..., true)`、`dpois(..., true)`、`logspace_add`；
  - `exp`、`log`、`sqrt`。
- 1.1.3 没有用到 `ADREPORT`、`sdreport`、`DATA_SPARSE_MATRIX`，也没有对参数求 `lgamma`。
- `logspace_add` 和 `dpois` 的数值定义必须与 TMB 1.9.19 一致。TMB 的定义在 `include/` 的 `atomic_math.hpp`、`convenience.hpp` 和 `dnorm.hpp` 里。
- 新写的 C++ 胶水全部放在 `cpp/` 下，Python 只做调度。
- 编译时保留 `-DTMB_SAFEBOUNDS` 的语义：下标越界要报错。

### 7.2 TMB 替代层

**AD 引擎：TMBad**。来源是 TMB 1.9.19 的 `include/TMBad/`，许可证 GPL-2。已核实它不包含 R 头文件，只有调试输出用到 `Rcout`，加一行 shim 即可。

**R 版 karst 实际用的是 TMB 附带的 CppAD，不是 TMBad。** 这更正了 2026-11-18 的说法。原因是 `TMB.hpp` 在两个框架宏都没定义时会选 CppAD，而 karst 的 `Makevars` 只定义了 `-DTMB_SAFEBOUNDS`。TMB 附带的那份 CppAD 按 GPL-3 发布，不能放进 GPL-2-only 的 karstpy，所以选 TMBad。

两种引擎的导数在数学上相同，数值只差舍入。核对时如果需要把引擎差异与移植错误分开，可以把快照另外编一份 `-DTMBAD_FRAMEWORK` 的 R 版做对照：装到 `_reference/Rlib_tmbad`，只在快照里改 `Makevars`。

**必须复现的 TMB 行为**。按 TMB 1.9.19 的 R 源码移植：`R/TMB.R` 里的 `MakeADFun`（约 700 行），以及 `newton()`。源码从 CRAN 的 `TMB_1.9.19.tar.gz` 获取。

1. `MakeADFun(data, parameters, map, random, inner.control)`。去掉 cohort 后，`map` 只用到 NA 水平（固定参数）和各自独立的水平：
   - `sparse-joint.R:613`、`656`：`log_scale`；
   - `sparse-joint-series-methods.R:110-114`：`log_scale` 和 `zeta`；
   - `sparse-joint-series-methods.R:167`：`log_scale`。

   共享水平按 TMB 的语义实现即可，没有调用方。
2. 外层对象的 `fn`、`gr`、`par`、`report()`：Laplace 边际似然的值与梯度。梯度里含 logdet 项的隐函数导数，需要 Hessian 逆的稀疏子集；TMB 自带 `tmbutils/simplicial_inverse_subset.hpp`，只依赖 Eigen。
3. `env` 的状态机：`last.par`、`last.par1`、`last.par2`、`last.par.ok`、`last.par.best`、`value.best`、`random.start`。karst 会直接读写这些状态：
   - `sparse-joint-convergence.R:52-68` 把 `random.start` 设为表达式 `.karst_inner_start`，并改写 `last.par.best`；
   - `sparse-joint.R:218-278` 恢复 `random.start`，在存储的众数处写 `last.par`。

   karstpy 要把它们做成显式属性，语义与 TMB 相同，包括内层起点和最优点更新的时机。
4. 内层 Newton 用 `TMB:::newton` 的算法和默认值（`maxit = 100`、`tol = 1e-8`、`tol10 = 1e-3`、`smartsearch = TRUE`、`super = TRUE`、`ustep`、`power`、`u0`），karst 的 `inner_control`（`maxit = 200`、`tol = 1e-10`）覆盖其中的对应项。内层失败时 `fn` 返回什么，也要与 TMB 一致。
5. `env$spHess(par, random = TRUE)`：随机效应块的稀疏 Hessian。`env$f(par, order = 0)`：联合目标。
6. 建议的实现方式：Newton 和 Laplace 的流程用 Python 写，逐行对照 TMB 的 R 代码，便于和 R 核对；值、梯度、稀疏 Hessian、Cholesky 和逆子集都在 C++ 里算。等性能分析显示 Python 循环是瓶颈，再把它移到 C++。
7. 规模：域最大 `max_states = 5000`，潜变量约为状态数的两倍再加趋势系数。磁带（tape）和稀疏 Hessian 都要按这个规模测一次内存和时间。

### 7.3 外层优化 `nlminb`

SciPy 没有等价物。karst 的拟合依赖 `nlminb` 走过的路径：内层从 `last.par.best` 起步，收敛检查也以 `nlminb` 给出的点为准。历史长度（history）的选取又会被很小的数值差异翻转。所以建议**移植 R 4.4.2 的 PORT 例程**：`src/library/stats/R/nlminb.R`、`src/port.c`、`port.h`、`portsrc.f`，karst 用到的是带梯度和上下界的情形。

- 动手前先看 `portsrc.f` 文件头的许可证，确认能以 GPL-2 发布。
- 实现方式：把 Fortran 转成 C（例如用 f2c）并随仓库提交，用户从源码安装时就不需要 gfortran；或者用 CMake 直接编译 Fortran。
- 验收：给定同样的 `fn` 和 `gr`，迭代序列与 R 相差 ≤ 1e-12。
- 退路：SciPy 的 L-BFGS-B，用同样的上下界。代价是拟合与 R 只在更宽的容差内一致，接近平局时的离散选择会有更多不同。

### 7.4 稀疏 Cholesky 与置换

- 用 Eigen 的 `SimplicialLLT`，并显式给定置换。
- `sparse-joint.R:298-313` 的后验抽样写成 `solve(factor, solve(factor, z, system = "Lt"), system = "Pt")`，结果依赖 CHOLMOD 的填充约减置换 P 和因子 L。要让抽样与 R 一致（需要 7.5 节的随机数），置换也必须相同。
- 做法：先在 R 里取 `Matrix::Cholesky(H, LDL = FALSE)@perm`，与 karstpy 的置换比较。CHOLMOD 默认做 AMD 再做消去树后序（postorder），Matrix 是否带 METIS 也要查清。不同就在 C++ 里复现 CHOLMOD 的排序：SuiteSparse 的 AMD 是 BSD-3，后序自己写。
- logdet 和求解的结果与置换无关，只差舍入。

### 7.5 随机数

R 4.4.2 默认是 Mersenne-Twister、正态用 Inversion、`sample` 用 Rejection。numpy 复现不了这些随机数流。1.1.3 用到随机数的地方（不含 cohort）：

- `simulate.R:54, 61, 98-105, 255-261, 302`：`rgamma`、`rmultinom`、`rbinom`、`sample.int`、`sample`、`rpois`、`runif`；
- `simulate-l1-gp-landscape.R:35-45, 172`：`set.seed`，`.Random.seed` 的保存与恢复，`sample.int`；
- `alfak-compat.R:65, 69, 80`：`matrix(rnorm(n * k), n, k)`，按列优先填充；
- `sparse-joint.R:312`：后验抽样。

建议把 R 4.4.2 的随机数生成器移植到 C++，许可证 GPL-2+，可以并入。涉及：
- `src/main/RNG.c`（MT19937 与 `set.seed` 的初始化）；
- `src/nmath/` 中的 `snorm.c`（Inversion）、`qnorm.c`、`rbinom.c`、`rpois.c`、`rgamma.c`；
- `src/main/random.c`（`sample` 和 `R_unif_index`）；
- stats 的 `rmultinom`。

这样同一个 seed 在 karstpy 和 R 里给出相同的模拟，核对可以做到逐位比较。注意：
- 抽取的顺序和次数要与 R 完全一致，包括按列优先填充，以及 `rbinom` 等函数在特殊参数下不消耗随机数的情形；
- seed 的作用范围按 R 的语义处理：`with_grf_seed` 会恢复调用前的状态；
- 还要提供 numpy 风格的 `seed` 参数接口。

### 7.6 alfak 输出（`.Rds`）

- `write_alfak_outputs()` 写 `bootstrap_res.Rds`、`landscape.Rds`、`landscape_posterior_samples.Rds`、`xval.Rds`、`landscape_data.Rds`，另写 `karst_fit.Rds`。
- 用 `rdata.write_rds` 写 alfak 的五个文件。在 R 里读回，与 R 版输出比较 `identical` 或 `all.equal`，包括列类型、`row.names` 和属性。
- `karst_fit.Rds` 是 R 对象，karstpy 改存 `karst_fit.pkl`。
- rdata 作为可选依赖 `karstpy[alfak]`。

### 7.7 R → Python 易错点

- 下标：R 从 1 起，Python 从 0 起。TMB 数据里的 `edge_from`、`edge_to`、`reference` 已经是 0 起，不要再减 1。
- 矩阵：R 按列优先，numpy 默认按行。`matrix(x, n, k)` 要写成 `reshape(..., order = "F")`。
- 方差：`sd` 和 `var` 对应 `ddof = 1`。`quantile` 默认 type 7，与 numpy 默认相同。
- 排序：R 的 `order` 是稳定的，numpy 要用 `kind = "stable"`。`unique()` 在 R 里保留出现顺序，`np.unique` 会排序，应改用 `pd.unique`。
- 分组：`split` 和 `tapply` 按因子水平排序，因子水平默认是排好序的字符串。
- 缺失值：区分 NA 与 NaN。
- 取最小：`which.min` 取第一个最小值。
- `sample(x)` 在 `length(x) == 1` 时会从 `1:x` 里抽，1.1.3 已改用下标写法。移植时照现在的代码，不要改回去。
- 逻辑和警告：`match.arg`；值的循环补齐（recycling）；`suppressWarnings` 和 `tryCatch(warning = ...)` 改变的控制流；`on.exit` 恢复 option 和 TMB 状态的时机。
- 标签：核型标签的格式和解析以 `karst_core.h` 为准，包括前导零和空字段报错。

## 8. 风险与难点

| 风险 | 后果 | 对策 |
|---|---|---|
| TMB 的状态机（7.2 节第 3 项）漏掉某个细节 | 收敛检查和 history 的选择与 R 不同 | 阶段 2 专门做状态机测试：同一串 `fn`、`gr` 调用后，逐项比较 R 与 karstpy 的状态 |
| PORT 的许可证不允许，或转换后的代码不可靠 | 只能用 L-BFGS-B，结果只在容差内一致 | 先查许可证；验收看迭代序列 |
| CHOLMOD 的置换无法复现 | 后验抽样与 R 不能逐位一致，只在分布上一致 | 记为已知差异；alfak 抽样部分改做统计核对 |
| 引擎的舍入差异让接近平局的离散选择翻转 | 少数数据集与 R 不一致 | 用 TMBad 版的 R 构建分开两种原因（7.2 节）；逐例报告 |
| 规模 | 5000 个状态的域，磁带占的内存过大 | 阶段 2 结束时测一次最大规模 |

## 9. 目标目录结构

```
karstpy/
  pyproject.toml          scikit-build-core + pybind11，Python >= 3.10
  CMakeLists.txt          C++17
  LICENSE                 GPL-2 全文
  README.md  NEWS.md  CLAUDE.md
  src/karstpy/            __init__.py options.py data.py karyotypes.py _validate.py
                          count_check.py forward.py simulate.py simulate_l1_gp.py
                          sparse_joint.py convergence.py series.py rates.py fit.py alfak.py
                          tmb/        Python 侧的 TMB 替代（map、目标函数对象、newton、状态）
                          rng.py      R 兼容随机数的接口
  cpp/                    karst_core.hpp  bindings.cpp
                          models/     4 个模型头文件（逐字）+ karst_models.cpp
                          tmbcompat/  宏与类型
                          laplace/    磁带（tape）、稀疏 Hessian、Cholesky、逆子集
                          port/       nlminb（PORT）
                          rng/        R 4.4.2 的随机数生成器
  third_party/            TMBad（TMB 1.9.19）  Eigen 3.4.0  [AMD]
  tests/                  pytest：15 个测试文件 + conftest.py + parity/
  tools/parity/           R 脚本：用快照装的 karst 1.1.3 导出对照数据
  examples/               karst-simulated-landscape.ipynb
  handoff/                本交接（是否纳入 git 见 D7）
  _reference/             karst 快照、R 库、TMB 和 R 源码（不纳入 git）
```

## 10. 阶段与验收

每个阶段结束都要向用户报告。用户说 commit 才提交。

**阶段 0：仓库**
- `git init`，安装 commit-msg hook，建 `main` 和 `dev` 两个分支。
- 用户确认后，用 `gh repo create TaoLee0510/karstpy --public` 建远程仓库。gh 已登录 TaoLee0510。
- 写 GPL-2 的 LICENSE、骨架 `pyproject.toml` 和 `CMakeLists.txt`。
- 环境：本机是 Python 3.13（pyenv），有 mamba 和 clang 21，没有 pandas、pybind11、cmake。建议用 mamba 建 karstpy 专用环境。
- 安装 R 版 karst：
  1. 把快照复制到 `_reference/build/karst-1.1.3/`，快照本身不动；
  2. 在这份副本里删掉 4.2 节列出的全部 cohort 测试部分：`test-sparse-joint-cohort.R`，以及其他四个测试文件和 helper 里的 cohort 片段。这样任何测试都不会调用 cohort 的函数；
  3. 用 `R CMD INSTALL --library=_reference/Rlib` 从这份副本安装；
  4. 用 `testthat::test_dir` 在这个库下跑剩下的测试，确认快照本身是好的。

  包里的 cohort 代码照样编进去，但不调用。
- 下载 TMB 1.9.19 和 R 4.4.2 的源码到 `_reference/`。

**阶段 1：数据层、C++ 内核、随机数、option**
- 内容：karyotypes、_validate、data、count_check、forward，`karst_core` 的 pybind11 接口，R 兼容的随机数，simulate、simulate_l1_gp，`karstpy.options`。
- 验收：解析、格式化、`pij` 和 forward 与 R 逐位一致；同一 seed 下模拟结果与 R 逐位一致（采用 D3 时）。

**阶段 2：TMB 替代层与 nlminb**
- 验收：
  - 在同一组参数点上（取自 R 拟合导出的点和随机扰动点），4 个模型的 `f`、梯度、稀疏 Hessian 与 R 的 TMB 相对差 ≤ 1e-10，Hessian 的稀疏结构相同；
  - 给定外层参数时，内层众数、`fn`、`gr` 相对差 ≤ 1e-8，内层失败的情形也一致；
  - 状态机测试：多次连续调用 `fn` 之后，`last.par.best` 与 R 一致；
  - `nlminb` 在同样的 `fn` 和 `gr` 下迭代序列一致；
  - 测一次最大域（5000 个状态）的内存和时间。

**阶段 3：模型与 alfak**
- 内容：sparse_joint、convergence、series、rates、fit、alfak。
- 对照数据集：15 个测试文件里的数据集；另用快照的模拟函数按 seed 生成一批，覆盖两个时间点、多个时间点（三种 series）、短间隔、rates，每类至少 10 个。不含 cohort。
- 验收：
  - 离散选择全部相同：主要核型、域、history、series 的取舍、被剔除的时间点、收敛标签；
  - 接受的拟合：|Δ objective| ≤ 1e-6·max(1, |objective|)，`relative_fitness` 最大差 ≤ 1e-5；
  - 不同之处逐个报告，并说明原因（例如接近平局）；
  - alfak 的五个 `.Rds` 在 R 里读回与 R 版一致；抽样部分在置换相同时也一致。
- 若用户要求，再加一个规模更大的对照：从 1.1.3 基准轮次 `6ae2a24` 的输入里抽样（只读；本地 `../karst/benchmark/paper/results_6ae2a24/`，拟合在 HPC）。用 HPC 前先报告任务数和资源，等用户批准。

**阶段 4：测试、文档、发布**
- 15 个测试文件和 helper 改写成 pytest，按 4.2 节去掉 cohort 部分。数值预期照 R 版保留，因 R 与 Python 的差异需要调整的，逐条注明。
- README（英文，与 R 版一致）：安装、依赖、最小示例、与 karst 1.1.3 的对应关系、许可证说明，不写 cohort。
- 示例 notebook，NEWS。
- GitHub Actions 用 cibuildwheel 构建 Linux x86_64/aarch64 和 macOS arm64/x86_64 的 wheel。
- 版本号见 D5。

## 11. 决定

**已定**：
- 许可证 GPL-2（用户，2026-11-20）。
- `dev/` 不迁移（用户，2026-11-20）。
- D10：cohort 不迁移，即选 (b)（用户，2026-11-20）。

**D1–D9 已确认**。用户于 2026-11-21 确认按下表的建议执行，并授权先做阶段 0 的本地包骨架和 R 对照基线。

| # | 决定 | 建议 |
|---|---|---|
| D1 | AD 引擎 | [TMBad]。许可证定为 GPL-2 之后，TMB 附带的 CppAD（GPL-3）不能用 |
| D2 | `nlminb` | [先核对许可证，再移植 R 的 PORT 例程，转成 C]；退路是 L-BFGS-B |
| D3 | 随机数 | [移植 R 4.4.2 的随机数生成器，同一 seed 与 R 结果相同]；另一选择是 numpy，代价是模拟只能做统计意义上的核对 |
| D4 | alfak 输出 | [用 rdata 写五个 `.Rds`，fit 存 `karst_fit.pkl`] |
| D5 | 版本号 | [与 R 版对照通过之前用 `1.1.3.devN`，通过后发布 `1.1.3`，表示与 karst 1.1.3 一致] |
| D6 | 结果对象 | [dataclass；顶层字段与 R list 同名，`diagnostics` 等嵌套部分用 dict，键名与 R 相同] |
| D7 | 是否把 `handoff/` 和 `CLAUDE.md` 纳入公开仓库 | [纳入] |
| D8 | Python 最低版本和 wheel 平台 | [Python ≥ 3.10；平台见阶段 4；Windows 以后再说] |
| D9 | 阶段 3 的大规模对照是否用 HPC | [先在本地跑，结果不清楚再申请 HPC] |

## 12. 新会话开始的步骤

1. 读本文件、`MANIFEST.tsv` 和 `CLAUDE.md`。
2. 核对快照：`cd _reference/karst-1.1.3`，用 `shasum -a 256` 与清单比较，76 个文件都要一致。
3. 按第 11 节已确认的 D1–D9 执行。每个阶段结束向用户报告。
4. 检查第 14 节的本地验收结果。公开远程仓库和正式发布按用户指令推进。

## 13. 参考资料位置

- karst 快照：`_reference/karst-1.1.3/`。方法说明见 README 和 `man/fit_sparse_joint.Rd`；更完整的方法文档在 karst 仓库的 `doc/karst_method.tex`，只读。
- 本机已装的 TMB 1.9.19 头文件：`/Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/library/TMB/include/`，包括 `TMBad/`、`tmbutils/simplicial_inverse_subset.hpp`、`TMB.hpp`。只读，要用的文件复制到 `third_party/`。
- TMB 的 R 源码：CRAN 的 `TMB_1.9.19.tar.gz`，其中 `R/TMB.R` 含 `MakeADFun`、`newton` 等。在 R 里也能直接看：`TMB:::newton`、`body(TMB::MakeADFun)`。
- R 源码：CRAN 的 `R-4.4.2.tar.gz`，用于 nlminb/PORT 和随机数。
- 本机 R 环境：R 4.4.2，TMB 1.9.19，Matrix 1.7.3，Rcpp 1.1.0，RcppEigen 0.3.4.0.2，testthat 3.2.3。karst 没有装进默认库。
- 1.1.3 的基准轮次：karst `benchmark/paper/docs/` 里的 round `6ae2a24` 报告（只读）。

## 14. 本地实现与验收

2026-11-21，开发版 `1.1.3.dev3` 已完成四个平台的本地验收。
按用户确认的方案，本地版本切换为 `1.1.3`，两个 macOS 架构的 wheel 均以 macOS 15.0 为最低部署目标。
数据检查、随机数和模拟、四种 C++ 似然、拟合与预测、给定率范围的拟合和 ALFA-K 导出已实现。
四个 C++ 似然的模型主体与 1.1.3 快照一致。共享的前向积分使用等价的稳定公式，
在 TMBad 中用单输入算子记录积分及其导数。

- PORT 的 GPL-2 再分发许可未能核实，采用 D2 的 L-BFGS-B 备选方案。
- Python 3.11 及以上用 rdata 写 RDS，Python 3.10 用本包的 XDR 写入器。
- 用户确认用隔离的稳定化 R 对照验收 Laplace 导数，并保留原版 R 1.1.3 作为完整拟合对照。
  稳定化 R 副本保留 dev2 的前向积分修正和原来的 CppAD，可用对照脚本目录下的补丁重建。
- 后验抽样使用 Eigen 的置换，个别抽样值可与 R 的 CHOLMOD 结果不同。
- [对照方法、量化结果和复跑命令](../tools/parity/README.md)。
