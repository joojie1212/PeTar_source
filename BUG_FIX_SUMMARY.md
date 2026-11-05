# Bug 最小修复汇总

修复日期：2026-11-04  
依据：`bugs/BUG_TABLE.md` 中 B01–B40，以及后续报告的近共振异常 B41、SIMD 静态 TLS 启动崩溃 B42，并按当前工作区源码重新核对。  
状态统计：**已修复 38 项；当前版本未复现/路径已串行 2 项；部分修复 2 项。**

状态说明：

- **已修复**：已实施最小代码修改，并至少通过编译或对应回归测试。
- **当前未触发**：报告所述并行路径在当前源码中不存在或为串行，未为不存在的路径引入额外同步。
- **部分修复**：已修复能明确定位的实例，但原报告范围仍有未完成部分。

## 修复表

| 编号 | 模块 | 状态 | 实际最小修改 | 验证与注意事项 |
|---|---|---|---|---|
| B01 | SDAR / petar | 当前未触发 | 当前 `integrateGroupsOneStep()` 按组串行执行，没有截图所述组级并行读可变扰动源的路径，因此未加入会改变算法语义的锁。 | 已静态核对当前调度；若以后恢复 `PETAR_PAR_AR`，仍需事件屏障和扰动源版本化。 |
| B02 | petar | 已修复 | `bse`、`bseEmp`、`mobse` 三套 `trdot.f` 的紧致天体分支改用已定义的 `age*10`，不再读取未初始化 `dt`。 | BSE 集成回归与完整构建通过。 |
| B03 | petar | 已修复 | `getTimeStepBinary()` 将 GR 和质量转移候选先各自换算为 N-body 时间，再比较并直接返回，去掉末尾第二次 `/tscale`。 | 回归覆盖 `tscale=0.5,1,2`，物理时间步一致。 |
| B04 | SDAR | 当前未触发 | 当前邻居/校正路径对 `initial_step_flag` 的读写为串行循环，未发现截图中的共享并行写路径，未额外修改。 | 已静态核对；未来并行化该循环时仍需线程局部归并。 |
| B05 | petar | 已修复 | 两个 `bse_manager.merge()` 入口纳入与其他 BSE Fortran 调用相同的未命名 OpenMP critical 区。 | 默认 OpenMP 构建通过；事件先后顺序仍由调度决定。 |
| B06 | petar | 已修复 | BH 自旋投影方位角由无符号 `acos` 改为按 `S1·(S2×L)` 约定的 `atan2`；Python 进动输出沿用初始相位分支符号。 | 正负分支回归通过；负初始相位保持为负分支。 |
| B07 | petar | 已修复 | ASCII 粒子行增加可选 `#BHMERGER_V1 time_merger` 后缀；读取器兼容旧行，旧格式将 deadline 初始化为 `-1`。 | 新格式双行往返和旧格式读取回归通过。 |
| B08 | petar | 已修复 | `modify_branch` 初始化为 `{0,0}`，避免 BSE 关闭或叶分支未赋值时读取垃圾值。 | 默认及宏组合编译通过。 |
| B09 | petar | 已修复 | SIMD 邻居计数改为用原始双精度坐标做相对距离判定；组合力路径保留 SIMD 力计算，仅用双精度复核邻居数。 | `0` 与 `1e8` 整体平移回归均得到相同邻居数；代价是邻居计数多一层标量复核。 |
| B10 | FDPS | 已修复 | `getEpjFromId()` 直接遍历已压紧、只含 EP 的 `epj_sorted_`，不再用同下标的混合 `tp_glb_` SP 标志过滤。 | EP/SP 交错映射回归通过；尚未做多 rank 专项压力测试。 |
| B11 | petar | 已修复 | `--id-offset` 从 `atoi` 改为 `strtoll`，并检查非法尾字符、`ERANGE`、负值及 `PS::S64` 范围。 | 默认构建通过，64 位解析路径不再经 32 位截断。 |
| B12 | SDAR | 已修复 | `HERMITE_DEBUG` 检查数组改为按 `groups.getSize()` 分配的 `std::vector`，索引域与组槽位一致。 | `HERMITE_DEBUG` 宏分支语法编译和 hard-debug 完整构建通过。 |
| B13 | petar | 已修复 | 所有活动 AVX2/AVX-512、单/双精度 EP/SP 预取路径在末轮将 `j+1` 夹到当前 `j`。 | `USE_QUAD+AVX_PRELOAD` 编译通过；尚未做真实 `NJMAX` 力值压力对照。 |
| B14 | petar | 已修复 | `evolveBinary()` 与 `merge()` 的 `menv/renv` 输出数组零初始化，覆盖 Fortran 提前返回未赋值情形。 | BSE 集成回归和完整构建通过。 |
| B15 | FDPS（已定位部分） | 部分修复 | 已定位的 `calc_force__core__walk_tree/dispatch` 改为 OpenMP reduction 后在循环外统一累计。 | 已消除该处共享计时写；原截图称“若干”但未列出的 petar/SDAR profiling 点仍需 Helgrind 清单才能逐项确认。 |
| B16 | petar | 已修复 | 孤立粒子并行演化的动能变化使用 reduction，循环后一次更新四个能量累计字段。 | OpenMP 完整构建及集成 sanitizer 回归通过。 |
| B17 | petar | 已修复 | 逃逸粒子能量损失使用 reduction，循环后统一更新参考能量与累计修正。 | 4096 粒子、1/4 线程结果一致。 |
| B18 | petar / SDAR | 已修复 | 递归父节点仅在左右分支均返回 destroy 时传播整组 destroy；单侧销毁降为 merge/change。 | 相关模板及默认构建通过；未做长时多层事件序列复现。 |
| B19 | petar | 已修复 | BH remnant 接口质量改用 `getMass(star)` 的 N-body 质量，不再把 BSE 太阳质量与 N-body 位置速度混用。 | BH/Python 桥接构建和异常回归通过。 |
| B20 | petar | 已修复 | `evolveBinary()` 在任何同步/预演化前保存两颗星入口质量，所有返回路径的 `dm` 都相对该入口快照计算。 | 异步年龄 BSE 回归核对 `dm == mt_final-mt_entry`。 |
| B21 | petar | 已修复 | `RECORD_CM_IN_HEADER` 下，逃逸判定和删除能量记账都先恢复 `stat.pcm.vel`。 | 1/4 线程回归按恢复后的速度得到相同能量。 |
| B22 | petar | 已修复 | `EnergyAndMomentum::writeBinary()` 从 `this` 起始地址写完整对象，不再从中间字段 `&ekin` 越界读取。 | 非首字段、`ekin`、角动量字段二进制往返通过。 |
| B23 | SDAR | 已修复 | destroy 结算同步更新 `de_sd_change_interrupt_`，与 slowdown 总累计保持同一事件分类。 | 默认与 SDAR 宏分支编译通过。 |
| B24 | petar / SDAR | 已修复 | BSE 后轨道锥曲线类型变号时，用当前分离距离和径向速度符号重算异常角；不相容的分离距离显式拒绝，并重算全部轨道缓存。 | 默认构建通过；近抛物线和双向转换仍建议用物理场景做数值连续性验收。 |
| B25 | SDAR | 已修复 | 两个 `bin_host[ki]` 循环把 `ki` 边界判断移到数组访问之前。 | 3–6 粒子树在 ASan/UBSan 下通过。 |
| B26 | petar | 已修复 | `INTEGRATED_CUTOFF_FUNCTION` 下构造、`setR` 和二进制读取后均建立 `q_=r_in/r_out`；读取要求恰好两个值。 | 宏分支编译、势函数有限性和二进制往返通过。 |
| B27 | petar | 已修复 | `AR_CHANGEOVER` 分支的拼写 `invr` 改为已有局部变量 `inv_r`。 | `AR_CHANGEOVER` 宏分支语法编译通过。 |
| B28 | petar | 已修复 | 四个 remnant 调用及 inspiral 字典先检查 NULL/类型，再校验数组形状、长度、字节序、连续性和有限值；异常直接失败，不提交部分结果。 | 注入四类 Python 异常均不崩溃，测试通过。 |
| B29 | petar | 已修复 | 为临时 float、空参数 tuple、NumPy item/array 和结果补齐引用释放；字典插入后释放本地引用。 | Python 循环回归保留内存增长约 10 KB，低于 64 KB 门限；沙箱内关闭 LeakSanitizer。 |
| B30 | petar | 已修复 | `DEBUG Dynamic_merge` 输出置于 `AR_DEBUG_PRINT` 条件下。 | 正常构建不再无条件输出该诊断。 |
| B31 | petar | 已修复 | `PtclHard::print()` 的 routing 字段统一写入调用者传入的 `_fout`。 | `ostringstream` 捕获回归通过。 |
| B32 | petar / SDAR | 已修复 | `acos` 输入统一夹到 `[-1,1]`，非有限/退化几何显式返回不稳定或最大不稳定度。 | sanitizer 回归及默认构建通过。 |
| B33 | SDAR | 已修复 | `copyDataIter()` 根据前一子树实际复制节点数推进下一写入位置，避免左右子树覆盖。 | 6 叶、5 内部节点拓扑及指针归属回归通过。 |
| B34 | SDAR | 已修复 | `AR::Information` 用 `sizeof(ds)` 读写浮点步长，并加入 `ARI2` 格式标记；旧损坏格式显式拒绝。 | 非整数 `ds` 往返通过；旧 AR debug dump 不再静默误读。 |
| B35 | SDAR | 已修复 | `calcNextDtLimit()` 检查负值、非有限和 64 位转换溢出；`bitmap==0` 返回 `dt_min`，避免无限右移。 | `0`、`nextafter(dt_min,0)`、`dt_min` 和块边界回归通过。 |
| B36 | SDAR | 已修复 | `List::operator=` 改用现有单参数 `reserveMem()`，增加自赋值保护，并保留 copy/local 模式语义。 | copy/local 赋值和自赋值回归通过。 |
| B37 | SDAR | 已修复 | `ParticleH4::print()` 逐分量打印 `acc0/acc1`，且 `acc1` 不再误用 `acc0`。 | 不同向量值输出回归通过。 |
| B38 | SDAR | 已修复 | 修正 `_bini`、`calInvR`、错误累积变量/返回名、函数大小写，以及 slowdown 数组错误成员访问。 | `AR_TIME_FUNCTION_MULTI_R`、`AR_SLOWDOWN_MASSRATIO` 分别及组合编译通过。 |
| B39 | SDAR | 已修复 | `SymplecticStep` 只保存 `SYS2 + order + type`，不再持久化裸指针；读取先校验、`clear()` 旧数组再重建。 | 2/4/6/-6/-8 阶各重复读取 50 次通过；旧裸指针格式被拒绝。 |
| B40 | petar | 部分修复 | 去掉进动输出相位分支的 `np.random.choice`，改为保持初始有向相位；未把 BSE Fortran 与所有 BH 随机源的完整状态加入常规检查点。 | 相位分支已确定；连续运行与中途重启的完整随机序列仍不能保证逐事件一致。 |
| B41 | petar | 已修复 | 当完整 `inspiral_precav` 因近共振的非有限进动周期而失败时，C++ 桥接层先仅请求演化后的守恒量和三个 `deltachi` 根，避开随机进动相位采样；仅在两个物理解根按尺度相关容差确认合并后，调用库内已有 `deltachiresonance` 取共振极限，再由 `conserved_to_angles` 重建末端自旋角并保持初始相位分支。Python 源码未改。 | 注入同型 `OverflowError` 后成功恢复有限角度；根未合并时不会误用共振解，仍回退到遭遇时刻自旋。正常非共振路径不变。 |
| B42 | petar | 已修复 | x86 SIMD 力核的多个大型 `static thread_local PhantomGrapeQuad` 改为线程局部 `unique_ptr` 首次使用时堆分配，并为缓冲类型显式增加 64 字节对齐；暗物质 SIMD 力核同步处理。保留每线程独立实例和 SIMD 算法。 | 发布程序 `.tbss` 由约 18.1 MiB 降至约 5.1 KiB；默认 8 MiB 栈下 `petar -h` 返回 0。回归检查堆缓冲对齐、SIMD/标量力与邻居数一致，并限制测试程序 `.tbss <= 1 MiB`。 |

## 验证结果

- `make -C PeTar_new_version -j2 all`：通过；生成默认 MPI/OpenMP/AVX2/BSE 主程序、hard-debug 和格式转换工具。
- `bash PeTar_new_version/test/bugfix/run.sh`：通过：
  - 核心结构/序列化 ASan+UBSan 回归；
  - PeTar/FDPS/BSE/OpenMP 集成 ASan+UBSan 回归；
  - SIMD 工作区 64 字节对齐、堆生命周期、标量数值对照及静态 TLS 大小回归；
  - Python remnant 异常、近共振恢复/非共振拒绝恢复与引用保留量回归；
  - Python 进动相位符号回归。
- 可选宏语法编译均通过：`AR_TIME_FUNCTION_MULTI_R`、`AR_SLOWDOWN_MASSRATIO`（含二者组合）、`AR_CHANGEOVER`、`INTEGRATED_CUTOFF_FUNCTION`、`HERMITE_DEBUG`、`USE_QUAD+AVX_PRELOAD`。
- 现有 Peters 合并时间测试通过：最大采样相对误差 `0.00113518`；渐近切换跳变 `0.020585`。该功能及其测试属于修复前已存在的用户工作，不计入本表 38 项修复。

## 已知限制与兼容性

1. 当前托管沙箱中 LeakSanitizer 无法枚举进程，因此回归脚本默认 `ASAN_OPTIONS=detect_leaks=0`；在非受限环境可运行 `BUGFIX_DETECT_LEAKS=1 bash PeTar_new_version/test/bugfix/run.sh`。
2. 尚未完成真实多 rank MPI 数值对照、Helgrind 全局竞争扫描、`NJMAX` 边界力值压力测试及连续/重启随机序列对照。
3. B34 的 `ARI2` 和 B39 的 `SYS2` 是新的调试二进制格式；旧文件包含不完整浮点或进程指针，现选择明确拒绝，而不是危险迁移。
4. `1013_new_4` 中现存的 hard-debug dump 是旧 `SymplecticStep` 裸指针格式，当前程序会按 B39 的安全校验明确拒绝，因此无法用新版二进制直接复放该历史事件。B41 已用同型 `OverflowError` 和近/非近共振根的注入回归验证；仍建议用兼容当前格式的新生产运行复核总体 remnant 统计分布。
5. `inspiralZJ.cpp` 与 bundled `eccentricity.py` 原文件使用 CRLF；本次保留其换行风格，未做无关的全文件格式化。B41 未修改 Python 源码。

## 源码版本说明

本地 `main`、`HEAD` 与本地保存的 `origin/main` 均为：

`06fe527376ecf6aab627155b570e0c2b7de86a68`（2026-10-21，`a new function, refind tree time step when cluster calculation too expansive`）。

2026-11-04 在线执行 `git ls-remote origin refs/heads/main` 成功，远端返回同一提交 `06fe527376ecf6aab627155b570e0c2b7de86a68`；因此已确认当前源码基线与 GitHub 仓库最新 `main` 一致。安装的程序还包含本表所列、尚未提交的工作区修复。

工作区在修复前已有用户修改（包括 Peters 合并时间相关文件、`README.md`、`Makefile.in` 和 `ar_interaction.hpp` 的部分内容）；本次修复保留了这些修改，没有重置或覆盖。
