# Bug 最小修复汇总

首次修复日期：2026-11-04；最近补充修复：2026-11-10（后续事件记账和警告修改见文末）。

依据：`bugs/BUG_TABLE.md` 中 B01–B40，以及后续报告的近共振异常 B41、SIMD 静态 TLS 启动崩溃 B42、hard 能量并发汇总丢账 B43，以及 BSE GW 与动力学并合时钟/日志混用 B44，以及逃逸双星拆分与移除能量错误 B45。

状态统计：**已修复 41 项；当前版本未复现/路径已串行 2 项；部分修复 2 项。**

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
| B09 | petar | 按用户要求回退 | 2026-11-08：恢复 SIMD 邻居搜索及组合力内核的邻居计数，移除两处标量双精度复核；保留堆分配与对齐等其他修复。 | 优先效率；不再保证大绝对坐标下的单精度近邻判定平移稳健性。 |
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
| B43 | petar / hard 能量汇总 | 已修复 | 为 `SystemHard::driveForMultiClusterOMP()` 和 `finishIntegrateInterruptClustersOMP()` 中共享的 `energy += ...` 增加同名 `omp critical(petar_hard_energy_accumulate)`，保护完整 `HardEnergy` 累加，防止十个能量字段被并发覆盖。 | 完整构建通过；实际四线程积分中的累计不一致次数由修复前 755 次降为 0，末尾残差与单线程参考一致。已备份并更新安装主程序；历史模拟尚未重跑，中断恢复分支尚未专项运行验证。 |
| B44 | petar / 并合预测与 change | 已修复 | 分离绝对并合截止时刻与动力学轨道检查点；纯 BSE GW 不重置截止时刻；BSE 前检测动力学，BSE 后仅同步检查点；跨 BSE 保留小扰动累积；change 排除内部耗散更新。 | 完整构建、扩展 orbit_change 及 tide/frozen/Peters 回归通过；6 组实际 BSE 回调零误报，6 组按固定截止时刻并合；未重跑完整星团；后续已按用户要求备份并替换安装程序。 |
| B45 | petar / 逃逸删除 | 已修复 | 当前态直接势能；互认且束缚的双星按质心成组删除；全局删除集合的内部势能只计一次；跨MPI一致修正参考能量。 | 1/2/4 MPI进程各13场景通过；3个N500异常区间短重跑误差大幅降低，目标双星同刻删除。详见 DEBUG_ESCAPE_REMOVAL.log。 |

## B45 修复与 debug（2026-11-07）

N500 AiywPI 在约40.4 trh输出出现 +1.982011 |E0| 跳变。逃逸记录确认78/212仍束缚，删除时缓存势能滞后；结合整对删除时的势能重复计数，预测跳变+1.982012，与实测吻合。另两个区间涉及双星拆分/分时删除。

新增 `src/escape_removal.hpp`，重写 `PeTar::removeParticles()` 的选择和记账阶段：收集跨rank候选与伴星，验证互认ID和当前束缚能；以质心半径/外部能量成对决定，避免内部速度和负束缚能混入星团逃逸判据。只计算候选粒子的当前全局势能；全局删除能量使用 `-sum(K_i+m_i phi_i)+U_SS`。显式正质量删除同样纳入集合，零质量残骸不重复记账。

实际生产路径回归在MPI=1/2/4、OMP=4下各13项通过，初末直接能量差与记账一致至2e-12阈值内。测试覆盖跨rank配对、空rank、两成员跨半径边界、已解体双星、强制删除、软化、两对双星同时移除等。

三个ASCII快照短重跑的最大相对残差（以重跑初始能量归一化）：

| 起始快照 | 修复前 | 修复后 |
|---|---:|---:|
| data.159 | 0.449672 | 0.000184791 |
| data.201 | 2.931106 | 0.000007990 |
| data.495 | 0.215870 | 0.000018478 |

修复后500/233、78/212、499/64均在各自同一时刻删除。源码构建通过。
首次调试中解决了FDPS数组归约接口不兼容，以及测试夹具未初始化FDPS树导致的析构错误；没有修改FDPS本身。MPI因沙箱socket限制而申请权限后运行。

完整过程见根目录 `DEBUG_ESCAPE_REMOVAL.log`；自动化验证与详细输出在 `diagnostics/escape_removal_fix/`；可复现测试与说明在 `PeTar_new_version/test/escape_removal.cxx/.md`。
2026-11-07 16:30（Asia/Shanghai）已备份并替换安装主程序；SHA256为 `68fbde12bd3c44996a0ebc68c3ad44c35d47652272dd931ed96177c9ae5b53a5`，备份路径为 `/home/joojie/bin/petar.mpi.omp.avx2.bse.tt2nd.before-escape-removal-fix-2956d79602ac`。安装记录见 `diagnostics/escape_removal_fix/installation.json`。

只处理互认伴星且当前束缚的双星，不将整个层级多体系统自动作为一个删除对象。可选Galpy/DM分支未专项验证，完整N500长程模拟未重跑，原运行保持不变。

## B44 修复过程（2026-11-07）

原代码使用 `abs(Δa)>1e-3*a || abs(Δe)>1e-3` 决定 BSE 后是否重建预测，
既可能遗漏高偏心率下显著的小 Δe，又可能在反复重建 `now + lifetime` 时移动截止时刻。
此外，轨道基准混用使 BSE 的小改变量可能进入 change，大更新也可能覆盖此前未记录的动力学变化。

修改 `ar_interaction.hpp`、`particle_base.hpp`、`bhmerger_output.hpp`：

- 将绝对时钟与每次动力学采样的轨道检查点分开。纯 BSE 更新保留时钟，不使用旧的绝对 Δe 阈值。
- BSE 前采样动力学变化，内部更新后仅同步检查点；拒绝的 BSE Roche/coalescence 提议不重启时钟。
- 以完整 Peters 因子的 `log(a^4 F(e))` 变化计算有符号累积量，净寿命比小于 0.9 或大于 1.1 才写 change 和重算时钟；真实动力学仍可提前或延后并合。
- 小动力学增量跨 BSE 保留；成员交换、质量/类型变化、解体正确重建或清理。双曲 GW 捕获建立预测但不写动力学 change。
- change 保留前 18 列，追加 `dynamics_dlog_tgw`。旧/新轨道现在是相邻动力学检查点，不再保证旧时钟等于旧检查点时间加该轨道的孤立寿命，不能用此关系标定新日志的等时线。

验证：`make -C PeTar_new_version -j2 all build/petar.orbit_change.test` 成功。
扩展回归包含 2000 次内部更新无重置/误报、高 e 小变化、动力学提前与延后、跨内部更新的累积、交换和内部解体/捕获。
实际 Fortran BSE + AR 回调覆盖 e=0、0.9、0.9999，各 10/100 次初始回调：6 组纯 GW 均零 change，施加速度扰动后各一条；另外 6 组持续回调直到并合，均在保留的截止时刻触发。
`tide_collision`、`frozen_binary`、`peters_merger_time` 回归均通过。
实际并合测试中 Python remnant 模块出现既有 near-resonance 恢复信息与 sqrt RuntimeWarning，测试退出码均为 0。
详细日志及构建哈希位于 `diagnostics/dynamics_only_change/`（该目录被 Git 忽略），可复现测试在 `PeTar_new_version/test/orbit_change.cxx`。

限制：未做完整 N=4000 重跑或多进程星团验证；动力学诊断仍包含积分误差，不能直接证明某个第三体的贡献。
新增状态改变原始二进制粒子布局，旧二进制快照/hard dump 不兼容；ASCII 重启保持原有重新建立时钟行为。
新旧 change 不能混合追加。源码和工作区构建已更新；随后于 2026-11-07 15:44（Asia/Shanghai）按用户要求替换 `/home/joojie/bin/petar.mpi.omp.avx2.bse.tt2nd`，`petar` 软链接保持指向此文件。安装文件与已测试构建的 SHA256 均为 `2956d79602ac75058360d233a0101e1ae6ed78718e6fafafba0cf872876bb1fe`。旧程序备份为 `/home/joojie/bin/petar.mpi.omp.avx2.bse.tt2nd.before-dynamics-clock-fix-a78dd82f1ca4`；安装记录见 `diagnostics/dynamics_only_change/installation.json`。

## B43 修复过程（2026-11-07）

### 发现与定位

在比较 N=500、1000、2000、4000、10000 的能量曲线时，发现 N=1000 在 ID 39/944 并合所在输出区间出现约 −44.33% 的相对残差跳变。该残差已经扣除了代码记录的 interrupt 等物理修正，不能直接解释为引力波辐射损失。

先复放原始 `dump_merger.0.3.1.1793873237`：该步骤从 `0.0003262832760810852 Myr` 开始，包含五次轨道修改和一次并合。局部 interrupt 净修正约为 `+1.9823437109e9`；加入质心动能修正后约为 `+1.9822877063e9`，与直接计算的初末粒子动能加势能之差一致，相对差约 `3.6e-16`。因此，这次局部并合步骤没有显示出修正公式本身的错误。

进一步检查发现，[hard.hpp](PeTar_new_version/src/hard.hpp) 的两个 OpenMP 并行循环均直接写共享成员 `SystemHard::energy`：

- `driveForMultiClusterOMP()`：各 cluster 完成积分后汇总能量。
- `finishIntegrateInterruptClustersOMP()`：中断 cluster 恢复并完成积分后汇总能量。

`HardEnergy::operator+=` 对十个浮点字段逐一执行普通读、加、写，没有同步保护。不同线程可覆盖彼此的贡献，且 total modify、interrupt、hard error 和 slowdown 字段可能丢失不同的量。循环末尾的 barrier 不能避免循环内部的数据竞争。这与 B01 所检查的单个 Hermite 系统内部按组串行积分不是同一条路径。

### 最小修复

在上述两处 `energy += ...` 外使用同名临界区：

```cpp
#pragma omp critical(petar_hard_energy_accumulate)
{
    energy += /* 当前 cluster 的 HardEnergy */;
}
```

锁保护的是整笔十字段汇总，未改变并合物理公式、积分参数或已有模拟数据。也未将其他 PROFILE 计数的并发问题混入此次修复。

### 验证证据

1. **并发压力复现：** 使用实际 `HardEnergy` 类型，四线程累计 200000 个可精确表示的正负整数贡献，正确和为 100000。无保护版本五次试验的 total modify 为 73990、83993、80203、79830、80256，各能量字段还出现不同差额；加锁后五次试验的各检查字段均为 100000。该压力测试比例不代表原星团运行的丢账比例。
2. **真实积分路径：** 在诊断副本中另设线程独立的旁路累计器，dispatch 结束后串行求和，与原共享累计比较。从原 N=1000 的 `data.27` 短程运行，四线程无保护版本出现 755 次不一致；单线程和加锁后的四线程均为 0 次。此段当前代码重跑未重现原并合历史，用于核对汇总路径，不用于解释全部历史残差。
3. **正式修复验证：** 应用两处修补后，`make -C PeTar_new_version -j2 all` 成功。再从正式源码制作旁路核对版本，四线程短程运行完成，累计不一致为 0；末尾相对残差为 `−8.53631093197093e−6`，与单线程参考一致。
4. **覆盖范围：** 普通多 cluster 汇总路径已运行验证；中断恢复路径采用相同保护，已编译和静态检查，但本次默认 interrupt mode=1 的短程运行未专项触发该恢复分支。未将 2026-11-04 的整套 sanitizer 回归声称为本次重新执行的结果。

### 安装与结果解释

已备份并更新 `/home/joojie/bin/petar` 指向的主程序 `/home/joojie/bin/petar.mpi.omp.avx2.bse.tt2nd`。旧程序备份为：

```text
/home/joojie/bin/petar.mpi.omp.avx2.bse.tt2nd.before-energy-fix-aefbd62761
```

安装文件与工作区构建文件的 SHA256 一致：

```text
a78dd82f1ca4f60161046876b341d3e696692ba29ca7f590662306ec2147a4a6
```

这项修复确认解决共享汇总的数据竞争，但尚未定量证明原 −44.33% 跳变全部由此造成。`error_hard_cum` 同样经过有缺陷的汇总路径，因此旧日志中总误差与 hard 误差的相关性不能单独证明积分器精度不足。旧输出缺少每个线程每步的贡献，不能可靠地事后修补；需用修复后的程序重跑，再更新能量曲线和物理解读。

本地证据见 [核查报告](diagnostics/merger_energy_audit/REPORT.md)、[修复报告](diagnostics/merger_energy_fix/REPORT.md)、[验证数据](diagnostics/merger_energy_fix/validation.json) 和 [安装记录](diagnostics/merger_energy_fix/installation.json)。这些诊断目录被 `.gitignore` 排除，不随源码提交；关键发现与验证结果已在本节保留。

## 初次验证结果（2026-11-04，B01–B42）

- `make -C PeTar_new_version -j2 all`：通过；生成默认 MPI/OpenMP/AVX2/BSE 主程序、hard-debug 和格式转换工具。
- `bash PeTar_new_version/test/bugfix/run.sh`：通过：
  - 核心结构/序列化 ASan+UBSan 回归；
  - PeTar/FDPS/BSE/OpenMP 集成 ASan+UBSan 回归；
  - SIMD 工作区 64 字节对齐、堆生命周期、标量数值对照及静态 TLS 大小回归；
  - Python remnant 异常、近共振恢复/非共振拒绝恢复与引用保留量回归；
  - Python 进动相位符号回归。
- 可选宏语法编译均通过：`AR_TIME_FUNCTION_MULTI_R`、`AR_SLOWDOWN_MASSRATIO`（含二者组合）、`AR_CHANGEOVER`、`INTEGRATED_CUTOFF_FUNCTION`、`HERMITE_DEBUG`、`USE_QUAD+AVX_PRELOAD`。
- 现有 Peters 合并时间测试通过：最大采样相对误差 `0.00113518`；渐近切换跳变 `0.020585`。该功能及其测试属于修复前已存在的用户工作，不计入当时的 38 项修复或本次更新后的 39 项修复。

## 已知限制与兼容性

1. 当前托管沙箱中 LeakSanitizer 无法枚举进程，因此回归脚本默认 `ASAN_OPTIONS=detect_leaks=0`；在非受限环境可运行 `BUGFIX_DETECT_LEAKS=1 bash PeTar_new_version/test/bugfix/run.sh`。
2. 尚未完成真实多 rank MPI 数值对照、Helgrind 全局竞争扫描、`NJMAX` 边界力值压力测试及连续/重启随机序列对照。
3. B34 的 `ARI2` 和 B39 的 `SYS2` 是新的调试二进制格式；旧文件包含不完整浮点或进程指针，现选择明确拒绝，而不是危险迁移。
4. `1013_new_4` 中现存的 hard-debug dump 是旧 `SymplecticStep` 裸指针格式，当前程序会按 B39 的安全校验明确拒绝，因此无法用新版二进制直接复放该历史事件。B41 已用同型 `OverflowError` 和近/非近共振根的注入回归验证；仍建议用兼容当前格式的新生产运行复核总体 remnant 统计分布。
5. `inspiralZJ.cpp` 与 bundled `eccentricity.py` 原文件使用 CRLF；本次保留其换行风格，未做无关的全文件格式化。B41 未修改 Python 源码。

## 初次修复的源码版本说明（2026-11-04）

当时本地 `main`、`HEAD` 与本地保存的 `origin/main` 均为：

`06fe527376ecf6aab627155b570e0c2b7de86a68`（2026-10-21，`a new function, refind tree time step when cluster calculation too expansive`）。

2026-11-04 在线执行 `git ls-remote origin refs/heads/main` 成功，远端返回同一提交 `06fe527376ecf6aab627155b570e0c2b7de86a68`；因此确认了当时的源码基线与 GitHub 仓库 `main` 一致。当时安装的程序还包含本表 B01–B42 中已实施、尚未提交的工作区修复；此处不是 2026-11-07 的远端状态核对。

工作区在修复前已有用户修改（包括 Peters 合并时间相关文件、`README.md`、`Makefile.in` 和 `ar_interaction.hpp` 的部分内容）；本次修复保留了这些修改，没有重置或覆盖。

## 2026-11-07：slowdown / 三体负时间步稳定性修复

已完成原 N=2000 large-step dump 的回归。三体及以上系统自动启用持续的局部误差控制，
保留 slowdown；另加 TTL 分母恢复和 64 次非正时间步的诊断退出保护，避免无限缩步。
运行日志新增 `AR_ADAPTIVE_MULTIPLE`、`AR_SLOWDOWN_RECOVERY`、
`AR_SLOWDOWN_RECOVERY_FAILED`、`AR_TIME_STALL`。
默认/严格容差重放的相对能量残差为 6.89e-8 / 6.17e-8，均正常完成。
详细修改、发生率限制、故障注入和部署状态见 [slowdown debug 记录](bugs/SLOWDOWN_DEBUG.md)。

## ASTRA 登记册修复进展（2026-11-07）

此节为本轮增量，保留先前修复日志及用户既有更改。用户授权所有 B 类修复，直接影响并合的 A 类逐项确认；已收到 IO-01、BSE-06/07、BSI-01/03/13/14 的批准。以下记录包含已实施但尚待专项验证的条目，不能将“已改源码”一律解读为完成生产验证。全表当前状态及未确认项保留在 [astra_bug_analyz.md](astra_bug_analyz.md)。

| 原编号 | 本轮修改与边界 | 验证级别 |
|---|---|---|
| ASTRA-FD-01 | 拥有的裸指针初始化为 nullptr，未 initialize 的析构路径不再读取未初始化指针。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-16 | 析构补齐 adr_* 工作区及 MPI status 数组释放。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-18 | 粒子交换进程计数 n_proc_exch_ 初始化为 0。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-19 | TimeProfile 构造函数初始化五个 calc_force__core 子计时器；新增构造回归。 | 针对性回归 |
| ASTRA-FD-17 | 线程工作区内存计数改为累加；尚未完成所有遗漏数组的完整核算，本项仍有残留。 | 部分修复，内存核算仍待补齐 |
| ASTRA-FD-11 | MemoryPool::clearAll 保留已分配底址，将 top 重置为 bottom，避免仍有容量但基址为空。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-10 | reserve 先用局部变量检查新容量，再交给 ReallocInner 扩容，保留旧容量供拷贝；新增保值回归。 | 针对性回归 |
| ASTRA-FD-09 | 重采样前恢复已有本地样本的逻辑长度，再追加第二种粒子系统的样本。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-21 | Matrix2 的 xy 乘积修正；新增四分量矩阵乘法数值回归。 | 针对性回归 |
| ASTRA-FD-05 | 周期镜像直接力缓冲索引修正为 j*n_recv+k。 | 源码核对/默认构建；待专项验证 |
| ASTRA-IO-04 | FileHeaderNoOffset 按自身 sizeof 读写。 | 源码核对/默认构建；待专项验证 |
| ASTRA-H4-07 | HermiteInformation 的 time_origin 构造和 clear 均初始化为 0。 | 源码核对/默认构建；待专项验证 |
| ASTRA-H4-05 | used_mask 改用覆盖 index_offset_group_+n_group 的动态存储，消除旧栈数组长度不足。 | 源码核对/默认构建；待专项验证 |
| ASTRA-AR-13 | 潮汐张量最近邻索引为 -1 时不再解引用。 | 源码核对/默认构建；待专项验证 |
| ASTRA-DG-16 | 粒子和 hard dump 的读取数错误提示分别修正为 5 和 4。 | 源码核对/默认构建；待专项验证 |
| ASTRA-DG-09 | Tprofile::end 只取一次时间用于 time/tbar。 | 源码核对/默认构建；待专项验证 |
| ASTRA-DG-06 | 邻居诊断将平方距离与平方外半径比较。 | 源码核对/默认构建；待专项验证 |
| ASTRA-EN-06 | 删除遮蔽成员的局部 mass_modify_list；多 rank 中断场景尚未专项运行。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-09 | StarParameter::initial 显式初始化 lum。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-07 | bse_test 的第二颗星备份使用当前 bin[i]。 | 源码核对/默认构建；待专项验证 |
| ASTRA-TR-10 | Kepler 求解使用有符号 cos 导数，并加入区间保护与浮点收敛判据；e=0..0.99 的 505 个样例核对方程残差。 | 针对性回归 |
| ASTRA-TR-11 | 邻居速度/近心判据传入运行时 gravitational_constant；G 非 1 时会改变分组，尚未做星团轨迹数值对照。 | 源码核对/默认构建；待专项验证 |
| ASTRA-AR-12 | R4 扰动指标分母改为 G²，单位缩放结果会改变。 | 源码核对/默认构建；待专项验证 |
| ASTRA-TR-17 | 相对加速度改为 ai-aj，消除 ai-ai 恒零。 | 源码核对/默认构建；待专项验证 |
| ASTRA-TR-19 | quad 分块预装使用块内 i_tmp 索引。 | 源码核对/默认构建；待专项验证 |
| ASTRA-H4-06 | 力与粒子字段初始化，clear/clearForce 清除 KDKDK4 的 acorr。 | 源码核对/默认构建；待专项验证 |
| ASTRA-TR-15 | FPSoft 的邻居计数、rank、地址等瞬态字段显式初始化。 | 源码核对/默认构建；待专项验证 |
| ASTRA-IO-06 | FPSoft 二进制尾部逐字段写入 acc、势和 S64 邻居数，避免 acorr 插入改变布局；旧的错误 KDKDK4 快照不能恢复丢失字段。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-03 | 参数文件读取只受 rank==0 控制，print_flag 仅控制打印。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-04 | 数据格式使用运行时范围检查，拒绝 0..3 以外取值。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-09 | 与 BL-04 共用一个范围校验修复，不重复计为独立补丁。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-12 | regularTimeStep 拒绝非有限或非正输入；单粒子自动步长经过同一校验。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-08 | bse_test 双星演化串行调用共享 Fortran 状态，避免 OpenMP 公共块及容器竞争。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSE-08 | 双 kw=15 的提前返回设置 event 结束索引为 0。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-06 | BSEEMP 金属丰度说明与校验统一为有限正数。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-10 | mlwind 的 VALUE1 公共块移除不匹配的 mxns 字段。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-01 | 经用户批准，接触分支潮汐因子改为 (ttid/(2tc))²；编译通过，尚未完成长期双星潮汐轨迹对照。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-07 | hard dump 初始化先 clear，线程数及索引验证后再访问；全局对象采用 C++17 inline。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-10 | dump 数量在 resize/read 前检查为非负。 | 源码核对/默认构建；待专项验证 |
| ASTRA-DG-07 | dump_number 使用原子递增。 | 源码核对/默认构建；待专项验证 |
| ASTRA-DG-12 | 与 BL-10 共用负数量前置检查。 | 源码核对/默认构建；待专项验证 |
| ASTRA-DG-14 | ASSERT 宏改为 do/while(0)。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BL-11 | 静态变量定义改为 C++17 inline，消除多翻译单元重复定义。 | 源码核对/默认构建；待专项验证 |
| ASTRA-AR-10 | 首个过冲步使用该步开始时间计算缩步比例。 | 源码核对/默认构建；待专项验证 |
| ASTRA-TR-14 | 写回 n_ngb 包含自身，与力计算计数约定一致。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-06 | LOOP_TREE 的七组固定 12 线程工作区改为按 n_thread 动态分配。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-08 | Fisher–Yates 抽样每步使用区间 [i,n_src-1]，避免固定区间引入偏差。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-13 | 访问 adr[n_cnt-1] 前显式检查 n_cnt>0，无候选时诊断退出。 | 源码核对/默认构建；待专项验证 |
| ASTRA-FD-02 | 零质量 MomentMonopole 不再执行 pos/mass。 | 源码核对/默认构建；待专项验证 |
| ASTRA-IO-01 | 经用户批准，ASCII #BHMERGER_V2 保存 deadline、a/e、检查点时间及累计动力学量；兼容 V1/无后缀，旧文件缺失历史仍无法补回。新增重启后首次动力学检查保持 deadline 的回归。 | 针对性回归 |
| ASTRA-BSI-03 | 经用户批准，bse/bseEmp 的三个 CE kick 调用设置 SNEXPLOSION，仅无超新星时用旧能量更新 ECC；MOBSE 原有守卫保留。编译及控制流核对通过，尚未构造真实 CE 内超新星存活样例。 | 源码核对/默认构建；待专项验证 |
| ASTRA-BSI-14 | 经用户批准，bse/bseEmp/mobse 在抽样长度为 0 时输出零 kick，保留 RNG 调用次数；现有处方均对抽样长度作缩放。sigma=0/完全回落/普通非零 kick 回归通过；旧代码在同一零 kick 测试触发 SIGFPE。 | 针对性回归 |
| ASTRA-BSI-13 | 经用户批准，bse/bseEmp 的 rapid/delayed 分支共用 bounded_fallback：先计算回落质量，再限制到 [0,包层质量]；零包层返回零质量/比例，代数消去薄包层除法。实际 Fortran 边界测试覆盖 84 组输入，普通公式两例一致。 | 针对性回归 |
| ASTRA-BSE-06 | 经用户批准，同 BSE-07，在 CE/MIX 后刷新 menv/renv 并清空消失星的结构量；真实 MIX/CE 回归核对新包层参数。 | 针对性回归 |
| ASTRA-BSE-07 | 经用户批准，新增 merge_structure 返回刷新后的结构、光度、寿命；保留旧 merge ABI，自旋按 CE 返回角动量或 MIX 保留星原角动量更新。BH 拟合质量覆盖后同步 BH 核质量/半径。修正 MOBSE merge 中 hrdiag 参数签名。真实 MIX、CE 和 BH 合并结构回归通过；MIX 自旋角动量保持。 | 针对性回归 |

### 数值影响与回归记录

- fallback 测试：7 个核质量 × 6 个包层质量 × 2 个处方，共 84 组；包层包含 0、1e-300、0.01、0.1、1、10。检查有限性、比例/质量范围和质量一致性。薄包层的回落质量会从原公式的超额值降至可用包层质量，故 remnant mass 与 kick 均可能改变；正常未截断两例保持原公式。
- kick 测试：sigma=0 和完全回落得到四分量全零；普通固定种子样例为 `(6.28253038533489, 92.1658691558563, 44.1394354019848, 102.383140130117)` km/s。用修复前 kick.f 跑同一零弥散测试，在启用浮点异常捕获时退出 136（SIGFPE）；修复后通过。
- 恒星并合测试：真实 MIX 后 mass/radius/lum 为 `1.8 / 1.07144 / 17.249`，真实 CE 后为 `4.74121 / 25.0182 / 4254.69`（BSE 太阳单位）；与零时长重新演化的结构参数一致。另检查 MIX 保留星自旋角动量及 BH 合并后的质量/半径一致性。这里是测试夹具结果，不是生产并合率预测。
- Kepler 测试：e=0、0.3、0.5、0.9、0.99，各 101 个平均近点角；检查 `E-e*sin(E)-M` 残差 <1e-10。只修 cos 符号的中间版本仍出现不收敛，最终加入区间保护。
- IO-01：ASCII 往返后首次调用动力学检查，deadline 保持不变，检查点时间前进，累计动力学量保留。V1/无后缀快照仍可读，但不能恢复缺失历史；本项不补全 IO-03 的 RNG 状态契约。
- 通用集成测试的逃逸夹具从重合点改为 256 个不重合粒子，以解析两两势能校核删除记账，并对照 1/4 线程。旧夹具因 B45 正确拒绝无软化重合点而中止；修改夹具没有放宽能量断言。
- BSEEMP 修改文件和 MOBSE merge/kick 单独 Fortran 编译通过；尚未执行这两种后端的完整演化回归。MOBSE 原 hrdiag 调用参数与实际签名不符，本轮在刷新结构时一并对齐。
- CE 内超新星后仍存活的真实事件、长时潮汐双星、非默认力核/LOOP_TREE、完整多 rank 和星团长期历史仍待专项验证。TR-11、AR-12、TR-17 等变化已在表中注明，尚无可靠生产轨迹差值，不虚构数值变化幅度。

可重复测试入口：`bash PeTar_new_version/test/bugfix/run.sh`。新增 `sn_kick.f90` 使用实际生产 Fortran 例程及 invalid/zero/overflow 异常捕获；`stellar_merge.cxx` 调用实际 BSEManager。没有安装到用户 bin 目录或提交代码。

### 后续逐项授权的增量

| 原编号 | 状态 | 修改与验证 |
|---|---|---|
| ASTRA-BH-05 | 已修复（专项回归） | 经逐项批准，仅圆轨道允许 e 的 (1,1) 零常量输出；其他轨道字段读取末端样本。注入合法圆轨道结果后，remnant 拟合收到末端 deltaphi=0.375，而不是遭遇时角度。 |
| ASTRA-BH-08 | 已修复（专项回归） | 经逐项批准，e² 仅在 32 epsilon 的消去舍入容差内裁为 0；明显负值或非有限项返回失败，不误走直接并合。边界和明显错误输入回归通过。 |
| ASTRA-BH-09 | 已修复（专项回归） | 经逐项批准，构建嵌入绝对源码包/安装包路径，支持 PETAR_PRECESSION_PATH 覆盖；新增安装包目标。从 /tmp 启动实际导入源码包和暂存包通过，未安装到用户 bin。 |
| ASTRA-BH-10 | 部分修复 | 经批准只处理角度端点舍入：Python arccos 共享有限性/32 epsilon 边界检查，明显越界报错；原来的无条件 clip 同步收紧。端点和相位符号回归通过。零自旋及未证明为重根的退化问题仍未关闭。 |
| ASTRA-BH-11 | 已修复（专项回归） | 经逐项批准，解释器存在与 NumPy API/路径初始化分开记录，在 GIL 下初始化桥接。宿主预初始化与独立解释器两种测试通过；finalizePython 不关闭宿主解释器。宿主自行反复销毁/重建 NumPy 解释器不在支持范围内。 |
| ASTRA-AR-08 | 已修复（专项回归） | 经逐项批准，ARRAY 单存活粒子漂移使用 dt/κ，物理时间仍增加 dt。实际 integrator 中注入并合事件，κ=1/8 时的位移分别为 (2,4,6)/(0.25,0.5,0.75)，结束时间均为 2。 |
| ASTRA-AR-09 | 已修复（专项回归） | 经逐项批准，零质量人工节点清空轨道缓存，重算质量/质心/当前分离并标为不稳定；保留正质量子树。周期更新跳过人工节点，Hermite 拆组不再消费其虚假轨道。4 正质量→2 正+2 零、单存活、全零及重合节点回归通过；完整生产并合/拆组仍待长程验证。 |
| ASTRA-DG-10 | 已修复（专项回归） | SysCounts/SysProfile 用显式成员指针表遍历，消除跨独立成员的指针算术；保持列顺序及 DARKMATTER 条件列。计数 clear 回归通过。 |
| ASTRA-DG-15 | 已修复（专项回归） | Min/Max 与普通行均统计 time，barrier 列独立归约 tbar；单 rank 对照通过。多 rank 的 time/tbar 极值可以来自不同 rank，不代表同一 rank 的联合样本。 |
| ASTRA-DG-05 | 已修改，待专项验证 | Large cluster 提示受 HARD_DEBUG_PRINT 控制，不再无条件输出。 |
| ASTRA-FD-12 | 已修改，待专项验证 | CalcRank1D 沿单坐标方向搜索时检查 rank_axis 边界，越出该方向首末域则明确退出，避免先越界访问 domain。 |
| ASTRA-FD-20 | 已修改，待专项验证 | 周期边界的无 cutoff 拒绝检查同时覆盖 LONG、LONG_SYMMETRY、LONG_SCATTER，避免进入不支持 shift 的路径。 |
| ASTRA-MC-03 | 已修复（专项回归） | 经逐项批准，送入 BSE 的周期使用实际太阳质量，返回半长轴使用演化后质量。链接包装实际 get_binaries 验证 1 AU、1+1 Msun 的周期为 365.25/sqrt(2) 天；质量减半且周期不变时 a 缩小为 cbrt(0.5) AU。 |
| ASTRA-MC-14 | 已修复（专项回归） | 经逐项批准，main 在 SSE 构建中首次演化前调用 instar。碰撞表关键项验证为 ktype(0,0)=1、ktype(10,10)=15，SSE 初态生成冒烟测试通过。 |
| ASTRA-MC-16 | 已修复（专项回归） | main.c 在系统头前定义 _GNU_SOURCE，main.h 补 output5 原型；当前编译器的 gnu2x + implicit-function-declaration 设为错误的检查通过。未声称已使用 GCC15。 |
| ASTRA-MC-05 | 已修复（专项回归） | 经逐项批准，G 的单位明确为 pc³/(Msun Myr²)，逃逸速度转为 km/s，估算式补 G；直接比较 kick，任一分量超限即重抽。每次重抽恢复原始双星质心和恒星记录。仅第一组分超限的注入测试验证恰好重抽一次、质心不漂移。 |
| ASTRA-MC-06 | 部分修复 | 去掉 generate_m2 对 evolv1 的 if(epoch) 限制，零年龄幂律质量函数的 R/L 均为有限正数。f=2 单质量快捷分支仍未走 SSE 结构初始化，因此本项不整体关闭。 |
| ASTRA-MC-07 | 已修复（专项回归） | f=2 单质量分支以实际单星质量填 star[][0] 和初始质量字段，随后统一归一化；N=20、单质量 2 Msun 输出全部质量为 2，未再触发半质量半径除零。 |
| ASTRA-MC-10 | 已修复（专项回归） | 解析器接受帮助及实现支持的 0..4；使用 -M 20 -f 3 和 -N 20 -f 4 的实际生成测试通过。 |
| ASTRA-MC-19 | 已修改，待专项验证 | fractalize radial 分支逐个粒子归一化速度，避免循环结束后访问 star[N]；平均速度为 0 时不除零。 |
| ASTRA-MC-20 | 已修复（专项回归） | 输出文件名按输入长度配置缓冲，使用 snprintf，不再修改 argv；创建文件失败时明确报错。超过 90 字符前缀的 info/output3/4/5 测试通过。 |
| ASTRA-MC-21 | 已修改，待专项验证 | 最终能量检查补齐 i<j 的全部配对，新增 --check-energy 入口；与 MC-08 气体重标逻辑分开，后者尚未修复。 |
| ASTRA-MC-28 | 已修改，待专项验证 | 源文件的 MC-16 声明修复同时覆盖测试编译入口；现代 C 严格声明检查通过，自带 DM 测试结果另记。 |

MC-04 按用户明确选择暂不修改；所有未确认项继续留在分析表。以上条目中的“部分修复”与“待专项验证”均不是完全关闭。

BH-09 的安装目标为 `install_bhmerger_python`，将 bundled 包放入 `<prefix>/share/petar/python/precession`；`make install` 依赖该目标。测试通过 `DESTDIR=/tmp/astra-package-stage` 暂存。路径覆盖变量指向包含 `precession` 目录的父目录，例如 `PETAR_PRECESSION_PATH=/path/to/python-packages`。没有修改用户 bin 目录。

MC-03 的修复不仅补上输入周期的总质量因子，还在返回时使用演化后的质量。MC-05 的重抽恢复避免旧位置更新叠加；此改变有独立质心回归。MC-04 的并合/解绑产物重建仍未实施，不保证上述修复能让该未支持路径运行成功。

验证入口：`bash PeTar_new_version/test/bugfix/run.sh`、`bash mcluster/tests/run_astra.sh`。当前运行环境的 GNU C2x 严格隐式声明检查已通过，但没有借此宣称完整 GCC15 平台验证。Python 源文件与桥接 C++ 保留原 CRLF；检查空白使用 `git -c core.whitespace=cr-at-eol diff --check`。

### FDPS 与输出增量（2026-11-07）

- **ASTRA-FD-03（FDPS）：部分修复。** REUSE_LIST 在粒子数与旧排列长度不符时明确中止，防止读取新增未初始化索引；相同数量但成员/顺序变化仍由调用方遵守复用契约。
- **ASTRA-FD-04（FDPS）：已修改，专项回归通过。** 零 opening-angle 列表按 tp 地址解码 EP/SP，并加入 SP；混合 SP+两 EP 的集成测试通过，多 rank 生产运行未重跑。
- **ASTRA-FD-14（FDPS）：部分修复。** 旧交互表变体修正 EP/SP 地址解码及跳过 SP 时的循环推进；MakeIPGroupLongGLBTreeCellAsIPGBox 的地址空间问题仍保留。
- **ASTRA-FD-15（FDPS）：部分修复。** 未实现的 LongCutoffP2P 辅助接口现在明确中止，避免静默返回空发送列表；完整 P2P 功能仍未实现。
- **ASTRA-FD-22（FDPS）：已修改，专项回归通过。** Morton 坐标在整数转换前检查有限性、范围和正 half_len；仅容忍边界舍入。half_len=1e-300、1、1e300 的端点测试通过；非默认周期域尚待专项验证。
- **ASTRA-FD-23（FDPS）：已修改，专项回归通过。** 矩检查按 tp 地址解码 EP/SP；混合列表质量检查通过。
- **ASTRA-IO-08（PeTar）：已修改，专项回归通过。** 初始输出检测目标与输入是否为同一文件（包括硬链接），跨 rank 广播保护标志；保留状态输出及后续编号，新前缀仍可写初始快照。
- **ASTRA-DG-04（PeTar）：已修改，专项回归通过。** status 表头延后至首次 output、粒子数量已知时写出；空 append 文件也补标题，非空文件不重复。

数值影响：FD-04 正确计入原先遗漏的 SP 力贡献，FD-22 对合法根盒内坐标维持 Morton 排序、拒绝明显越界输入；其余地址修复影响被触发的旧辅助路径。IO-08/DG-04 只修改输出行为。MC-04 按用户意见暂不修改，未确认条目保留。

本轮输出专项与完整 `test/bugfix/run.sh` 已通过：硬链接输入快照保留原内容，下一编号照常输出，`-w 2` 表头包含粒子列。默认 PeTar 构建通过。mcluster 普通构建在固定种子 L3 测试中间歇 SIGSEGV，ASan/UBSan 构建连续 20 次未复现；仍在定位，不能宣称 mcluster 全部回归通过。

- **ASTRA-GEN-03（PeTar）**：DomainDecision 非 4 倍进程数显式中止，避免 NDEBUG 下静默无操作；未修改域划分算法，生产默认路径不调用此辅助接口。

### PeTar 增量及范围调整（2026-11-07）

用户要求停止处理 mcluster，保留既有修改和验证记录，不再继续该模块。

- **ASTRA-BH-21（PeTar）：已修改，专项回归通过。** 互斥分支使用 elif，非法参数显式抛 TypeError，eval_q 接受 0<q<=1；平行/反平行自旋和等质量接口测试通过，不代表所有演化式的等质量奇点已解决。
- **ASTRA-BH-22（PeTar）：已修改，专项回归通过。** 几何边界改为 ±(chi1+q*chi2)/(1+q)，极点及 101×101 角度网格包络检查通过。
- **ASTRA-BH-23（PeTar）：已修改，专项回归通过。** m<=0.01 使用五阶级数，m=0 返回 3/8；保留上端 regularization。100 位高精度参考覆盖 0、1e-30 至 0.99 及切换点，最大绝对误差 4.70e-14。
- **ASTRA-BH-13（PeTar）：已修改，当前环境回归通过。** 单元素数组改用 np.asarray(Ssav).item() 提取标量；当前进动回归通过，未在 NumPy 2.x 单独验证。
- **ASTRA-BH-14（PeTar）：已修改，当前环境回归通过。** 两处 np.in1d 改为 np.isin，输出筛选语义保持；当前进动回归通过，未在 NumPy 2.x 单独验证。
- **ASTRA-TR-08（PeTar）：已修改，专项回归通过。** CM 预测位移范数修正为 dx²+dy²+dz²；非默认 SMOOTH_CM_USING_VEL_PRED 完整轨迹尚未运行。
- **ASTRA-TR-18（PeTar）：已修改，专项回归通过。** continue 首尾步使用对称分解首列系数，DKD 时对应漂移；保持中间 kick/drift 状态切换。KDK/DKD 和普通/KDKDK 系数两轮循环均验证两类步长总和。
- **ASTRA-TR-20（PeTar）：已修改，专项回归通过。** 读步长或推进前显式检查 setStep 已生成表，避免空表下标/取模；未初始化调用已验证明确中止。

数值影响：BH-22 扩展到正确几何可达区间，BH-23 修正小 m 平均自旋系数；TR-08 改变非默认 CM 平滑阈值，TR-18 修正 DKD 首尾漂移量。BH-21/13/14 的正常输入数学表达式不变。GEN-03 NDEBUG 守卫探针通过。

### 日志、生成器与构建诊断增量

- **ASTRA-BL-08（PeTar）：已修改，待专项验证。** 软化长度先检查有限非负；自动/显式 changeover 尺度确定后与实际 base r_in 比较，不再与无量纲 ratio_r_cut 比较。
- **ASTRA-BH-19（PeTar）：已修改。** 修正五处接口文档：上下界/根返回顺序、hatJ 向量及 kick 默认开关；未修改这些接口的返回值和默认值。
- **ASTRA-BH-27（PeTar）：已修改，专项回归通过。** 写失败明确报告文件名；持锁回滚普通文件的半行，回滚失败另报错；检查 close 错误。注入短写后 EIO 验证原文件保持完整，随后正常写入成功。硬件/文件系统拒绝回滚时不能保证无半行。
- **ASTRA-BH-28（PeTar）：已修改，串行专项回归通过。** rank0 在 append 新/空文件时持锁补表头，非空不重复；初始化状态广播给全部 rank，失败集体中止。新文件/重复 append/失败路径已测试，MPI 多 rank 故障注入未运行。
- **ASTRA-GEN-01（PeTar）：已修改，专项回归通过。** Kepler 抽样助手共享调用方按 rank/seed 初始化的单一随机流，去掉各自同种子的静态流；相同 seed 重复调用可重现，序列相对旧版会改变。
- **ASTRA-GEN-04（PeTar）：已修改，专项回归通过。** Plummer 增加默认 G=1 的显式参数，尺度 r_scale 包含 G，速度包含 sqrt(G*M/r_scale)；PeTar 主入口传运行时 G。保留现有截断半径处方。
- **ASTRA-BL-13（PeTar）：部分修复／澄清。** 补充 HARD_DUMP/AR_DEBUG_DUMP 功能提示。复核发现 configure.ac 已通过 AC_CONFIG_FILES 生成 get_version.hpp，版本来自 PeTar/SDAR VERSION；不能把缺 Makefile 规则单独判为未生成。已生成头不会随 VERSION 自动刷新这一构建依赖问题仍保留。

数值影响：GEN-01 更改 Kepler 初态随机序列；GEN-04 在非单位 G/M 时修正尺度与速度，默认 G=M=1 不变。BL-08 对不合法软化/半径组合明确拒绝。日志与文档修复不改变并合判据。

本批收尾验证：`make -C PeTar_new_version -j2` 和 `bash PeTar_new_version/test/bugfix/run.sh` 均退出 0；额外纯 z 方向 CM 预测探针通过，并已加入测试入口。日志中的 “Failed to write merger event log” 是故障注入测试的预期报错，测试已验证半行回滚。`git diff --check`（兼容原 CRLF 文件设置）通过。

本轮累计涉及 102 个原编号，包含部分修复、文档澄清及仍待专项验证项，**不代表 102 项全部关闭**。未确认项、未获批的 A 类及剩余已确认项仍在分析表；mcluster 已按最新指示停止处理。MPI 多 rank 日志故障、全部非默认构建及长期星团演化未全量验证。

### 新版编译安装（2026-11-07）

已按用户要求强制重编译并安装到 `/home/joojie`。配置：MPI/OpenMP/AVX2、BSE、BHMERGER、二阶潮汐张量；`petar` 指向 `petar.mpi.omp.avx2.bse.tt2nd`。主程序、hard.debug、format.transfer 的 SHA-256 与编译产物一致，安装后的 `petar -h` 在沙箱外退出 0（MPI 本机 socket 在沙箱内受限）。

旧版备份与构建、安装、校验记录：`diagnostics/astra_install_20261107_b4f4mfpc`。本次安装不代表此前待专项验证项已关闭；未编译或安装 mcluster。

### 三阶潮汐张量与四极矩安装（2026-11-07）

按用户确认配置 `--with-tidal-tensor=3rd --enable-quad`，保留 MPI/OpenMP/AVX2、BSE、BHMERGER 等选项，强制重编译并安装。`/home/joojie/bin/petar` 现指向 `petar.mpi.omp.avx2.bse.tt3rd`；安装后帮助输出确认三阶潮汐张量和四极矩。mcluster 普通版及 SSE 版亦重新编译安装，保留普通版原有指向源码构建产物的链接，SSE 版位于 `/home/joojie/bin/mcluster_sse`。5 个可执行文件校验值均与构建产物一致。本次按要求不备份，不修改 mcluster 算法。构建与校验记录：`diagnostics/petar_tt3_quad_install_wi7wfczo`。

## 2026-11-08：AR 固定步长异常回退与 change 日志优化

按用户要求，取代此前对全部多体 `later->none` 的无条件切换：正常保留
`later`，试算误差超限或时间/TTL 非法时恢复 BSE 回调之后的步级备份、减半 ds，
本区间自适应重算；不重放 BSE、随机事件或日志。删除 AR_ADAPTIVE_MULTIPLE 输出。
生产 change 改用 `data.change.<rank>.<thread>` 独占缓冲流，快照时刷新；保留
bhmerger 共享日志。提供 merge_change_logs.py 合并供原分析脚本读取。
验证与限制见 PeTar_new_version/test/ar_fixed_retry.md；故障 dump 重放通过，但
紧容差下子步数增加，未声称全星团提速。未安装或重启已有运行。

### 高频扰动小 N 性能对照（2026-11-08）

以 0d1dd21 为修改前基线，对 16/32 个 30 Msun 黑洞的高密度子系统，各取
3 个初态、每版交替重复 5 次，共 60 次成功运行；固定 BSE/Python 随机种子，
单核计时包含积分及事件输出。六组中位耗时合计下降约 1.85%，单组变化
为 -9.08% 至 +4.71%，没有一致明显加速，未进行正式统计显著性检验。
各组并合次数及有序成员配对一致，但部分 AR 子步、BSE 事件行及 change 行数改变。
五组在两版均有 hard-energy 告警，因此仅作为性能压力测试，不作为物理精度验证。
单线程未覆盖日志锁争用，不能外推 N4000/16 核收益。
完整输入、构建、原始记录、图表及限制见 diagnostics/ar_ab_20261108/REPORT.md。

### AR 热路径等价优化（2026-11-08，项目 1/2/4/5）

- `calcAccPert` 的两个扰动者指针数组改为线程独占复用缓冲区；每次重新填充，不缓存位置预测。同步力函数无递归回调，线程之间不共享缓冲。
- 并合预测拆分为更新轨道及使用现有轨道两层，删除检查/同步后立即刷新时的重复 `calcOrbit`；BSE 改变状态后的轨道刷新保留。
- 并合寿命形状函数采用两项、精确输入匹配的线程局部缓存；保留原有数学表达式、阈值和浮点运算，先读取旧检查点再计算当前值。不增加粒子字段或改变重启格式。
- 多体内力改为每对一次距离计算，保留每个目标按邻居编号累加以及原能量归约顺序。changeover 外半径相同时保留原有方向相关的选择规则。纯两体专用核未改。
- 未修改项目 3 的日志策略、BSE 调用时机、步长及物理判据。生产构建通过；带/不带 changeover 的力与旧实现对照、并合预测及缓存专项通过（ASAN/UBSAN）。
- 6 个高密度小 N 算例，每版 3 次，共 36 次：AR/H4/同步计数与能量残差一致；最终粒子及并合/change/BSE/SSE 文件逐字节一致。初测中位耗时合计下降 4.16%，部分计时与构建任务重叠，不能作为严格提速结论。原有能量告警仍在，不代表物理精度认证。
- 详情与原始输出：`diagnostics/ar_hotpath_20261108/REPORT.md`。未安装或重启运行中的模拟。

收尾：完整 `test/bugfix/run.sh`（新增带/不带 AR_CHANGEOVER 两种集成测试）退出 0；`git diff --check` 通过。日志中的写失败信息来自预期故障注入。构建和回归日志已保存至上述诊断目录。


## 2026-11-08：并合成员级势能记账与残骸 changeover 同步更新

- 接续中断补丁：增加事件作用域、线程局部回调，逐对结算前身/残骸的簇内 hard/soft 势能变化，包括不同半径之间的质量转移及外半径以外常数项；扣除 AR 内部和 Hermite CM 已记贡献。簇内 soft 补偿同时扣除原 `dm*pot_soft` 的对应贡献，避免重复记账。未修改 AR/TTL 能量参考和并合物理处方。
- 未分组正质量粒子按自身质量计算目标半径，稳定分组按组质量覆盖；目标在同步点提交，搜索前扩展范围，提交后重建软力和势能缓存。也允许原组成员半径缩回自身质量尺度。旧快照无需新增字段，按原有 mean_mass_inv 初始化规则恢复目标。
- 默认生产配置编译通过。新增连续并合、延迟半径应用、实际软力修正入口的总势能/合力闭合回归；完整 ASan/UBSan bugfix 回归覆盖带/不带 AR_CHANGEOVER。
- N2000 真实 21 粒子事件：hard 残差从 -77172292.148693472 降至 +45450.217162102461，绝对值减少约 99.941%；AR/H4/同步步数及最终粒子文件不变。这仅证明该 hard 事件的主要记账跳变得到修复，并非全局误差清零。
- 真实 N2000 data.26 的 16 树步重启，1 MPI rank×2 threads、2 ranks×2 threads 均通过，无 hard-energy 告警；该区间没有再次并合。未安装或重启正在运行的模拟。
- 限制：簇外软场仍是树步冻结近似，其他 AR 组内部坐标未额外同步；事件残余误差、长期收敛及 N4000 尚未验证。每个 hard 区间增加 O(N_hard) 状态备份；半径提交额外一次树重建，不能沿用此前“开销很小”的估计作为实测结论。
- 原始回放、重启、计时及详细说明：`diagnostics/merger_boundary_20261108/REPORT.md`。


### 新版编译安装（2026-11-09）

按用户要求编译并安装当前工作区新版到 `/home/joojie`，保留 MPI/OpenMP/AVX2、BSE/BHMERGER、TT3 与四极矩配置，包含并合势能记账与同步 changeover 更新。`petar` 指向 `petar.mpi.omp.avx2.bse.tt3rd`；主程序、hard.debug、format.transfer 的安装文件 SHA-256 均与构建产物一致，安装后的 `petar -h` 在沙箱外退出 0。未启动模拟。构建、安装、帮助输出及校验记录保存在 `diagnostics/petar_install_20261109/`。此前列明的长期收敛与性能验证限制仍保留。


## 2026-11-09：确认同步树重建的质量表示错误（尚未正式修复）

排查 SsIDvu/N2000 的 -12.3999% 全局误差，确认上一轮 changeover 补丁新增的树重建在 kick 恢复 hard 成员真实质量之后调用，而簇势能修正要求人工粒子替代成员的源质量仍为零。误用 `_pj.mass` 产生额外线性截断势能项，污染 pot_soft/pot_tot 和加速度；随后 dm*pot_soft 抬高参考能量，导致负误差跳变。

- 两份真实并合 dump 的质量/位置核算分别解释大跳变的 99.6926%、99.9063%，合计约 12.5544 个原初能量百分点。
- 固定状态使用真实生产力修正路径复现：总势能凭空增加 19915812.8008；恢复正确源质量约定后全部 soft/total 势缓存偏差为 0。
- 相同初态和种子的完整程序诊断 A/B，在共同同步时刻原逻辑总势能增加 29327108.3331，临时修正版为 0。之后轨道会分开，未声称逐位重现原并合或完成长期收敛验证。
- 仅在 `diagnostics/global_merger_audit_20261109/` 修改诊断副本；生产源码及安装程序尚未修复。后续需覆盖 kick、输出/writeback、group_data 转为 CM、MPI 连接簇及终止边界的生命周期，不能只补最终能量数值。
- 证据、命令、原始日志与限制见 `diagnostics/global_merger_audit_20261109/REPORT.md`。此前“修复完成”的说法应限定为已验证的局部 hard 事件；全局路径存在本项已确认缺陷。


## 2026-11-09：同步重建源质量表示正式修复（源码及工作区构建）

- 引入 TreeSourceMassGuard：重建期间将人工粒子替代的真实组成员设为零源质量，覆盖 soft 数组、本地及 MPI 连接 hard 数组；完成后恢复调用前质量。无人工粒子的双星不屏蔽。保留分组元数据与并合物理处方，不靠事后抵消能量数值修复。
- `rebuildSoftForceWithChangeover()` 统一执行受保护的软树和势能重建；退出边界先提交/重建，再将 group_data 转为 CM，修复字段解释次序。DM 耦合力在恢复真实质量之后计算。
- 新增真实 FDPS 树生命周期测试及 Makefile.in 构建入口。1 rank×1 thread、2 ranks×2 threads 的 kick 前、kick 后、输出写回后，势能/加速度偏差均为 0，原质量完整恢复；绕过 guard 的负对照捕获相对误差 24.4038，按预期退出 255。生产程序的双进程短程退出测试通过。
- 默认生产编译和完整 ASan/UBSan bugfix 回归通过。含一次 541+1549 并合的 N2000 短程重启通过，无 hard-energy 告警；最后残差 -0.22982556%（相对该重启初始 |E|），最大单输出增量 0.06691068 个百分点。日志存在原有自旋进动失败回退；不是原历史轨道逐位重放，也不是长期全局误差清零证明。
- 本次尚未安装；`/home/joojie/bin/petar` 仍是上一版。详细报告、命令、退出码与校验值见 `diagnostics/changeover_source_fix_20261109/REPORT.md`。


### 同步源质量修复版安装（2026-11-09）

按用户要求将上述修复版安装到 `/home/joojie`，保留 MPI/OpenMP/AVX2、BSE/BHMERGER、TT3/quad 配置。`petar` 指向 `petar.mpi.omp.avx2.bse.tt3rd`；安装主程序 SHA-256 为 `fa82d91970512aeed4428b7c359c14bef5cbaf409d3083450e3c3efbd3565d23`，与工作区构建产物一致，hard.debug 与 format.transfer 亦校验一致。安装后的 `petar -h` 验证通过。未启动模拟。记录：`diagnostics/changeover_source_fix_20261109/install/`。此前“尚未安装”的状态由本条更新；长期收敛验证限制仍保留。


## 2026-11-09：Hzkiqy 人造轨道粒子 changeover 提交修复

- 确认真实成员提交新半径后，人造轨道源遗漏 `updateWithRScale()`，导致重建软树与hard使用不同分割。同步边界现将人造待提交半径计入MPI全局重建判断，并在重建前统一提交。保留源质量guard，不修改能量参考或并合物理。
- 新增真实FDPS边界的目标半径断言和恢复旧半径的负对照；1 rank×1 thread、2 ranks×2 threads均通过，负对照按预期退出134。生产构建、完整ASan/UBSan bugfix回归及独立生产短跑通过。
- 更正上一轮诊断：`isNextStart()` 在结束kick取步长后不是速度同步标志，应使用 `stat.time-time_kick`。旧“解释48%”比例失效；同历史单步从-0.03226649降至+0.00001752个原初能量百分点，绝对增量下降99.9457%。原data.status及演化图不受此诊断问题影响。
- 剩余：修复版约0.5 trh短跑仍有-0.0900个百分点全局残差。同步追踪将最大剩余跳变定位到新组形成、旧drift关闭的首half-kick，早于本次修复的提交后路径；其源/CM/潮汐旧分割一致性还需因果验证。不能宣称整个累计能量问题已修完。
- 本次仅更新源码及工作区构建，未安装、未改动原模拟。证据与限制见 `diagnostics/Hzkiqy_changeover_fix_20261109/REPORT.md`。

### 人造轨道半径修复版安装与区间重跑（2026-11-09）

按用户要求完成安装，`/home/joojie/bin/petar`主程序哈希为 `d7de30fe453f743844446ff724fb07e8481359bfa7425d85f72df17f81d672f8`，与构建一致，帮助启动通过。旧安装程序备份在 `diagnostics/Hzkiqy_installed_replay_20261109/petar.before_install`。

同一data.76→77区间、1rank×4threads、原输出节奏：旧版重跑末残差-0.18257555%、新安装版-0.08994678%（均按原模拟初始能量）；新版按重启能量为-0.27017359%。密集输出诊断峰值为原初能量的0.07562475%，仍存在跳变。三轮退出0，首末时刻匹配，独立快照能量复核通过。确认局部修复已安装，不能宣称全部能量误差修完；事件历史和输出节奏差异的限制见该目录REPORT.md。

原Hzkiqy主运行已于22:32停止（status退出143），最后177.988918 trh0，旧版累计-8.12514645%；本次未续跑该轨迹或启动完整长程模拟。

## 2026-11-10：新版残余首half-kick误差定位（尚未生产修复）

确认quad模式的三个等质量伪粒子仍在 `hard.hpp:2935` 按j%2复制真实成员的旧半径，沿用了与轨道采样节点有关的成员对应关系。对30+1028.732 Msun双星，三个352.911 Msun节点采用“小、大、小”半径，导致小半径soft分量对应705.821而非30 Msun。上一项提交修复在首half-kick之后生效，不能消除已施加的错误冲量。

两个新组形成的固定状态对照中，保持伪粒子位置、每位置总质量和四极矩不变，只按真实成员质量混合旧cutoff权重，残余跳变从-0.0172375、-0.0129636降到+0.0001458、-0.0000725个原初能量百分点。用原伪源逐对重建与树力的能量差仅约1e-8个百分点；无并合或真实质量变化。主体落在组外邻星ID1693、950的soft kick，不能归因于树开角或漏扣并合/逃逸能量。

仅诊断验证，安装程序及生产源码未改；正式修复还需处理旧分割、CM/潮汐和MPI表示的一致性，不能把局部对照当作全程收敛。完整证据见 `diagnostics/Hzkiqy_residual_force_20261109/REPORT.md`。

## 2026-11-10：修复 quad 旧截断半径的源质量权重

- 三个四极矩位置各拆成 m1/3、m2/3 两个共点源，匹配 hard.hpp 的交替成员旧半径赋值；保留每个位置总质量、质心与四极矩，提交后仍等价于原源分布。
- 处理共点人工源在软力、势能、changeover 和四阶力梯度修正中的有限极限；四阶局部源参数从无 acc 的 Ptcl 修正为实际 FPSoft。
- 增加质量/解析二阶矩/共点源回归及真实分组边界旧半径承载质量检查；单/双 MPI 边界、完整 ASAN/UBSAN bugfix 套件通过；旧三节点负对照按预期失败。
- 固定状态两点误差从 −0.017237518%、−0.012963574% 降至 +0.000145773%、−0.000072459% 原始 |E0|，局部降低 99.15%、99.44%。
- 同一修复二进制原频率重跑两次末端分别 −0.00416724%、−0.04967395% 原始 |E0|；并合前相同，反冲随机速度不同，不能把第一次数值当作稳定保证。密集输出末端 +0.01824325%、峰值 0.02741386%；2 MPI × 2 threads 末端 +0.00913614%。仍有其他残余误差。
- 已安装到 /home/joojie/bin，main SHA256 a260ead926e9e2ff06e4998887445210790079db103f0c70fc2b8d5add2c83d6；旧程序已备份。普通重启快照兼容，旧人工粒子布局的 hard 调试转储需匹配旧程序，新 hard.debug 会拒绝不匹配布局。
- 完整证据：diagnostics/Hzkiqy_quad_cutoff_fix_20261110/REPORT.md。

## 2026-11-10：非并合轨道更新的 hard/soft 事件记账与小分母警告保护

### 问题与修改

- DMlRMi/N2000 的 6 粒子事件中，1757+1669 在并合前更新轨道，质量不变但位置重建；AR 已计入内部能量变化，而 Hermite 的 `pot*dm` 为零。原事件补偿仅覆盖真正并合，遗漏与组外粒子的势能变化。三次更新缺失 hard 项合计约 −17031.927，主要来自邻星 59。
- `ar_interaction.hpp` 在两成员质量为正的 BSE `postProcess` 开始前保存事件状态；`merger_partition.hpp` 在质量或位置实际改变时结算，包括非并合更新。纯速度变化仍由现有 AR/CM 动能记账处理，无状态变化不调用结算。
- 复用 `finishMergerPartition` 的逐对 hard/soft 分区计算：扣除 AR 已计入的内部势能变化、Hermite CM 质量项和边界已计入的 soft 质量项，避免重复记账。soft 补偿仅进入全局参考能量，不进入 hard 残差或 AR 时间变换参考。
- 非并合更新不重新安排 changeover 半径或修改搜索半径；实际并合仍沿用原半径更新规则。增加事件头文件的 Makefile 依赖，确保增量编译能感知修改。
- 两处 hard-energy 警告共用 `HardEnergyWarning`，分母为 `max(1, |E_SD|, 1e-3*(|K_SD|+|U_SD|))`。保留 `--energy-err-hard` 阈值、原始残差和累计误差，只降低净能量相消时的警告敏感性。NaN/Inf 仍视为异常；分量先乘 1e-3 再求和，避免尺度计算溢出。
- 中间调试检查使用已扣除事件补偿的残差及修正参考能量；最终检查沿用当前净能量。警告补充 normalized_error、energy_scale、cancellation_floor；hard.debug 对每次终态输出检查尺度，即使警告被抑制仍可审计。README 已解释新判据。

### 验证

- 默认主程序与 hard.debug 完整构建成功。
- 完整 `test/bugfix/run.sh` 回归退出 0；ASan/UBSan 默认及 `AR_CHANGEOVER` 两种集成配置通过，包括新增非并合记账与小分母检查。日志中的 merger log 写入失败是既有故障注入测试的预期输出。
- 原 7 个 hard dump 按生产阈值 1e-3 回放均成功，过程和最终检查均无 hard-energy 警告。最大事件残差从 −16576.7 降至 +455.273，归一化误差 3.23199e-5，此事件没有启用小分母下限；其改善来自事件记账。
- 其余 6 个事件的原始残差不变，仅由小分母保护消除警告。7 个回放在可用输出精度下，最终动能、势能及 H4/AR/同步步数均与旧回放一致；不是粒子二进制终态逐位比较。
- 新增回归覆盖质量与质心不变的非并合位置重建、组外势能、hard/soft 分配及不更改待更新半径；原连续并合测试保留。警告测试覆盖近零净能量、真实大残差、零能量、NaN/Inf 和尺度溢出。

### 范围与产物

- 正式源码和工作区构建产物已更新；未覆盖 `/home/joojie/bin` 安装程序，未重启现有模拟。已有未提交修改保留。
- 仍有积分残差；本次未验证长期全局能量收敛，也未改变/验证大幅瞬时轨道重建的物理近似。新增事件快照及逐对结算有额外开销，尚未做长期性能基准。
- 本次构建、测试和 7 个回放日志：`diagnostics/DMlRMi_energy_fix_20261110/`；前期因果分析：`diagnostics/DMlRMi_hard_audit_20261110/REPORT.md`。

## 2026-11-10：SD 残差定位及紧双星局部坐标势能修复

### Debug 定位记录

- 对象为 DMlRMi/N2000 `hard_large_energy.0.3.3.1794281331`：11 粒子、1 个 AR 双星组，原 `dE_SD≈−1.29338`。关闭恒星演化结果不变，排除此前 BSE 非并合事件遗漏项。
- AR slowdown 因子更新前后残差连续：基线 30 次、收紧参数 398 次更新，诊断输出的残差跳变量均为零。最终 AR 内部 SD 残差约 7e-11，远小于 hard 层的单位量级残差。
- 根因：成员坐标从 AR 局部系加回质心后，Hermite 又用两原点坐标相减计算极小间距，丢失相对位置精度。基线同一终态局部间距 9.7098568545044537e-11 pc，写回坐标得到 9.7098568542023842e-11 pc；相对误差约 −3.11e-11 被约 −4.17e10 的内部势能放大，产生 +1.2970962524 的势能差。
- 诊断只替换能量计算中的组内距离，SD 残差由 −1.2933319462 降至 +0.0037643063，步数、动能及 SD 修正累计一致。该事件主要是能量测量的坐标相消误差，不是 SD 因子漏账或 AR 内部积分失控。

### 正式修改

- `hermite_interaction.hpp::calcEnergy` 为活动 AR 成员建立 group ID 和局部位置索引。在原有唯一一轮粒子对循环中，同组 pair 直接使用 AR 局部坐标差；组外 pair 保持原点坐标差。
- 保留全局 pair 顺序、质量、hard changeover、软化及外部扰动能量计算；不改变动力学力、AR 能量参考、粒子状态、步长控制或 SD 因子。
- 索引为线程局部可复用 vector，避免共享 manager 下的数据竞争和反复分配；不改变参数二进制布局。无 AR 组时不初始化索引。没有重复求势能、没有新增 pair 数、开方次数或 changeover 求值次数。新增成本为 O(N) 索引初始化及每对的组归属判断；缓冲区保留每线程峰值容量。

### 验证与性能

- 默认主程序和 hard.debug 完整构建成功。新增测试覆盖紧双星共同平移、多个活动组/稀疏组槽位/非连续成员编号、三成员组、组外单星、不同 cutoff、非零软化及取消分组后不使用旧索引。
- 正式 hard.debug：基线 `dE_SD` 从 −1.29338 降至 +0.00376429，H4/AR/同步步数仍为 62/104/62；eta4=0.005、AR 容差 1e-12、ds_scale=0.03 时，从 +0.342675 降至 +3.78317e-7，步数仍为 300/896/383。诊断 O2 与正式调试 O0 的末位和步数有差别，各对照均采用相同构建。
- 原 7 个 hard dump 全部回放退出 0，过程及最终 hard-energy 警告均为 0。
- 完整 `test/bugfix/run.sh` 退出 0，默认和 AR_CHANGEOVER 两种配置的 ASan/UBSan 集成回归及新增局部坐标测试通过；日志中的 merger log 写入失败属于既有故障注入测试的预期输出。
- O3 能量内核微基准（5 次交替测量取中位数）：11 粒子单组 0.246722→0.263251 微秒/次，约 +6.7%；64 粒子单组 6.62976→6.79666 微秒/次，约 +2.5%；256 粒子单组 99.4177→103.591 微秒/次，约 +4.2%。其他测试配置耗时比约 0.91–1.10。测量有系统噪声，使用简化粒子存储且不含外部潮汐计算，仅评估该内核，不代表整个模拟耗时；尚未做端到端长期性能基准。
- 剩余物理能量变化包含 slowdown 动力学效应，不能把物理能量残差直接当作 SD Hamiltonian 积分误差；本次未证明长期全局能量收敛。

### 产物与生效范围

- 定位报告和逐阶段日志：`diagnostics/DMlRMi_sd_audit_20261110/REPORT.md`。
- 正式构建/回归/回放/性能日志及微基准源码：`diagnostics/DMlRMi_sd_coordinate_fix_20261110/`。
- 正式源码及工作区构建产物已更新；未安装到 `/home/joojie/bin`，未重启运行中的模拟。原有未提交修改保留。

## 2026-11-10 13:21（北京时间）：安装最新版能量修复

- 已将通过回归和事件回放验证的主程序、hard.debug 及格式转换程序安装到 `/home/joojie/bin`；`petar` 入口指向最新 `petar.mpi.omp.avx2.bse.tt3rd`。
- 包含非并合轨道更新 hard/soft 记账、小分母警告保护及 SD 局部坐标势能修复。
- 采用临时文件加原子替换，旧版备份至 `diagnostics/petar_install_20261110_energy_fixes/backup_20261110_132111/`。
- 安装后主程序 SHA256：`4db04db603d9ea6d97b84644854ede28276d0fba11a5c59b393ff4d8333ba9a7`；hard.debug：`b17c15272ec184dfbcec0bbee2b749ef9fc9fcbc3cd4a5994e2f5670f10a2a3f`，均与已验证构建相同。
- 主程序及 hard.debug 的 `-h` 启动检查均退出 0。未停止或重启正在运行的任务；已启动进程继续使用其原有二进制，新启动进程使用新版。
- 安装清单及启动日志：`diagnostics/petar_install_20261110_energy_fixes/`。此前“未安装”描述为当时状态，本条记录后续安装结果。

## 2026-11-10：GOlAQ0 能量跳变定位到 AR 邻组扰动力表示不一致（仅诊断验证）

- 20.50–21.50 trh0物理累计残差增加0.827429个百分点；hard累计占约90%。四个116–127粒子大残差转储单线程复现，合计0.628549个初始能量百分点，约占总跳变76%。
- 确认 `ARInteraction::calcAccPert` 对邻组仅用质心位置，而Hermite已解析同组成员。近距离组间相遇时，质心运动与内部运动使用不同源分布，产生实际非守恒做功，不是单纯能量记账。
- 首119粒子事件收紧Hermite10倍、收紧AR、关闭soft潮汐或恒星演化均保留约+5.62e6残差。四份回放内真实粒子changeover半径均未改变，不能归因于本区间半径提交瞬间。
- 仅在diagnostics副本将AR邻组源展开到成员位置，四项残差+5.61950e6、+3.12304e6、+5.92897e6、+9.81155e6分别降到+67.1554、-131.713、+115.839、-490.210，均下降超过99.99%。独立积分力差做功解释原残差至约0.07%–2.6%的差异。
- 尚未修改/安装生产程序。诊断实现未完成异步成员预测、slowdown和跨MPI等正式验证；全局剩余约24%跳变尚未完全定位，旧cnpcJS逃逸及与最近changeover补丁的因果关联仍未证实。
- 证据、独立副本、回放日志与可复现脚本：`diagnostics/GOlAQ0_energy_jump_20261110/REPORT.md`。

### AR邻组成员解析诊断副本性能初测

四个强遭遇转储均以-O3、单线程、关闭密集逐步打印测试，预热后各7次，积分区间CPU中位数分别增加8.9%、-0.7%、26.2%、35.0%，四例中位数之和增加19.1%。计时包含求力成本和轨迹/步数变化；保留布局兼容宏与调试断言，不是全程序生产性能测试。未测正式异步成员预测实现或整场N2000减速。详见 `diagnostics/GOlAQ0_energy_jump_20261110/PERFORMANCE.md`。

## 2026-11-10：AR邻组成员扰动力正式修复并安装

- `ARInteraction::calcAccPert` 的邻组源改为逐成员位置、质量和changeover权重，保持Hermite与AR的空间源表示一致。预测每调用计算一次，线程局部缓冲区复用容量。
- SDAR邻居记录绑定源积分器的只读预测回调，正常邻居生成及新组插入均覆盖。局部成员按源AR自身时间和slowdown后的速度线性预测到请求时刻，支持向前/向后预测；质心继续使用Hermite预测。嵌套TREE与ARRAY分支沿用各自速度分解，文件格式未变。
- 新回归额外发现无邻居时旧扰动力残留，补齐该分支acc_pert/pot_pert清零。此前其他工作区改动保留。
- 四个生产转储的百万量级残差分别降至-402.348、-311.428、-164.104、-665.526；归一化约4.97e-8至2.41e-7。首事件/最大事件收紧Hermite10倍进一步降至-23.2224/-4.51925，支持剩余误差的收敛性。
- 完整ASan/UBSan bugfix脚本通过（AR_CHANGEOVER开/关两种集成测试）；新增回归覆盖AR/CM异步时刻、反向预测、slowdown=4、直接成员力、质心净扰动、邻居回调复制及4线程缓存隔离/清空。1 MPI×4线程、2 MPI×2线程生产短跑退出0、无hard告警。短跑并非全场N2000验证。
- 2026-11-10 14:41 CST已备份并原子替换`/home/joojie/bin/petar.mpi.omp.avx2.bse.tt3rd`和`petar.hard.debug`；format.transfer校验一致。`petar`入口指向新版，主程序SHA256：`63fde223f289a3c157b91bbaf9e405a2adfa856b063562947b192277c6e2e648`。main与hard.debug帮助启动退出0；format.transfer帮助按源码返回255，属于既有行为。
- 备份：`diagnostics/GOlAQ0_group_force_fix_20261110/backup_20261110_144116/`。报告、补丁、回放、测试、安装哈希见同目录REPORT.md、validation.json、installation.json。此次没有停止或重启已有主模拟；已运行进程仍使用原二进制。
- 限制：成员预测仍为时间离散近似；未证明原全局跳变剩余部分全部消除，未完成整场性能测量或将旧cnpcJS的逃逸归因于本缺陷。

## 2026-11-10：x8ymxY/N2000 hard 修复后 soft 耦合与记账诊断

- 原始 data.47→48（23.49854→23.99851 trh0）物理残差增加 3858331.9714，即 0.099054216 个百分点；hard 累计变化仅 −1998.9559。全粒子对重算得到缓存势能误差变化只有 729.5446；四个逃逸粒子的删除能量重建与总参考修正变化相差 −1646.1214，无同量级漏记账证据。
- 当前源码单线程重启回放捕获两个大单步事件：直接物理 Hamiltonian 分解中，hard 漂移误差仅 +4.6148、+1.6443，而 soft 踢步做功加 soft 势能变化分别为 +283666.6618、−183723.2771；闭合误差小于 5e−6。两步真实质量和 changeover 半径未变、无人工粒子，参考修正接近零。固定状态改用直接 soft 力仅改变端点虚拟踢步能量约 0.08–1.8。
- 正向事件中，30 M☉ 的 ID 1290 在 400.2178 M☉ 的 ID 1754 近旁快速运动，瞬时 r/v 约为树步长的 1/30。hard 可解析内部运动，但外层 soft 端点踢步未充分解析沿快速轨迹的做功；这属于时间分裂误差，不能通过给 etot_ref 补账消除。
- 固定同一同步局部初态、相同区间，只缩小 soft 步长：直接物理能量变化 h 时 +283670.93，h/8 时 +12882.64，h/32 时 +812.80，h/64 时 −65.24，h/128 时 −33.69，h/256 时 −8.12；hard 残差一直约 3–5。中间变号，不能宣称全范围严格单调二阶收敛；细步结果支持积分时间分辨率不足。h/256 仅作诊断，不是全场推荐参数。
- 与上次 AR 邻组空间源修复不同，此次主导项在外层 soft–hard 时间耦合。未证明由上次修复新引入，也不能将穿越固定 changeover 区间误称为半径提交漏账。重启轨迹与原运行不同，尚未逐步重现原始全部跳变或完成全场误差归因。
- 此次仅增加诊断副本/脚本/报告，未更改或安装生产程序，未重启主模拟。报告及原始证据：`diagnostics/x8ymxY_soft_coupling_20261110/REPORT.md`；闭合与收敛数据分别为 `split_audit.json`、`convergence.json`。

### x8ymxY 正向大残差：slowdown 与动态分组补查

- 对同一局部树步加只读探针，相关 hard cluster 有 11 成员，最大误差对的邻星已全部包含。hard cluster 成员间依然按 changeover 分配力；1147–1290 的 soft 力占比约 93.84%→100%，因此 cluster 大小不是 hard 积分覆盖完整力的保证。
- 树步内短暂形成 1754–1290 两成员 AR 组，17 个可见 Hermite 边界采样中 slowdown κ 全为 1，9 个 hard 单星扰动邻居，soft 潮汐指针为空、soft_pert_min=0。人工粒子为零不代表树步内没有 AR，profile 初始 AR_group_N 也不统计所有动态组。
- 代码支持 slowdown 微扰判据及 soft 潮汐，但该动态组没有预建人工粒子的相容张量可继承；soft 仍由外层踢步承担。不能将此直接称为全局漏力/漏账，或据此强制 slowdown。未做 hard 半径控制试验，不能断言 hard 覆盖半径有缺陷。
- 仅诊断副本，能量结果与无分组探针一致。记录见同目录 REPORT.md 补充节、group_audit.json 和 group_dt_div1/stderr。

### x8ymxY 潮汐耦合影响尺度补查

- 本事件成员接受各自的 soft 踢步，故 soft 张量为空不等于外部潮汐物理作用被删除；直接加内部张量而保留完整差分踢步会重复计力。正式人工组采用共同 CM 踢步加内部潮汐的不同分工，不能只复制张量指针。
- 按实际 changeover 对两个端点直接求和，1754–1290 的 soft 相对加速度/内部牛顿相对加速度约 0.0147% 和 1.13%。成员 cutoff 不同，这个 soft 分力差不等于共同光滑场的纯潮汐梯度，也不是轨道变化的预测。
- 当前证据证明外层时间分裂需改善，未证明必须采用张量或补张量即可解决整步残差。未做一致的张量有无对照，不能推断全场并合/逃逸变化。脚本和数值见 `diagnostics/x8ymxY_soft_coupling_20261110/soft_relative_force_scale.py/json`，报告已补充；生产源码未改。
