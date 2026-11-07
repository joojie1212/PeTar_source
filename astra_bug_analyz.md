# astra bug 登记册当前源码复核

复核日期：2026-11-07（Asia/Shanghai）。输入：[astra_bug_registry.md](astra_bug_registry.md)。范围：该表正文的 **192 个稳定编号**，逐项给出当前状态和维护模块；初次审计未修改源码；后续按用户授权修复，当前状态已在逐项表及第 8 节更新，不重新编号原登记册。

**修复前审计结论：124 项仍保留原问题机制，4 项部分修复仍有残留，38 项原机制已修复，3 项在当前实现中不成立；另有 18 项待确认、4 项不宜作为缺陷、1 项重启契约缺口。** “仍存在”包括非默认宏、辅助接口、退化输入和无生产调用者的潜伏问题，绝不表示 124 项都在当前模拟中触发。124 也不是去重后的独立根因数。

## 1. 基线、判断口径和归属

本次以 **当前工作区** 为准：Git HEAD `5b9a35a`（`debug on slowdown`），包括用户已有未提交修改，特别是 B44/B45 的并合时钟和逃逸删除实现。登记册声称的旧基线 `fbebd9c` 在本仓库对象库中不存在，不能进行完整基线 diff；采用逐项源码核对、现存修复提交 `e4f066a`、[修复汇总](BUG_FIX_SUMMARY.md) 和有限回归测试交叉判断。修复汇总是线索，不替代当前源码。旧报告中的速度差、Q、Valgrind/ASan 报错、发生比例和生产配置并非本次复现实测，本文不把它们作为新的运行结果。

当前 [Makefile](PeTar_new_version/Makefile) 配置为 MPI+OpenMP、x86 SIMD/F32、BSE、BHMERGER、KDK、AR_SLOWDOWN_TREE、ext_mode=off、use_quad=no；CLUSTER_VELOCITY、HARD_DUMP、AR_DEBUG_DUMP 开启。该配置只代表工作区构建，不证明任意外部安装二进制或历史模拟采用相同宏。原登记册中的 PETAR_INNER_THREADS/team/task 代码不能套到当前串行 Hermite 组循环上。

| 状态 | 含义 |
|---|---|
| 仍存在 | 当前仍有对应错误表达式、状态缺失或接口风险；是否生产可达看逐项条件，部分只完成静态确认 |
| 已修复 | 原条目所指机制已被消除；不保证所在模块没有其它问题，也不表示每项都做了本次动态复现 |
| 部分修复 | 原条目的某一子问题已消除，但行为/状态链尚未闭合 |
| 当前不成立 | 原报告依赖的代码/并行路径已不存在；不等于已实现一般性并行安全方案 |
| 待确认 | 需要模型契约、合法输入可达性或独立数值证据，现有材料不足以判为确定 bug |
| 非缺陷 | 属正常行为、冗余代码或额外工作量；不列入修复队列 |
| 契约缺口 | 逐事件连续重启等承诺未实现，是否阻断依赖项目约定 |

模块按实现及修复责任划分，不按 ASTRA 的 AREA 字母猜测：

- **SDAR**：`SDAR/src/AR`、`Hermite`、`Common`；如 `TR-03/04/10` 实际属于 SDAR。
- **FDPS**：`FDPS/src`；包括 `TR-02` 的 ID 映射。`DG-01` 主修复点在 FDPS，另涉及已移除的 SDAR 并行路径。
- **PeTar**：`PeTar_new_version` 的主程序、软力/SIMD、AR/Hermite 物理适配、BSE Fortran 与接口、BH/precession。`H4-06/07` 实际在 PeTar，`BL-06` 是 PeTar 对可选 SDAR 接口的适配缺口。
- **PeTar 配套（mcluster）**：为满足三大模块口径归入 PeTar 生态，但它位于独立的 `mcluster/`，并非 PeTar 积分器代码；全部 28 项保留这一标识。

跨模块项只按主责计数一次：`AR-01` 主责 SDAR、PeTar 回调联动；`AR-02/03` 主责 PeTar、SDAR 消费状态；`TR-07` 涉及 PeTar 与 SDAR 两处；`BL-06` 主责 PeTar、SDAR 泛型接口联动。

## 2. 状态统计

| 主责模块 | 仍存在 | 部分修复 | 已修复 | 当前不成立 | 待确认 | 非缺陷 | 契约缺口 | 合计 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SDAR | 6 | 0 | 10 | 3 | 0 | 0 | 0 | 19 |
| FDPS | 23 | 0 | 2 | 0 | 1 | 0 | 0 | 26 |
| PeTar | 74 | 4 | 26 | 0 | 12 | 2 | 1 | 119 |
| PeTar 配套（mcluster） | 21 | 0 | 0 | 0 | 5 | 2 | 0 | 28 |
| 合计 | 124 | 4 | 38 | 3 | 18 | 4 | 1 | 192 |

如将 mcluster 合并，PeTar 总数为上表两个 PeTar 行之和（147 项）。统计保留原编号：`BL-04/BL-09` 同一格式校验根因，`BL-10/DG-12` 部分重叠，`MC-16/MC-28` 是构建兼容问题的两个入口；不要按行数估算独立补丁数。

## 3. 修复前的优先级及残留问题（当前进展见第 8 节）

| 优先级 | 条目 | 判断及下一步 |
|---|---|---|
| 高 | TR-11（PeTar） | [treeNeighborSearch](PeTar_new_version/src/petar.hpp:971) 仍把 G=1 传入速度/近心判据；应统一运行时 G，并用单位变换一致性与连通分量回归验收 |
| 高 | IO-01（PeTar，部分修复） | [FPSoft ASCII](PeTar_new_version/src/soft_ptcl.hpp:112) 保存 deadline，但不保存动力学检查点；[checkBHMergerDynamicalChange](PeTar_new_version/src/ar_interaction.hpp:171) 因 record_time<0 重新预测。需整条重启—首次事件调用回归，单做 write/read 不够 |
| 高 | EN-06（PeTar） | [局部 mass_modify_list](PeTar_new_version/src/petar.hpp:1989) 遮蔽成员，连接簇中断后的质量变化可能漏入全局能量修正；需要多 rank 中断场景 |
| 高，恒星演化条件 | BSE-07、BSI-01、BSI-03 | 并合后恒星参数未刷新、接触潮汐因子量纲错误、CE 超新星后偏心率被旧轨道量覆盖；优先做实际 Fortran 事件回归 |
| 高，使用该初态链时 | MC-01/02/03/04/05/06/07/08 | 单位、Q、双星预演化和归一化问题仍在；须按实际 mcluster 命令及下游 petar.init 参数判断，不能推断每个生产初态都受污染 |
| 中 | BH-08/09、BH-05/10（部分修复） | 圆轨道负舍入、模块路径、返回形状和奇点恢复仍可能改变拟合分支；区别“安全回退”和“物理进动已正确完成” |
| 中 | FD-01、FD-19 | 析构未初始化指针和 TimeProfile 漏初始化仍在；FDPS 原计时竞争修复不覆盖这些问题 |
| 按启用条件 | AR-08/09/10、H4-05/06、TR-17/19、FD 系列潜伏项 | 先确认 ARRAY、KDKDK4、QUAD、大规模分块、复用列表或特殊 MPI/边界配置，再安排相应回归 |

上述优先级针对本次代码审计，不代表已经识别出某条现存模拟异常的原因。

### 四项不能简单写“已修复”

1. **IO-01**：修复汇总 B07 的“保存 deadline”仍成立，但 B44 后新检查点缺失会导致首次动力学检查重建 deadline。既有 ASCII 测试只核对 `time_merger` 往返，不检查恢复后的事件时刻，不能封闭此项。直接用 `ParticleBase::writeAscii` 与真实快照 `FPSoft::writeAscii` 也不能混淆。
2. **BSE-06**：数组零初始化解决 UB，不能保证 CE/MIX 后包层参数是新恒星结构的正确值。
3. **BH-05**：形状检查解决越界，未证明 `(1,1)` 圆轨道输出会正确演化，而非安全回退。
4. **BH-10**：近共振恢复仅覆盖验证为根合并的子集，不能代替整个 precession 库的端点/零自旋修复。

### 原表中应调整的判断

- `AR-01/H4-01/AR-11` 的并行前提在当前实现不成立；不能继续沿用 D1/D2 的生产可达结论。
- `BL-05` 的 G=1 默认、`MC-12` 的共同质心平移，以及 `BSI-02/04/11/12` 等模型相关项需要契约/公式证据。特别是经历失质量的恒星，等效初始质量大于当前质量不自动构成 bug。
- `GEN-02` 的无调用者、`MC-25` 改变抽样次数后随机流变化、`MC-26` 冗余守卫，本次不作为缺陷；`TR-09` 维持原非缺陷结论。
- `MC-19` 的 `i==N` 访问可能越界，`MC-20` 为缓冲区安全问题；不能只凭原 D4 标签断言只影响打印。`FD-20` 若用于对应周期力计算，也不是纯诊断问题。
- `IO-08` 覆盖输入还需要输出前缀/路径一致；`BL-12` 的0输入通常下溢后返回0，真正确定不终止的是正无穷分支，后续零步长停滞另论。

## 4. 逐项结论

下表覆盖全部192项。**原等级仅用于对照，不是本次重新授予的严重度**；以“当前状态/依据与触发条件”为准。源码链接指向当前文件，不复用已漂移的旧行号。复杂数值结论若未重跑，已经在说明中限定。

### SDAR

| 编号 | 原等级 | 问题 | 当前状态 | 当前依据与条件（含源码） |
|---|---|---|---|---|
| `ASTRA-AR-01` | D1 | AR 组并行下，他组质量事件对扰动源的可见性依赖调度次序 | **当前不成立** | [hermite_integrator.h](SDAR/src/Hermite/hermite_integrator.h)：当前 integrateGroupsOneStep 是普通 for 循环；原 task/team 组级并行实现不存在。跨 cluster 并行仍在，但不是原报告的同一共享扰动源路径。 |
| `ASTRA-AR-04` | D4 | destroy 分支漏记 slowdown 事件分类项 | **已修复** | [symplectic_integrator.h](SDAR/src/AR/symplectic_integrator.h)：destroy 结算已同时扣减 de_sd_change_cum_ 与 de_sd_change_interrupt_。 |
| `ASTRA-AR-06` | D3 | `AR_TIME_FUNCTION_MULTI_R` / `AR_SLOWDOWN_MASSRATIO` 分支不可编译 | **已修复** | [symplectic_integrator.h](SDAR/src/AR/symplectic_integrator.h)：_bin、calcInvR、gt_drift_inv、calcGtDriftInvIter 和 binary_slowdown[i]->mass 均已修正；不等于已验证这些非默认公式的全部物理精度。 |
| `ASTRA-AR-07` | D3 | `SymplecticStep` 写裸指针，读回前先清零 size 跳过 `clear()` | **已修复** | [symplectic_step.h](SDAR/src/AR/symplectic_step.h)：SymplecticStep 保存 SYS2/order/type，读取先 clear 再重建；核心 sanitizer 往返回归通过。 |
| `ASTRA-AR-08` | D3 | `AR_SLOWDOWN_ARRAY` 下合并后单存活粒子直接漂移缺少 `1/κ` | **已修复（专项回归）** | [symplectic_integrator.h](SDAR/src/AR/symplectic_integrator.h)：经逐项批准，ARRAY 单存活粒子漂移使用 dt/κ，物理时间仍增加 dt。实际 integrator 中注入并合事件，κ=1/8 时的位移分别为 (2,4,6)/(0.25,0.5,0.75)，结束时间均为 2。 |
| `ASTRA-AR-09` | D3 | 零质量成员追加为人工树根后轨道量/`r`/`stab` 未计算，复用上一棵树残值 | **已修复（专项回归）** | [information.h](SDAR/src/AR/information.h) [binary_tree.h](SDAR/src/Common/binary_tree.h) [symplectic_integrator.h](SDAR/src/AR/symplectic_integrator.h) [hermite_integrator.h](SDAR/src/Hermite/hermite_integrator.h)：经逐项批准，零质量人工节点清空轨道缓存，重算质量/质心/当前分离并标为不稳定；保留正质量子树。周期更新跳过人工节点，Hermite 拆组不再消费其虚假轨道。4 正质量→2 正+2 零、单存活、全零及重合节点回归通过；完整生产并合/拆组仍待长程验证。 |
| `ASTRA-AR-10` | D3 | 时间同步首个子步用绝对时间 `_time_end/time_table[k]` 缩放步长 | **已修改，待专项验证** | [symplectic_integrator.h](SDAR/src/AR/symplectic_integrator.h)：首个过冲步使用该步开始时间计算缩步比例。 |
| `ASTRA-AR-11` | D3 潜伏 | team 并行 AR 组循环的中断记录按临界区到达顺序选取 | **当前不成立** | [hermite_integrator.h](SDAR/src/Hermite/hermite_integrator.h)：原 team/task 分支已不在当前 Hermite 调度中；现为按索引串行积分，原临界区到达顺序选中断的机制不存在。 |
| `ASTRA-H4-01` | D2 | `initial_step_flag` 在力循环内部被写，跨 item 读取无保护 | **当前不成立** | [neighbor.h](SDAR/src/Hermite/neighbor.h)：Neighbor 中标志传播仍在，但 Hermite 力计算和校正当前无 OpenMP 并行循环；不能把普通跨对象读写直接定为数据竞争。 |
| `ASTRA-H4-02` | D3 | `table_check_mask` 按活组数分配、按全局槽位索引 | **已修复** | [hermite_integrator.h](SDAR/src/Hermite/hermite_integrator.h)：table_check_mask 已按 groups.getSize() 分配 vector，检查也遍历组槽位范围。 |
| `ASTRA-H4-03` | D2 | `calcNextDtLimit` 的 `bitmap==0` 导致 `while` 不终止 | **已修复** | [block_time_step.h](SDAR/src/Hermite/block_time_step.h)：ticks 有有限性/范围检查，bitmap==0 返回 dt_min；nextafter(dt_min,0) 核心回归通过。 |
| `ASTRA-H4-04` | D4 | `ParticleH4::print` 把 `acc1` 打成 `acc0`，数组按地址输出 | **已修复** | [hermite_particle.h](SDAR/src/Hermite/hermite_particle.h)：print 分别输出 acc0/acc1 三分量；本次集成测试到逃逸测试前已通过打印检查。 |
| `ASTRA-H4-05` | D3 潜伏 | `checkNewGroup` 的 `used_mask` 按 `n_particle+n_group` 定尺却用 `index_offset_group_` 索引 | **已修改，待专项验证** | [hermite_integrator.h](SDAR/src/Hermite/hermite_integrator.h)：used_mask 改用覆盖 index_offset_group_+n_group 的动态存储，消除旧栈数组长度不足。 |
| `ASTRA-TR-03` | D2 | `generateBinaryTree` 的 `bin_host[ki]` 先比较后判界 | **已修复** | [binary_tree.h](SDAR/src/Common/binary_tree.h)：两个 bin_host 循环均先检查 ki 边界再访问；本次 3–6 粒子 ASan/UBSan 回归通过，主责 SDAR。 |
| `ASTRA-TR-04` | D2 | `copyDataIter` 每子树只把 `inow` 减 1，右子树覆盖左子树 | **已修复** | [binary_tree.h](SDAR/src/Common/binary_tree.h)：copyDataIter 用 _i-n_iter 推进子树位置；本次 6 叶树拓扑及指针归属回归通过，主责 SDAR。 |
| `ASTRA-TR-10` | D3 | Kepler 方程 Newton 用 `sqrt(1-sin²)` 代替 `cos`，丢符号 | **已修复（专项回归）** | [binary_tree.h](SDAR/src/Common/binary_tree.h)：Kepler 求解使用有符号 cos 导数，并加入区间保护与浮点收敛判据；e=0..0.99 的 505 个样例核对方程残差。 |
| `ASTRA-IO-02` | D2 | `AR::Information` 把 8 字节 `ds` 按 `sizeof(int)` 写读 | **已修复** | [information.h](SDAR/src/AR/information.h)：SDAR Information 用 sizeof(ds) 加 ARI2 标记；核心非整数 ds 往返通过，旧损坏格式拒绝。 |
| `ASTRA-IO-05` | D3 潜伏 | 组积分器 `writeBinary` 漏写 `gt_kick_inv_`、sd 能量与 `binary_slowdown` | **仍存在** | [symplectic_integrator.h](SDAR/src/AR/symplectic_integrator.h)：SDAR 组积分器 writeBinary 仍仅写部分能量/gt_drift 等字段，不写 gt_kick_inv 和完整 slowdown 数组；独立组序列化接口潜伏，不是常规 PeTar 快照路径。 |
| `ASTRA-BL-02` | D3 | `COMM::List::operator=` 调用不存在的 `reserveMem(nmax,mode)` | **已修复** | [list.h](SDAR/src/Common/list.h)：SDAR List 赋值使用单参数 reserveMem 并保护自赋值；核心 copy/local 回归通过。 |

### FDPS

| 编号 | 原等级 | 问题 | 当前状态 | 当前依据与条件（含源码） |
|---|---|---|---|---|
| `ASTRA-TR-02` | D3 | `getEpjFromId` 用混合 `tp_glb_[i]` 下标过滤压缩 `epj_sorted_` | **已修复** | [tree_for_force_impl.hpp](FDPS/src/tree_for_force_impl.hpp)：getEpjFromId 不再按 tp_glb 的 SP 标志过滤压紧 EP 表；本次交错 ID 映射回归通过，主责为 FDPS。 |
| `ASTRA-FD-01` | 条件 D2 | `TreeForForce` 默认构造不初始化裸指针，析构却无条件 `delete[]` | **已修改，待专项验证** | [tree_for_force.hpp](FDPS/src/tree_for_force.hpp)：拥有的裸指针初始化为 nullptr，未 initialize 的析构路径不再读取未初始化指针。 |
| `ASTRA-FD-02` | D3 | `MomentMonopole::set()` 对零质量 cell 做 `0/0`，NaN 沿 `accumulate` 爬到根 | **已修改，待专项验证** | [tree.hpp](FDPS/src/tree.hpp)：零质量 MomentMonopole 不再执行 pos/mass。 |
| `ASTRA-FD-03` | D3 | `REUSE_LIST` 复用上一轮排序排列与 `tc_*`，却按新粒子数 resize，新增下标未初始化 | **部分修复** | REUSE_LIST 在粒子数与旧排列长度不符时明确中止，防止读取新增未初始化索引；相同数量但成员/顺序变化仍由调用方遵守复用契约。 |
| `ASTRA-FD-04` | D3 | θ=0 交互表分支只填 EP、不填 SP；且用 `tp_glb_` 位置当 `epj_sorted_` 下标 | **已修改，专项回归通过** | 零 opening-angle 列表按 tp 地址解码 EP/SP，并加入 SP；混合 SP+两 EP 的集成测试通过，多 rank 生产运行未重跑。 |
| `ASTRA-FD-05` | D3 | `calcForceDirectParallel` 镜像展开用 `epj[j*n_image+k]`，步长应为 `n_recv` | **已修改，待专项验证** | [tree_for_force_impl.hpp](FDPS/src/tree_for_force_impl.hpp) [particle_system.hpp](FDPS/src/particle_system.hpp)：周期镜像直接力缓冲索引修正为 j*n_recv+k。 |
| `ASTRA-FD-06` | D3 | `-DLOOP_TREE` 下 per-thread 数组维度硬编码 12，初始化上界却是 `omp_get_max_threads()` | **已修改，待专项验证** | [tree_for_force_impl_force.hpp](FDPS/src/tree_for_force_impl_force.hpp)：LOOP_TREE 的七组固定 12 线程工作区改为按 n_thread 动态分配。 |
| `ASTRA-FD-07` | D3 | 以"缓存的 `omp_get_max_threads()`"当实际 team 大小做切分/索引 | **仍存在** | [ps_defs.hpp](FDPS/src/ps_defs.hpp)：切片使用 Comm::getNumberOfThread 的缓存线程数，未统一改为实际 team 大小；OMP_DYNAMIC/THREAD_LIMIT 等条件下风险保留。 |
| `ASTRA-FD-08` | D4 | `RandomSampling` 的 `j` 抽自 `[0,n_src)` 而非 `[i,n_src)`，子集分布偏斜 | **已修改，待专项验证** | [ps_defs.hpp](FDPS/src/ps_defs.hpp)：Fisher–Yates 抽样每步使用区间 [i,n_src-1]，避免固定区间引入偏差。 |
| `ASTRA-FD-09` | D3 | `collectSampleParticle` 首次分支 `freeMem()` 把 `size_` 归零，二次采样写回旧样本 | **已修改，待专项验证** | [domain_info.hpp](FDPS/src/domain_info.hpp)：重采样前恢复已有本地样本的逻辑长度，再追加第二种粒子系统的样本。 |
| `ASTRA-FD-10` | D3 | `reserve()` 先更新 `capacity_` 再 `ReallocInner`，后者按新容量 memcpy 旧段 | **已修复（专项回归）** | [reallocatable_array.hpp](FDPS/src/reallocatable_array.hpp)：reserve 先用局部变量检查新容量，再交给 ReallocInner 扩容，保留旧容量供拷贝；新增保值回归。 |
| `ASTRA-FD-11` | D3 | `MemoryPool::clearAll()` 置 `top_=NULL` 但保留 `cap_`，之后 `alloc` 直接返回 NULL | **已修改，待专项验证** | [memory_pool.hpp](FDPS/src/memory_pool.hpp)：MemoryPool::clearAll 保留已分配底址，将 top 重置为 bottom，避免仍有容量但基址为空。 |
| `ASTRA-FD-12` | D3 | `CalcRank1D` 的 `while` 对 `domain[]` 无界搜索 | **已修改，待专项验证** | [util.hpp](FDPS/src/util.hpp)：CalcRank1D 沿单坐标方向搜索时检查 rank_axis 边界，越出该方向首末域则明确退出，避免先越界访问 domain。 |
| `ASTRA-FD-13` | D3 | `LinkCell` 在 `n_cnt==0` 时用 `adr[-1]` 作 cell 下标 | **已修改，待专项验证** | [tree_for_force_utils.hpp](FDPS/src/tree_for_force_utils.hpp)：访问 adr[n_cnt-1] 前显式检查 n_cnt>0，无候选时诊断退出。 |
| `ASTRA-FD-14` | D3 | 旧 `MakeInteractionList*`/`MakeIPGroup*` 变体混用两套下标空间（SP 未 `ClearMSB`、EP 用 tp 下标） | **部分修复** | 旧交互表变体修正 EP/SP 地址解码及跳过 SP 时的循环推进；MakeIPGroupLongGLBTreeCellAsIPGBox 的地址空间问题仍保留。 |
| `ASTRA-FD-15` | D3 | `FindScatterParticleLongCutoffP2P` 函数体被掏空，所有出参原样返回 | **部分修复** | 未实现的 LongCutoffP2P 辅助接口现在明确中止，避免静默返回空发送列表；完整 P2P 功能仍未实现。 |
| `ASTRA-FD-16` | D4 | `~TreeForForce` 漏删 `adr_epj/adr_spj/adr_ipg_for_force_` 与 MPI `status_` | **已修改，待专项验证** | [tree_for_force.hpp](FDPS/src/tree_for_force.hpp)：析构补齐 adr_* 工作区及 MPI status 数组释放。 |
| `ASTRA-FD-17` | D4 | `getMemSizeUsed()` 循环体用 `tmp =` 而非 `tmp +=`，且漏计多组数组 | **部分修复** | [tree_for_force_impl.hpp](FDPS/src/tree_for_force_impl.hpp)：线程工作区内存计数改为累加；尚未完成所有遗漏数组的完整核算，本项仍有残留。 |
| `ASTRA-FD-18` | D4 | `ParticleSystem::n_proc_exch_` 未初始化即 `+=` | **已修改，待专项验证** | [particle_system.hpp](FDPS/src/particle_system.hpp)：粒子交换进程计数 n_proc_exch_ 初始化为 0。 |
| `ASTRA-FD-19` | D4 | `TimeProfile` 五个成员构造/`clear()` 前未初始化，却每步参与 `operator+` | **已修复（专项回归）** | [ps_defs.hpp](FDPS/src/ps_defs.hpp)：TimeProfile 构造函数初始化五个 calc_force__core 子计时器；新增构造回归。 |
| `ASTRA-FD-20` | D4 | 周期边界闸门只判 `SEARCH_MODE_LONG`，漏掉 `LONG_SYMMETRY`/`LONG_SCATTER` | **已修改，待专项验证** | [tree_for_force_impl_exlet.hpp](FDPS/src/tree_for_force_impl_exlet.hpp)：周期边界的无 cutoff 拒绝检查同时覆盖 LONG、LONG_SYMMETRY、LONG_SCATTER，避免进入不支持 shift 的路径。 |
| `ASTRA-FD-21` | D3 | `Matrix2::operator*` 的 xy 分量写成 `xx + rhs.xy`，应为 `xx * rhs.xy` | **已修复（专项回归）** | [matrix2.hpp](FDPS/src/matrix2.hpp)：Matrix2 的 xy 乘积修正；新增四分量矩阵乘法数值回归。 |
| `ASTRA-FD-22` | D3 | `MortonKey::getKey` 把负坐标 `(U64)` 转换（UB），越界点静默折叠到 key 0/nmax；`half_len==0` 时 `0*inf=NaN` | **已修改，专项回归通过** | Morton 坐标在整数转换前检查有限性、范围和正 half_len；仅容忍边界舍入。half_len=1e-300、1、1e300 的端点测试通过；非默认周期域尚待专项验证。 |
| `ASTRA-FD-23` | D3 | 矩检查函数用 `tp` 位置直取 `ep[i]`/`sp[i]`，未按 `adr_ptcl_`/`ClearMSB` 映射 | **已修改，专项回归通过** | 矩检查按 tp 地址解码 EP/SP；混合列表质量检查通过。 |
| `ASTRA-FD-24` | D3 | LET-P2P 的 `pos_domain_allgather_[adr_sp_tmp]` 只有下界检查、无上界 | **待确认** | [tree_for_force_utils_moment.hpp](FDPS/src/tree_for_force_utils_moment.hpp)：adr_sp_tmp 仅检查 >=0，上界断言被注释；需证明合法 LET 排列能超界，本次只确认防御性缺口。 |
| `ASTRA-DG-01` | D4 | Tree_NB 计时偏移与 `time_profile_` 累加竞争；AR 计步与 `noteArRegion` 裸累加 | **已修复** | [tree_for_force_impl_force.hpp](FDPS/src/tree_for_force_impl_force.hpp)：原 FDPS Tree_NB 位置已用线程局部 offset 和 walk/dispatch reduction；原 SDAR AR 计数现串行，noteArRegion/inner_parallel 实现已不在当前源码。仅关闭这条列明的机制，不泛称所有 profile 无竞争。 |

### PeTar

| 编号 | 原等级 | 问题 | 当前状态 | 当前依据与条件（含源码） |
|---|---|---|---|---|
| `ASTRA-AR-02` | 条件 D1 | 递归 `modify_return` 取子分支最大值，局部 destroy 被提升为整子树 destroy | **已修复** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：modify_branch 分开记录；仅左右分支都返回 3 才传播 destroy，否则降为局部修改。修复位于 PeTar 的事件回调，SDAR 是消费端。 |
| `ASTRA-AR-03` | 条件 D1/D2 | 束缚↔非束缚切换时复用旧 `ecca`，`calcParticles` 只按 `semi` 符号选公式 | **已修复** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：postProcess 检测 semi 变号，按当前分离距离及径向速度恢复异常角，并重算轨道缓存；潮汐交接代码也已重写。原复用旧 ecca 的分支已消除；未重跑原五轮数值反例。 |
| `ASTRA-AR-05` | D3 | `AR_CHANGEOVER` 分支引用未声明的 `invr` | **已修复** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：AR_CHANGEOVER 分支已使用 inv_r；原未声明 invr 已消除。 |
| `ASTRA-AR-12` | D3 | `AR_SLOWDOWN_PERT_R4` 下 `calcPertFromForce` 返回 `G·(m_p m_pert/r⁴)`，与 `calcPertFromMR`/`calcPertFromBinary` 的 `r⁻⁴` 归一化不一致 | **已修改，待专项验证** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：R4 扰动指标分母改为 G²，单位缩放结果会改变。 |
| `ASTRA-AR-13` | D4 | `findCloseSoftPert` 在 `_n_tt==0` 时读 `_tt[-1].group_id` | **已修改，待专项验证** | [ar_perturber.hpp](PeTar_new_version/src/ar_perturber.hpp)：潮汐张量最近邻索引为 -1 时不再解引用。 |
| `ASTRA-H4-06` | D3 | `ForceSoft::clear()` 与 `FPSoft::clearForce()` 不清 `acorr`，4 阶修正跨步累积 | **已修改，待专项验证** | [soft_ptcl.hpp](PeTar_new_version/src/soft_ptcl.hpp)：力与粒子字段初始化，clear/clearForce 清除 KDKDK4 的 acorr。 |
| `ASTRA-H4-07` | D4 | `HermiteInformation::time_origin` 无初值、`clear()` 不置位 | **已修改，待专项验证** | [hermite_information.hpp](PeTar_new_version/src/hermite_information.hpp)：HermiteInformation 的 time_origin 构造和 clear 均初始化为 0。 |
| `ASTRA-BSE-01` | D2 | 紧致天体分支读取调用方未初始化的 `dtm` | **已修复** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：bse/bseEmp/mobse 的 trdot 紧致天体分支均为 age*10，不再用未初始化 dt；本次整套集成测试在到达 BSE 单测前中止，未把旧测试记录当新结果。 |
| `ASTRA-BSE-02` | D2 | 双星步长候选混用 IN 与 Myr 并重复除 `tscale` | **已修复** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：GR 与 trflow 候选各除一次 tscale 后比较，最终直接 return dt。 |
| `ASTRA-BSE-03` | D2 | 动力学 `merge` 直接调 `merge_`，未持 Fortran 共同锁 | **已修复** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：AR 中两个 bse_manager.merge 入口均在未命名 omp critical 内，与其它生产 BSE 调用使用同一锁。测试工具的并发问题另见 BSI-08。 |
| `ASTRA-BSE-04` | D2 | runtime 关闭恒星演化时 `modify_branch[2]` 未初始化即参与判定 | **已修复** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：modify_branch[2] 已显式初始化为 {0,0}。 |
| `ASTRA-BSE-05` | D2 | 异步 `evolveBinary` 复用 `_out.dm` 作阶段性缓存，输出错的 `dm` | **已修复** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：evolveBinary 入口保存 mass_initial1/2，正常与提前返回路径按入口质量计算 dm；错误返回路径不承诺有效结果。 |
| `ASTRA-BSE-06` | D3 | MIX 分支不写 `menv/renv`；COEL 分支写入的是合并前包层值，C++ 无条件拷贝 | **已修复（专项回归）** | [merge.f](PeTar_new_version/bse-interface/bse/merge.f) [merge.f](PeTar_new_version/bse-interface/bseEmp/merge.f) [merge.f](PeTar_new_version/bse-interface/mobse/merge.f) [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：经用户批准，同 BSE-07，在 CE/MIX 后刷新 menv/renv 并清空消失星的结构量；真实 MIX/CE 回归核对新包层参数。 |
| `ASTRA-BSE-07` | D2 | `merge_` 回写的 `r/mc/rc/ospin` 是合并前遗留值，`lum` 完全不更新；`r` 被 `mergerCheckRadius` 消费 | **已修复（专项回归）** | [merge.f](PeTar_new_version/bse-interface/bse/merge.f) [merge.f](PeTar_new_version/bse-interface/bseEmp/merge.f) [merge.f](PeTar_new_version/bse-interface/mobse/merge.f) [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：经用户批准，新增 merge_structure 返回刷新后的结构、光度、寿命；保留旧 merge ABI，自旋按 CE 返回角动量或 MIX 保留星原角动量更新。BH 拟合质量覆盖后同步 BH 核质量/半径。修正 MOBSE merge 中 hrdiag 参数签名。真实 MIX、CE 和 BH 合并结构回归通过；MIX 自旋角动量保持。 |
| `ASTRA-BSE-08` | D3 | `evolveBinary` 双 `kw==15` 早退不写事件终止符，调用方读未初始化 `record[9][*]` | **已修改，待专项验证** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：双 kw=15 的提前返回设置 event 结束索引为 0。 |
| `ASTRA-BH-01` | 条件 D1 | `BSEManager::merge` 把 Msun 质量与内部单位 G/dr/dv 混用 | **已修复** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：calculateRemnant 输入质量改为 getMass(star)，与内部 G、dr、dv 一致；不再直接传太阳质量 mt。 |
| `ASTRA-BH-02` | 条件 D1 | remnant 方位角用无向 `acos`，丢绕角动量方向 | **已修复** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：投影方位角改 atan2，Python 端保持初始 cyclesign；本次桥接正负分支及 precession_sign 测试通过。 |
| `ASTRA-BH-03` | D3 | Python 调用返回 NULL 后直接 `PyArray_Check` | **已修复** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：validArray/arrayValue 先检查空指针与类型，再检查形状/数据；四个 remnant 异常注入均安全返回失败，本次桥接回归通过。 |
| `ASTRA-BH-04` | D4 | 每次 BH 并合泄漏 kwargs float、空 tuple 与一个 array 引用 | **已修复** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：kwargs 浮点、调用 tuple、结果和临时 array 引用有对应释放；本次 2000 次桥接调用保留内存增长 10805 B，小于测试阈值 65536 B。 |
| `ASTRA-BH-05` | D3 | `e==0` 时 `inspiral_precav` 的 `e` 返回 `(1,1)`，C++ 按 `[0,1]` 越界读 | **已修复（专项回归）** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：经逐项批准，仅圆轨道允许 e 的 (1,1) 零常量输出；其他轨道字段读取末端样本。注入合法圆轨道结果后，remnant 拟合收到末端 deltaphi=0.375，而不是遭遇时角度。 |
| `ASTRA-BH-06` | D3 | `remnantkick` 返回非数组时静默零 kick 且不 abort；无形状/有限性检查 | **已修复** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：remnantkick 使用 validArray(result,4)，错误类型/长度返回 false，不再静默提交零 kick；本次错误返回测试通过。 |
| `ASTRA-BH-07` | D4 | `PyArray_SimpleNewFromData` 包装栈数组并额外 `INCREF`，每次 inspiral 泄漏且留下悬垂指针 | **已修复** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：栈数组包装仍是同步调用，但多余 INCREF 已去除，kwargs/result/a_array 按序释放；原泄漏引起的遗留悬垂对象机制已消除，不承诺第三方恶意保留输入的安全性。 |
| `ASTRA-BH-08` | D3 | 束缚轨道 `e = sqrt(1+2EL²/(μG²m₁²m₂²))` 未裁剪，`e≲1e-8` 时舍入成负数 → `e=NaN` | **已修复（专项回归）** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：经逐项批准，e² 仅在 32 epsilon 的消去舍入容差内裁为 0；明显负值或非有限项返回失败，不误走直接并合。边界和明显错误输入回归通过。 |
| `ASTRA-BH-09` | 条件 D2 | 安装树不含 `precession` 包；`sys.path` 用相对 CWD 的 `./bhmerger/precession-master` | **已修复（专项回归）** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp) [Makefile.in](PeTar_new_version/Makefile.in)：经逐项批准，构建嵌入绝对源码包/安装包路径，支持 PETAR_PRECESSION_PATH 覆盖；新增安装包目标。从 /tmp 启动实际导入源码包和暂存包通过，未安装到用户 bin。 |
| `ASTRA-BH-10` | D3 | 退化/近对齐自旋：`χ=0` 的 0/0、`θ=0/π` 的 `arccos(1±ulp)`、`deltachiroots` 近重根被 `wraproots` 置 NaN、随后 `np.random.uniform(0,τ=∞)` | **部分修复** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py) [eccentricity.py](PeTar_new_version/bhmerger/precession-master/precession/eccentricity.py)：经批准只处理角度端点舍入：Python arccos 共享有限性/32 epsilon 边界检查，明显越界报错；原来的无条件 clip 同步收紧。端点和相位符号回归通过。零自旋及未证明为重根的退化问题仍未关闭。 |
| `ASTRA-BH-11` | D3 潜伏 | `initializePython` 在宿主已初始化 Python 时早退，跳过 `_import_array()` 与 `sys.path` 设置 | **已修复（专项回归）** | [inspiralZJ.cpp](PeTar_new_version/bhmerger/inspiralZJ.cpp)：经逐项批准，解释器存在与 NumPy API/路径初始化分开记录，在 GIL 下初始化桥接。宿主预初始化与独立解释器两种测试通过；finalizePython 不关闭宿主解释器。宿主自行反复销毁/重建 NumPy 解释器不在支持范围内。 |
| `ASTRA-BH-12` | D4 | `odeint` 容差守卫拼写为 `'rol'`/`'aol'`，用户 `rtol/atol` 被静默覆盖 | **仍存在** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：precession.py 容差守卫仍拼成 rol/aol，传 rtol/atol 可能被默认值覆盖；仅使用这些接口/显式容差时。 |
| `ASTRA-BH-13` | D4 现在 / 前瞻 D2 | `rhs_precav` 返回 `float(1-element ndarray)`，NumPy≥2 弃用 | **已修改，当前环境回归通过** | 单元素数组改用 np.asarray(Ssav).item() 提取标量；当前进动回归通过，未在 NumPy 2.x 单独验证。 |
| `ASTRA-BH-14` | D4 | `np.in1d`（NumPy≥2 弃用） | **已修改，当前环境回归通过** | 两处 np.in1d 改为 np.isin，输出筛选语义保持；当前进动回归通过，未在 NumPy 2.x 单独验证。 |
| `ASTRA-BH-15` | D3 | `precomputedroots[:-1]` 把 3 行根解包给 3 个名字 | **仍存在** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：仍用 deltachimin,deltachiplus,_ = precomputedroots[:-1]；三行输入被切成两行，指定快路径失败。 |
| `ASTRA-BH-16` | D3 | `return np.stack[(which, kappa, chieff)]` 索引 ufunc 而非调用 | **仍存在** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：widenutation_condition 仍 return np.stack[(which,kappa,chieff)]，将函数当下标对象。 |
| `ASTRA-BH-17` | D3 | `eval_costheta12` 使用未定义 `r`（签名无此参数） | **仍存在** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：eval_costheta12 签名无 r，守恒量分支仍使用 r；角度分支不受影响。 |
| `ASTRA-BH-18` | D3 | `inspiral_hybrid` 把 `deltachi`/`r` 传给不接受它们的 `inspiral_precav`；残留 `print("uswitch, eswitch",…)` | **待确认** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：hybrid 的可疑关键字传递和 eccentricity 中调试输出仍在；两套模块函数签名不同，原报告的完整吸积数值失败未重现，不把字面参数相似直接当已证实根因。 |
| `ASTRA-BH-19` | D4 | `remnantspindirection` 文档说返回 `thetaL` 实返 `Ĵ`；`remnantkick` 文档默认 `kms=False/full_output=False` 实为 True；`deltachiroots`/`kappalimits`/`deltachilimits` 的 Returns 顺序与代码相反 | **已修改** | 修正五处接口文档：上下界/根返回顺序、hatJ 向量及 kick 默认开关；未修改这些接口的返回值和默认值。 |
| `ASTRA-BH-20` | D3 | `rhs_precav` 的 `q==1` 分支丢弃第二项；近 q=1 分支本身病态 | **待确认** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：q==1 的分支确实缺交叉项，C++ 仍用 1-1e-10 避开；等质量奇异极限的正确物理表达需解析推导/数值复核，不能仅以近 q 的不稳定值证明。 |
| `ASTRA-BH-21` | D4 | `eval_S`/`eval_kappa` 的 `TypeError` 未 `raise`（第二分支写成 `if`）；`eval_q` 断言 `q<1` | **已修改，专项回归通过** | 互斥分支使用 elif，非法参数显式抛 TypeError，eval_q 接受 0<q<=1；平行/反平行自旋和等质量接口测试通过，不代表所有演化式的等质量奇点已解决。 |
| `ASTRA-BH-22` | D3 | `deltachilimits_definition` 用 `\|χ₁−qχ₂\|/(1+q)` 作 Δχ 上下限，与 `eval_deltachi` 的实际跨度 `±(χ₁+qχ₂)/(1+q)` 不符 | **已修改，专项回归通过** | 几何边界改为 ±(chi1+q*chi2)/(1+q)，极点及 101×101 角度网格包络检查通过。 |
| `ASTRA-BH-23` | D3 | `deltachitildeav2` 在 `m→0` 用 `tol=1e-7` 硬夹，但闭式在该区间灾难性相消 | **已修改，专项回归通过** | m<=0.01 使用五阶级数，m=0 返回 3/8；保留上端 regularization。100 位高精度参考覆盖 0、1e-30 至 0.99 及切换点，最大绝对误差 4.70e-14。 |
| `ASTRA-BH-24` | D3 | `kapparesonances` 的 up-down/down-up 递归分支不可达；数组输入被 `map` 截断到最短维 | **待确认** | [precession.py](PeTar_new_version/bhmerger/precession-master/precession/precession.py)：分支顺序和 map 广播截断代码未改；数组输入维度契约及端点递归失败率需专门复现，原 3458/20000 不是本次结果。 |
| `ASTRA-BH-25` | 条件 D2 | 并合事件日志写积分帧 `pos/vel`，快照写 `pos+pcm.pos`（仅 `ext_mode!=off`）；`d8fc857` 又删掉 `center.pos/vel` 六列，文件内无列可反推 | **仍存在** | [bhmerger_output.hpp](PeTar_new_version/src/bhmerger_output.hpp)：capture 直接保存粒子局部 pos/vel，未保存 CM offset；非零参考系偏移的跨文件比对仍缺信息，默认 ext_mode=off 不触发该条件。 |
| `ASTRA-BH-26` | D3 | 并合日志只记检测时刻 `time_now`，未抓 `time_merger/merger_semi/merger_ecc`（随后被清零） | **仍存在** | [bhmerger_output.hpp](PeTar_new_version/src/bhmerger_output.hpp)：并合 ParticleState 仍不含预测截止时刻和轨道检查点；新增 change 日志记录部分变化，不能替代每个 merger 的预测状态记录。 |
| `ASTRA-BH-27` | D3 | `write()` 忽略 `writeLocked` 的返回值；写循环中途失败会留下**半行**，事件被静默丢弃 | **已修改，专项回归通过** | 写失败明确报告文件名；持锁回滚普通文件的半行，回滚失败另报错；检查 close 错误。注入短写后 EIO 验证原文件保持完整，随后正常写入成功。硬件/文件系统拒绝回滚时不能保证无半行。 |
| `ASTRA-BH-28` | D3 | 表头只在 `!_append && _initialize_file` 时写，且 `_initialize_file` 仅 rank0 为真；rank0 表头写失败只把它自己的 `enabled_` 置 false | **已修改，串行专项回归通过** | rank0 在 append 新/空文件时持锁补表头，非空不重复；初始化状态广播给全部 rank，失败集体中止。新文件/重复 append/失败路径已测试，MPI 多 rank 故障注入未运行。 |
| `ASTRA-BH-29` | D3 | 表头/数据无单位与参考系标注；time/semi/pos 用代码单位，姊妹 `.bse` 日志同一事件乘 `tscale·year_to_day` 与 `rscale` | **仍存在** | [bhmerger_output.hpp](PeTar_new_version/src/bhmerger_output.hpp)：bhmerger 列标题仍没有单位/参考系说明；确有解释歧义，但原表“天与 Myr 相差约365”的数量级表述不能照用，应结合 tscale 和 year_to_day。 |
| `ASTRA-BH-30` | D4 | `capture()` 不记录 `binary_state`（pair id/中断态） | **仍存在** | [bhmerger_output.hpp](PeTar_new_version/src/bhmerger_output.hpp)：capture 仍不保存 binary_state/pair 状态；粒子 id 可追溯但不等价于保存中断态。 |
| `ASTRA-TR-01` | 条件 D1 | F32 邻居计数不具平移稳健性，真实连接被算成 `n_ngb==1` 而误归 single | **已修复** | [phantomquad_for_p3t_x86.hpp](PeTar_new_version/src/phantomquad_for_p3t_x86.hpp)：软力 wrapper 用原始 double 坐标复核邻居数，SIMD 力计算仍保留；本次整体平移邻居检查通过。不能扩展为 F32 力的所有精度问题都已解决。 |
| `ASTRA-TR-05` | D2 | quad + `AVX_PRELOAD` 下 `spjbuf[j+1]` 越过 `NJMAX` | **已修复** | [phantomquad_for_p3t_x86.hpp](PeTar_new_version/src/phantomquad_for_p3t_x86.hpp)：AVX_PRELOAD 末轮 j+1 已夹到 j，原 NJMAX 越界表达式已改；本次未重跑 NJMAX 极限力核压力测试。 |
| `ASTRA-TR-06` | D3 | `INTEGRATED_CUTOFF_FUNCTION` 下 `q_` 从不赋值，`readBinary` 允许半读 | **已修复** | [changeover.hpp](PeTar_new_version/src/changeover.hpp)：q_ 构造/两种 setR/读回均赋值，readBinary 要求恰好两项；本次势有限性与二进制往返通过。 |
| `ASTRA-TR-07` | D3 | 稳定性 `acos` 只夹上界，`ratio<-1` 得 NaN 后被判"稳定" | **已修复** | [stability.hpp](PeTar_new_version/src/stability.hpp)：PeTar stability 和 SDAR binary_tree 两边都夹 [-1,1] 并处理非有限比值；属于跨模块修复。 |
| `ASTRA-TR-08` | D3 | CM 预测位移写成 `dx²+2·dy·dz`，漏 `dz²` | **已修改，专项回归通过** | CM 预测位移范数修正为 dx²+dy²+dz²；非默认 SMOOTH_CM_USING_VEL_PRED 完整轨迹尚未运行。 |
| `ASTRA-TR-09` | 非缺陷 | wrapper 只装载 `n_ip_local` 个目标却给 kernel 传 `n_ip` | **非缺陷** | [soft_force.hpp](PeTar_new_version/src/soft_force.hpp)：wrapper 仍可向 kernel 传 n_ip，但有效结果只写回 n_ip_local；原登记册本就列为额外工作量，不计入未修复正确性缺陷。 |
| `ASTRA-TR-11` | 条件 D1 | 簇搜索两体判据把 `G` 写死 `1.0`；`-u 1` 生产配置下运行时 `G=4.4983e-3` → `peri` 小 222–250 倍，判据放宽 | **已修改，待专项验证** | [petar.hpp](PeTar_new_version/src/petar.hpp)：邻居速度/近心判据传入运行时 gravitational_constant；G 非 1 时会改变分组，尚未做星团轨迹数值对照。 |
| `ASTRA-TR-12` | D2 | `id_to_adr_pcluster_[...]` 用 `operator[]`，缺键静默插入 0，把边接到 `ptcl_cluster_[0][0]` | **待确认** | [cluster_list.hpp](PeTar_new_version/src/cluster_list.hpp)：operator[] 缺键自动插零仍在，但本次未证明合法建表/通信流程会缺键；记录为防御性缺口，不沿用“已造成错误连接”的强结论。 |
| `ASTRA-TR-13` | D3 | 判据对 NaN 无防护（`mass=0`/`r=0`/`e²<0`）；NaN 使 `peri<r_crit` 恒 false，静默丢边 | **仍存在** | [cluster_list.hpp](PeTar_new_version/src/cluster_list.hpp)：particleToSemiEcc 后未做有限性检查就比较 peri；退化质量/重合位置等输入可能丢边，合法生产输入触发率未测。 |
| `ASTRA-TR-14` | D4 | 覆写 `sys[i].n_ngb` 为「判据通过数、不含自身」，与 `:454` 单粒子分支保留的 1 语义漂移；只影响输出列 `N_b` | **已修改，待专项验证** | [cluster_list.hpp](PeTar_new_version/src/cluster_list.hpp)：写回 n_ngb 包含自身，与力计算计数约定一致。 |
| `ASTRA-TR-15` | D3 | `FPSoft(){}` 不初始化 `acc/pot/pos/vel`；转换构造器不补 `acorr` | **已修改，待专项验证** | [soft_ptcl.hpp](PeTar_new_version/src/soft_ptcl.hpp)：FPSoft 的邻居计数、rank、地址等瞬态字段显式初始化。 |
| `ASTRA-TR-16` | D4 | `group_data` 是 `artificial` 与 `cm` 的 16 B 完全重叠 union；`EPISoft::copyFromFP` 无条件按 `.artificial` 判 `type`，而建树时 mode 可能仍是 `cm` | **仍存在** | [soft_ptcl.hpp](PeTar_new_version/src/soft_ptcl.hpp)：group_data union 及 EPISoft 按 artificial 解释的路径未改；当前正常力调用顺序屏蔽后果，属于潜伏状态解释风险。 |
| `ASTRA-TR-17` | D3 | `const PS::F64vec da = _pi.acc - _pi.acc;` 恒为 0；KDKDK4 的 `acorr` 丢掉全部 `da` 项 | **已修改，待专项验证** | [hard.hpp](PeTar_new_version/src/hard.hpp)：相对加速度改为 ai-aj，消除 ai-ai 恒零。 |
| `ASTRA-TR-18` | D3 | `setDKDMode()` 把 `mode_` 置 1 后，`getDtStartContinue()` 返回 `coff_one_step_[0][1]`（drift 系数 `ds`），而调用方按 kick 使用；`next_is_kick_flag_` 也被翻反 | **已修改，专项回归通过** | continue 首尾步使用对称分解首列系数，DKD 时对应漂移；保持中间 kick/drift 状态切换。KDK/DKD 和普通/KDKDK 系数两轮循环均验证两类步长总和。 |
| `ASTRA-TR-19` | D3 潜伏 | `CalcForceEpSpQuadSimd` 用 j 全局下标 `i` 而不是槽位 `i_tmp` 调 `set_spj_one` | **已修改，待专项验证** | [soft_force.hpp](PeTar_new_version/src/soft_force.hpp)：quad 分块预装使用块内 i_tmp 索引。 |
| `ASTRA-TR-20` | D4 | `KickDriftStep` one-step API 无调用者，`%=coff_one_step_.size()` 在 `setStep` 前除零；`coff_[0]` 空表下标 | **已修改，专项回归通过** | 读步长或推进前显式检查 setStep 已生成表，避免空表下标/取模；未初始化调用已验证明确中止。 |
| `ASTRA-IO-01` | 条件 D1 | ASCII 检查点不保存 `time_merger/merger_semi/merger_ecc`，重启按 `now+lifetime` 重算 | **已修复（专项回归）** | [soft_ptcl.hpp](PeTar_new_version/src/soft_ptcl.hpp)：经用户批准，ASCII #BHMERGER_V2 保存 deadline、a/e、检查点时间及累计动力学量；兼容 V1/无后缀，旧文件缺失历史仍无法补回。新增重启后首次动力学检查保持 deadline 的回归。 |
| `ASTRA-IO-03` | 契约 | 常规重启未串接 Fortran `/RAND3/`+`/VALUE3/` 与 NumPy RNG 状态；`inspiralZJ.cpp` 连 seed 都不设 | **契约缺口** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：常规快照仍不保存 Fortran RAND3/VALUE3 及 NumPy 完整 RNG 状态；去掉相位随机分支、零自旋扰动改确定性并不等价于完整重启可复现。 |
| `ASTRA-IO-04` | D3 | `format_transfer` 的 `FileHeaderNoOffset` 用 `sizeof(FileHeader)` 读写 | **已修改，待专项验证** | [format_transfer.cxx](PeTar_new_version/src/format_transfer.cxx)：FileHeaderNoOffset 按自身 sizeof 读写。 |
| `ASTRA-IO-06` | D3 | 二进制快照 `fwrite(&acc, 6)` 依赖字段布局；`KDKDK_4TH` 在 `acc` 后插入 `acorr[3]` → 写出 acc+acorr，`pot_tot/pot_soft/n_ngb` 从未落盘 | **已修改，待专项验证** | [soft_ptcl.hpp](PeTar_new_version/src/soft_ptcl.hpp)：FPSoft 二进制尾部逐字段写入 acc、势和 S64 邻居数，避免 acorr 插入改变布局；旧的错误 KDKDK4 快照不能恢复丢失字段。 |
| `ASTRA-IO-07` | D3 | 快照头与粒子均为裸 `sizeof` 内存转储，**无 magic/version/构建签名**；`FileHeader` 在 `RECORD_CM_IN_HEADER` 下多 6 个 double，粒子结构随 `BHMERGER`/`BSE` 宏改变 | **仍存在** | [io.hpp](PeTar_new_version/src/io.hpp)：常规快照头/ParticleBase 仍原始内存写读，无全局版本/宏签名；ARI2/SYS2 仅修特定 SDAR 格式，不覆盖此项。 |
| `ASTRA-IO-08` | D4 | 重启后 `initialStep` 无条件 `file_header.nfile--; output()`，使初始快照写回**重启输入文件本身**（data.N 被覆盖） | **已修改，专项回归通过** | 初始输出检测目标与输入是否为同一文件（包括硬链接），跨 rank 广播保护标志；保留状态输出及后续编号，新前缀仍可写初始快照。 |
| `ASTRA-EN-01` | D2 | 孤立粒子演化阶段对四个共享能量字段裸 `+=` | **已修复** | [hard.hpp](PeTar_new_version/src/hard.hpp)：孤立演化以 de_kin_total reduction 收集，循环外一次写四个累计字段。 |
| `ASTRA-EN-02` | D2 | 逃逸删除在 OpenMP 循环内对四个共享能量字段裸 `-=` | **已修复** | [petar.hpp](PeTar_new_version/src/petar.hpp)：当前 removeParticles 以完整删除集合求 energy_change，四个账本字段在并行区外更新；原裸 -= 竞争已消除。 |
| `ASTRA-EN-03` | D3 | 删除动能只减 `0.5 m v_rel²`，未加回 `pcm.vel` | **已修复** | [petar.hpp](PeTar_new_version/src/petar.hpp)：RECORD_CM_IN_HEADER 下 energy_vel_offset=stat.pcm.vel，energyChange 用 vel+offset 算动能；代码已修，外场构建本次未专项验证。 |
| `ASTRA-EN-04` | D3 | `writeBinary` 从 `&ekin` 起写整个对象 | **已修复** | [energy.hpp](PeTar_new_version/src/energy.hpp)：EnergyAndMomentum::writeBinary 从 this 开始；本次非首字段与角动量往返回归通过。 |
| `ASTRA-EN-05` | D2 | 同一步删除多个逃逸粒子时逐粒子扣除 `m_i Φ_i`，把被删集合内部互势扣两次 | **已修复** | [petar.hpp](PeTar_new_version/src/petar.hpp)：energyChange 在 -sum(K+m*phi) 后对被删集合 i>j 互势加回一次；当前工作区 B45 补丁明确消除重复扣除。 |
| `ASTRA-EN-06` | D3 | `finishInterruptDrift` 用局部 `mass_modify_list` 遮蔽成员，中断收尾的 `dm` 不进能量账本也不清零 | **已修改，待专项验证** | [petar.hpp](PeTar_new_version/src/petar.hpp)：删除遮蔽成员的局部 mass_modify_list；多 rank 中断场景尚未专项运行。 |
| `ASTRA-BL-01` | D3 | `--id-offset` 用 32-bit `atoi` 解析 `PS::S64` | **已修复** | [petar.hpp](PeTar_new_version/src/petar.hpp)：id-offset 已 strtoll，并检查 ERANGE、尾字符和范围，不再 atoi 截断。 |
| `ASTRA-BL-03` | D3 | `--disable-print-info` 早于 `-p` 时整段参数文件读取被跳过 | **已修改，待专项验证** | [petar.hpp](PeTar_new_version/src/petar.hpp)：参数文件读取只受 rank==0 控制，print_flag 仅控制打印。 |
| `ASTRA-BL-04` | D3 | `checkParams` 的 `data_format` 条件用逻辑或连接、恒真；NDEBUG 下 `-i 4` 静默不写快照 | **已修改，待专项验证** | [petar.hpp](PeTar_new_version/src/petar.hpp)：数据格式使用运行时范围检查，拒绝 0..3 以外取值。 |
| `ASTRA-BL-05` | D3 | `ForceSoft::grav_const` 缺省 `1.0` 与 par 缺省同值；AU 数据若既不写 G 也不带 `-u 1`，引力静默弱 39.5× | **待确认** | [static_variables.hpp](PeTar_new_version/src/static_variables.hpp)：G=1 是无量纲/N-body 的合法默认值；把未声明单位的 AU 数据送入该默认模式，不能单凭 39.5 倍差异证明程序缺陷，需明确输入单位契约。 |
| `ASTRA-BL-06` | D3 | `HERMITE_PERT_FORCE` 分支调用的 `calcAccJerkPerturber` 在 `HermiteInteraction` 中不存在 | **仍存在** | [hermite_interaction.hpp](PeTar_new_version/src/hermite_interaction.hpp)：SDAR 的可选泛型接口调用 calcAccJerkPerturber，但 PeTar HermiteInteraction 未实现；主责归 PeTar 适配层，SDAR 联动，仅 HERMITE_PERT_FORCE。 |
| `ASTRA-BL-07` | D3 | `static HardDumpList hard_dump;` 写在头文件里（内部链接），每个包含它的 TU 各有一份；唯一 `initial()` 在 `petar.hpp:3084` | **已修改，待专项验证** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：hard dump 初始化先 clear，线程数及索引验证后再访问；全局对象采用 C++17 inline。 |
| `ASTRA-BL-08` | D3 | `checkParams` 的软化上界 assert 写成 `eps <= ratio_r_cut`（无量纲），意图是 `eps <= r_in = ratio_r_cut*r_out` | **已修改，待专项验证** | 软化长度先检查有限非负；自动/显式 changeover 尺度确定后与实际 base r_in 比较，不再与无量纲 ratio_r_cut 比较。 |
| `ASTRA-BL-09` | D4 | `checkParams` 对 `data_format` 的 assert 写成 `>=0 \|\| <=3`（恒真）；命令行路径 `case 'i'` 的 assert 才是正确的 `&&` | **已修改，待专项验证** | [petar.hpp](PeTar_new_version/src/petar.hpp)：与 BL-04 共用一个范围校验修复，不重复计为独立补丁。 |
| `ASTRA-BL-10` | D4 | `HardDumpList`：`dumpThread` 直接索引不校验；`initial()` 重复调用泄漏；`readOneCluster` 先按未校验 `n_group` resize 再查 `<0` | **已修改，待专项验证** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：dump 数量在 resize/read 前检查为非负。 |
| `ASTRA-BL-11` | D3 潜伏 | `static_variables.hpp` 头内定义 `Ptcl::`/`EPISoft::`/`ForceSoft::` 静态成员；`r_out=0`、`search_factor=0` 毒默认 | **已修改，待专项验证** | [static_variables.hpp](PeTar_new_version/src/static_variables.hpp)：静态变量定义改为 C++17 inline，消除多翻译单元重复定义。 |
| `ASTRA-BL-12` | D3 | `regularTimeStep` 对 `_dt=0`/`+inf` 死循环；`dt_soft` 的零/无穷回退无保护 | **已修改，待专项验证** | [petar.hpp](PeTar_new_version/src/petar.hpp)：regularTimeStep 拒绝非有限或非正输入；单粒子自动步长经过同一校验。 |
| `ASTRA-BL-13` | D4 | `get_version.hpp` 是 `@VERSION@` 模板产物但 `Makefile` 无生成规则；banner 版本恒为写死的 `1268_298`；`printFeatures` 不打印 `HARD_DUMP`/`AR_DEBUG_DUMP` | **部分修复／澄清** | 补充 HARD_DUMP/AR_DEBUG_DUMP 功能提示。复核发现 configure.ac 已通过 AC_CONFIG_FILES 生成 get_version.hpp，版本来自 PeTar/SDAR VERSION；不能把缺 Makefile 规则单独判为未生成。已生成头不会随 VERSION 自动刷新这一构建依赖问题仍保留。 |
| `ASTRA-DG-02` | D4 | `PtclHard::print` 的 `id_cluster`/`adr_org` 写到 `std::cerr` | **已修复** | [hard_ptcl.hpp](PeTar_new_version/src/hard_ptcl.hpp)：PtclHard::print 的 id_cluster/adr_org 改写 _fout；本次集成测试前段通过。 |
| `ASTRA-DG-03` | D4 | 并合路径无条件输出 `DEBUG Dynamic_merge invoked` 到 stderr 与 `fout_bse` | **已修复** | [ar_interaction.hpp](PeTar_new_version/src/ar_interaction.hpp)：Dynamic_merge 调试行现受 AR_DEBUG_PRINT 保护，默认构建不输出该行。 |
| `ASTRA-DG-04` | D4 | `-w 2` 的 `.status` 表头在读数据前写出，粒子列标题全缺 | **已修改，专项回归通过** | status 表头延后至首次 output、粒子数量已知时写出；空 append 文件也补标题，非空文件不重复。 |
| `ASTRA-DG-05` | D4 | `_n_ptcl>1000` 的“Large cluster”无条件写 stderr | **已修改，待专项验证** | [hard.hpp](PeTar_new_version/src/hard.hpp)：Large cluster 提示受 HARD_DEBUG_PRINT 控制，不再无条件输出。 |
| `ASTRA-DG-06` | D4 | `checkNeighborExist` 把平方距离 `r_min_sq` 与半径 `r_out` 直接比较 | **已修改，待专项验证** | [hard.hpp](PeTar_new_version/src/hard.hpp)：邻居诊断将平方距离与平方外半径比较。 |
| `ASTRA-DG-07` | D3 | `dump_number++` 在并行路径非原子；`hard_large_energy` 的 DATADUMP 位于并行 cluster 循环 | **已修改，待专项验证** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：dump_number 使用原子递增。 |
| `ASTRA-DG-08` | D4 | `calcEnergyPertOneGroup` 不使用传入的 `_perturber`，实取 `_group.perturber.soft_pert`；NULL 时静默 0 | **待确认** | [hermite_interaction.hpp](PeTar_new_version/src/hermite_interaction.hpp)：确实不使用 _perturber，而用 group.perturber.soft_pert；当前力/势取同一组对象，没有已证实不一致，更像接口冗余，不能仅凭未使用参数定 bug。 |
| `ASTRA-DG-09` | D4 | `Tprofile::end()` 连续两次 `PS::GetWtime()`，`tbar` 与 `time` 在不同时刻采样 | **已修改，待专项验证** | [profile.hpp](PeTar_new_version/src/profile.hpp)：Tprofile::end 只取一次时间用于 time/tbar。 |
| `ASTRA-DG-10` | D4 | `NumCounter* iptr = (NumCounter*)this+i` 跨成员指针算术（技术上 UB） | **已修复（专项回归）** | [profile.hpp](PeTar_new_version/src/profile.hpp)：SysCounts/SysProfile 用显式成员指针表遍历，消除跨独立成员的指针算术；保持列顺序及 DARKMATTER 条件列。计数 clear 回归通过。 |
| `ASTRA-DG-11` | D4 | dump 只存 4 个 `Ptcl::` 静态量，不存 `EPISoft::eps/r_out` 与 `ForceSoft::grav_const` | **仍存在** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：dump 保存 Ptcl 四静态量，未保存 eps/r_out/G；离开配套 par 文件的独立重放契约不完整。 |
| `ASTRA-DG-12` | D3 | `readOneCluster` 先 `resizeNoInitialize(n_group)` 再检查 `n_group<0` | **已修改，待专项验证** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：与 BL-10 共用负数量前置检查。 |
| `ASTRA-DG-13` | D4 | `readOneCluster` 覆写进程全局的 `Ptcl::search_factor/r_search_min/mean_mass_inv/r_group_crit_ratio` | **待确认** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：读 dump 确实恢复全局 Ptcl 参数；重放本来需要这些参数。若要求读取后恢复调用前全局状态才是缺口，当前单次工具流程不能直接定 bug。 |
| `ASTRA-DG-14` | D3 | `ASSERT` 宏展开为裸 `if(!(expr)){...}`，无 `do{}while(0)`、无尾分号 | **已修改，待专项验证** | [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：ASSERT 宏改为 do/while(0)。 |
| `ASTRA-DG-15` | D4 | profile Min/Max 行取 `time−tbar`，普通行取 `time`；`getMax/getMin` 只回填 `time`，`tbar` 留 0 | **已修复（专项回归）** | [profile.hpp](PeTar_new_version/src/profile.hpp)：Min/Max 与普通行均统计 time，barrier 列独立归约 tbar；单 rank 对照通过。多 rank 的 time/tbar 极值可以来自不同 rank，不代表同一 rank 的联合样本。 |
| `ASTRA-DG-16` | D4 | 错误信息中的读取数与实际不一致 | **已修改，待专项验证** | [soft_ptcl.hpp](PeTar_new_version/src/soft_ptcl.hpp) [hard_assert.hpp](PeTar_new_version/src/hard_assert.hpp)：粒子和 hard dump 的读取数错误提示分别修正为 5 和 4。 |
| `ASTRA-GEN-01` | D3 | Kepler 盘三个分布函数同用 `init_genrand(getRank())`，`a` 与 `e/i` 取自同一 MT 流的相邻位置；`seed` 只改角度/相位 | **已修改，专项回归通过** | Kepler 抽样助手共享调用方按 rank/seed 初始化的单一随机流，去掉各自同种子的静态流；相同 seed 重复调用可重现，序列相对旧版会改变。 |
| `ASTRA-GEN-02` | D4 | `DomainDecision` 全树无调用点（实际走 `dinfo.decomposeDomain()`） | **非缺陷** | [domain.hpp](PeTar_new_version/src/domain.hpp)：DomainDecision 未被使用是代码覆盖/死代码事实，本身不是运行 bug；按现状从未修复缺陷数中排除。 |
| `ASTRA-GEN-03` | D3 | `assert(n_proc%4==0)` + 同名 `if(n_proc%4==0){…}` 包住整个函数体 | **已修改，专项回归通过** | DomainDecision 对非 4 倍进程数显式报错并中止，不再依赖 NDEBUG 下被删除的 assert；当前生产路径无调用。 |
| `ASTRA-GEN-04` | D3 | Plummer 归一化 `coef=1/sqrt(r_scale)` 隐含 G=M=1，缺 `sqrt(G·mass_glb)` | **已修改，专项回归通过** | Plummer 增加默认 G=1 的显式参数，尺度 r_scale 包含 G，速度包含 sqrt(G*M/r_scale)；PeTar 主入口传运行时 G。保留现有截断半径处方。 |
| `ASTRA-GEN-05` | D3 | Kepler 盘生成器默认 `gravitational_constant=1.0`（AU/Msun/5e6s 体系）而 peTar 用 pc/Myr；且不生成中心恒星 | **待确认** | [particle_distribution_generator.hpp](PeTar_new_version/src/particle_distribution_generator.hpp)：Kepler 生成器的 G 默认1及不放中央星是事实；是否要求调用方自行加入中心星/换单位需要接口契约，当前入口已注释，不能直接升级物理错误。 |
| `ASTRA-BSI-01` | D2 | RLOF/接触分支对流阻尼因子写成 `ttid/(2tc)**2`（量纲 1/时间），detached 分支在下一行改正、接触分支未改 | **已修改，待专项验证** | [evolv2.f](PeTar_new_version/bse-interface/bse/evolv2.f)：经用户批准，接触分支潮汐因子改为 (ttid/(2tc))²；编译通过，尚未完成长期双星潮汐轨迹对照。 |
| `ASTRA-BSI-02` | D2 | `gntage` 对 `mc<mcagbf(Mhef)` 强制 `mmin=Mhef`，非单调 `mcheif` 的低质量根被排除，且无 `m0<=mt` 校验 | **待确认** | [gntage.f](PeTar_new_version/bse-interface/bse/gntage.f)：gntage 的 mmin=Mhef 限制仍在，但 m0 是演化等效初始质量，经历失质量时可大于当前 mt；原 m0>mt 不能单独作为错误证据，需证明漏根/年龄不符合模型。 |
| `ASTRA-BSI-03` | D2 | CE 非并合分支在 `kick()` 后用 kick 前 `EORBF/ECIRC` 重算 `ECC`；MOBSE 有 `SNexplosion` 守卫，`bse`/`bseEmp` 没有；源码注释自认 “should not be done in case of CE SN ... fix” | **已修改，待专项验证** | [comenv.f](PeTar_new_version/bse-interface/bse/comenv.f) [comenv.f](PeTar_new_version/bse-interface/bseEmp/comenv.f)：经用户批准，bse/bseEmp 的三个 CE kick 调用设置 SNEXPLOSION，仅无超新星时用旧能量更新 ECC；MOBSE 原有守卫保留。编译及控制流核对通过，尚未构造真实 CE 内超新星存活样例。 |
| `ASTRA-BSI-04` | D3 | 默认 `bse/instar.f` 缺 `ktype(11,11)=12`；`bseEmp`/`mobse` 都有 | **待确认** | [instar.f](PeTar_new_version/bse-interface/bse/instar.f)：bse 缺 ktype(11,11)=12，而 bseEmp 有 Tanikawa DD model 注释；变体模型差异本身不是 bug，需确定默认 BSE 应采用哪种双白矮星产物模型。 |
| `ASTRA-BSI-05` | D3 | BSE 解析器复用 PeTar `-p`，无条件打开 `<prefix>.bse`，缺失即 `abort()`；`stellar-evolution=0` 时也执行 | **仍存在** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：BSE 解析 -p 后仍读取追加 .bse 的文件，未受 stellar-evolution=0 控制；但并非绝对无条件，代码有 print_flag 守卫，应修正原描述。 |
| `ASTRA-BSI-06` | D3 | `BSEEMP` 帮助允许 `--bse-metallicity 0`，但 `checkParams()` 断言 `z>0`，`zcnsts` 对 `log10(0)` 无保护 | **已修改，待专项验证** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：BSEEMP 金属丰度说明与校验统一为有限正数。 |
| `ASTRA-BSI-07` | D3 | `bse_test` 错误路径备份写成 `bin[1].star[1]`，应为 `bin[i]`；`nbin==1` 时越界读 | **已修改，待专项验证** | [bse_test.cxx](PeTar_new_version/bse-interface/bse_test.cxx)：bse_test 的第二颗星备份使用当前 bin[i]。 |
| `ASTRA-BSI-08` | D3 | `bse_test` 在 `omp parallel for` 内对共享 `std::vector` 做 `push_back`，并并发调 Fortran BSE 无 `critical` 锁 | **已修改，待专项验证** | [bse_test.cxx](PeTar_new_version/bse-interface/bse_test.cxx)：bse_test 双星演化串行调用共享 Fortran 状态，避免 OpenMP 公共块及容器竞争。 |
| `ASTRA-BSI-09` | D4 | `StarParameter::initial()` 不初始化 `lum`；`recordEvent()` 在首次 `evolveBinary` 前读它写初始事件 | **已修改，待专项验证** | [bse_interface.h](PeTar_new_version/bse-interface/bse_interface.h)：StarParameter::initial 显式初始化 lum。 |
| `ASTRA-BSI-10` | D4 | `mlwind.f` 的 `/VALUE1/` 声明 4 个 double（`neta,bwind,hewind,mxns`），其余为 3 个；`hewind` 从未使用 | **已修改，待专项验证** | [mlwind.f](PeTar_new_version/bse-interface/bse/mlwind.f)：mlwind 的 VALUE1 公共块移除不匹配的 mxns 字段。 |
| `ASTRA-BSI-11` | D4 | `zcnsts` 的 `msp(81)` 三次项写成 `lzs*lzs*xr(83)`，相邻多项式均为 `c0+lzs*(c1+lzs*c2)` | **待确认** | [zcnsts.f](PeTar_new_version/bse-interface/bse/zcnsts.f)：msp(81) 的三次项确实保留，但“相邻多项式都是二次”不足以证明系数拟合错误；需权威公式/参考实现或物理基准复核。 |
| `ASTRA-BSI-12` | D3 | `dgcore` 的 `temp` 计算后未用：`if(temp.gt.the)`/`if(temp.gt.tc)` 被注释，`x=1`、`y=0` 写死；He→C 无条件、C/O→Ne 永不发生 | **待确认** | [dgcore.f](PeTar_new_version/bse-interface/bse/dgcore.f)：dgcore 固定 x=1,y=0 与注释条件确实存在；属于模型选择，需预期核转换规则才能判 bug，不能仅因变量未使用就定错。 |
| `ASTRA-BSI-13` | D3 | `FBFAC=FBTOT/(mt-mcx)` 无上限；`nsflag=3` 薄包层 `mt-mcx<0.2` 时 `FBFAC>1`，`kick` 的 `vk*(1-FBFAC)` 反向放大 | **已修复（专项回归）** | [hrdiag.f](PeTar_new_version/bse-interface/bse/hrdiag.f) [hrdiag.f](PeTar_new_version/bse-interface/bseEmp/hrdiag.f)：经用户批准，bse/bseEmp 的 rapid/delayed 分支共用 bounded_fallback：先计算回落质量，再限制到 [0,包层质量]；零包层返回零质量/比例，代数消去薄包层除法。实际 Fortran 边界测试覆盖 84 组输入，普通公式两例一致。 |
| `ASTRA-BSI-14` | D4 | `kick.f` 在 `DISP0=0`（`--bse-sigma 0`）时 `v(4)=0`，`v(k)*vk/v(4)` 产生 `0/0=NaN` | **已修复（专项回归）** | [kick.f](PeTar_new_version/bse-interface/bse/kick.f) [kick.f](PeTar_new_version/bse-interface/bseEmp/kick.f) [kick.f](PeTar_new_version/bse-interface/mobse/kick.f)：经用户批准，bse/bseEmp/mobse 在抽样长度为 0 时输出零 kick，保留 RNG 调用次数；现有处方均对抽样长度作缩放。sigma=0/完全回落/普通非零 kick 回归通过；旧代码在同一零 kick 测试触发 SIGFPE。 |

### PeTar 配套（mcluster）

| 编号 | 原等级 | 问题 | 当前状态 | 当前依据与条件（含源码） |
|---|---|---|---|---|
| `ASTRA-MC-01` | D1 | 输出速度实为 pc/Myr，表头/配套工具按 km/s 消费，差 1.0227 | **仍存在** | [main.h](mcluster/main.h)：G=0.00449850214 配 rvir/tscale 给 pc/Myr，而输出仍标 km/s；含 kms2pcmyr 的下游链会二次转换。原生产脚本是否仍用此参数本次未确认。 |
| `ASTRA-MC-02` | 条件 D1 | `-Q` 不参与速度标定，只写进 Nbody6 的 ZETA | **仍存在** | [main.c](mcluster/main.c)：-Q 只进入记录/部分 Nbody 参数输出，归一化速度未乘目标 Q 因子；-C 3/5 且 Q非0.5 时缺陷保留。 |
| `ASTRA-MC-03` | 条件 D1 | BSE 周期缺总质量因子：`P=sqrt(a³/(m1+m2))`，`m` 按 `M` 归一 | **已修复（专项回归）** | [main.c](mcluster/main.c)：经逐项批准，送入 BSE 的周期使用实际太阳质量，返回半长轴使用演化后质量。链接包装实际 get_binaries 验证 1 AU、1+1 Msun 的周期为 365.25/sqrt(2) 天；质量减半且周期不变时 a 缩小为 cbrt(0.5) AU。 |
| `ASTRA-MC-04` | D2 | BSE 并合/解绑（`P=-1`、两子星无质量）未处理：`0/0` 崩溃或假轨道 | **仍存在** | [main.c](mcluster/main.c)：evolv2 返回后未先分流并合/解绑/双零质量，直接 P²重建轨道和质量比速度；负 P 被平方掩盖，零总质量除零路径仍在。 **用户明确选择暂不修改**（2026-11-07）；除零/假轨道风险仍保留。 |
| `ASTRA-MC-05` | D2 | 逃逸速度缺 `G`；双星判据用 `sqrt(vkick)` 且 `&&` | **已修复（专项回归）** | [main.c](mcluster/main.c) [main.h](mcluster/main.h)：经逐项批准，G 的单位明确为 pc³/(Msun Myr²)，逃逸速度转为 km/s，估算式补 G；直接比较 kick，任一分量超限即重抽。每次重抽恢复原始双星质心和恒星记录。仅第一组分超限的注入测试验证恰好重抽一次、质心不漂移。 |
| `ASTRA-MC-06` | D2 | `-f 2` 在 `epoch=0` 时 `R=L=0`（`evolv1_` 被 `if (epoch)` 包住） | **部分修复** | [main.c](mcluster/main.c)：去掉 generate_m2 对 evolv1 的 if(epoch) 限制，零年龄幂律质量函数的 R/L 均为有限正数。f=2 单质量快捷分支仍未走 SSE 结构初始化，因此本项不整体关闭。 |
| `ASTRA-MC-07` | D3 | `-f 2` + 单个 `-m`：质量写成 `1/N`，`Rh3D` 累计不到 0.5M → 除零 | **已修复（专项回归）** | [main.c](mcluster/main.c)：f=2 单质量分支以实际单星质量填 star[][0] 和初始质量字段，随后统一归一化；N=20、单质量 2 Msun 输出全部质量为 2，未再触发半质量半径除零。 |
| `ASTRA-MC-08` | D2 | 气体势重标：epot 仍 `i<j-1`，外势只算 0.5 倍；仅 `-C 3` | **仍存在** | [main.c](mcluster/main.c)：气体重标的星际势循环仍 i<j-1，且外势能缩放保留半系数；仅 output3 路径，历史 Q 数字未重新实测。 |
| `ASTRA-MC-09` | D2 | 内置 `mlwind.f` 与 `bse/mlwind.f` 不同（SSE 原始风 vs Vink+Belczynski） | **待确认** | [mlwind.f](mcluster/mlwind.f)：两份 mlwind 模型确实不同，但预演化与运行采用不同风模型是否违背项目约定需明确；质量差异不能单独证明实现 bug。 |
| `ASTRA-MC-10` | D3 | `-f 3` / `-f 4` 被解析器拒绝（只收 0..2），实现与帮助都声称支持 | **已修复（专项回归）** | [main.c](mcluster/main.c)：解析器接受帮助及实现支持的 0..4；使用 -M 20 -f 3 和 -N 20 -f 4 的实际生成测试通过。 |
| `ASTRA-MC-11` | D3 | DM 输出缺 `--with-interrupt=sse` 的 24 列布局 | **仍存在** | [darkmatter.h](mcluster/darkmatter.h)：DM 输出仍只区分普通/BSE列，没有 SSE 布局分支；与 SSE 构建读入契约不匹配。 |
| `ASTRA-MC-12` | D3 | DM 输出在质心系而非密度中心；`\|cm\|=1.7–2.2 pc > Rh` | **待确认** | [darkmatter.h](mcluster/darkmatter.h)：代码明确有意对恒星+DM统一平移到联合 CM，保持两分量相对中心；仅以原点测 Rh 改变不能证明分布错误，需密度中心是否必须固定原点的契约。 |
| `ASTRA-MC-13` | D3 | `epochstar` 恒 0，与 `-e` 演化后的质量/`kw` 矛盾 | **仍存在** | [main.c](mcluster/main.c)：star[][13] 仍统一写 epochstar=0，与预演化后质量/类型可不一致；需要 -e>0 的初态消费链。 |
| `ASTRA-MC-14` | D3 | 缺少 `CALL instar`，`/TYPES/ ktype` 全 0 → 并合产物 `kw=0` | **已修复（专项回归）** | [main.c](mcluster/main.c) [main.h](mcluster/main.h)：经逐项批准，main 在 SSE 构建中首次演化前调用 instar。碰撞表关键项验证为 ktype(0,0)=1、ktype(10,10)=15，SSE 初态生成冒烟测试通过。 |
| `ASTRA-MC-15` | D3 | >100 M☉ 的 100 M☉ 夹紧被启用 | **待确认** | [hrdiag.f](mcluster/hrdiag.f)：>100太阳质量的夹紧代码仍在；它可能是旧 SSE 模型适用域限制，需明确项目承诺支持更高质量，不能单用另一模型输出判实现错误。 |
| `ASTRA-MC-16` | D3 | GCC 15 下 Makefile 构建失败（`feenableexcept` 缺 `_GNU_SOURCE`；`output5` 无原型） | **已修复（专项回归）** | [main.c](mcluster/main.c) [main.h](mcluster/main.h)：main.c 在系统头前定义 _GNU_SOURCE，main.h 补 output5 原型；当前编译器的 gnu2x + implicit-function-declaration 设为错误的检查通过。未声称已使用 GCC15。 |
| `ASTRA-MC-17` | D3 | King/EFF 二分查找在端点可越界 | **仍存在** | [main.c](mcluster/main.c)：King/EFF 搜索端点仍无可靠边界限制；当前随机采样是否触发未复现，保留条件性越界风险。 |
| `ASTRA-MC-18` | D4 | 双星在 Q/Rh 归一化之后插入，此后不重标 | **待确认** | [main.c](mcluster/main.c)：插入双星后不重标是事实，但应区分质心系统维里比和包含内部束缚能的维里比；登记册约2%偏差不能单独判为 bug。 |
| `ASTRA-MC-19` | D4 | `fractalize` 用陈旧 `i==N` 做速度归一化，实际未归一化任何星 | **已修改，待专项验证** | [main.c](mcluster/main.c)：fractalize radial 分支逐个粒子归一化速度，避免循环结束后访问 star[N]；平均速度为 0 时不除零。 |
| `ASTRA-MC-20` | D4 | `strcat` 拼到 argv 字符串；`output4/5` 用 50 字节栈缓冲 | **已修复（专项回归）** | [main.c](mcluster/main.c)：输出文件名按输入长度配置缓冲，使用 snprintf，不再修改 argv；创建文件失败时明确报错。超过 90 字符前缀的 info/output3/4/5 测试通过。 |
| `ASTRA-MC-21` | D4 | "final energy check" 仍 `i<j-1`；`check` 无开关 | **已修改，待专项验证** | [main.c](mcluster/main.c)：最终能量检查补齐 i<j 的全部配对，新增 --check-energy 入口；与 MC-08 气体重标逻辑分开，后者尚未修复。 |
| `ASTRA-MC-22` | D4 | `energy_order` 动能项带 `1/Nstars`，排序几乎只按 `-1/r` | **待确认** | [main.c](mcluster/main.c)：energy_order 动能仍含1/Nstars；排序启发式的物理定义未独立核对，本次不把“几乎只按半径”直接认定为实现错误。 |
| `ASTRA-MC-23` | D4 | 潮汐半径的 `omega`（km/s/pc）当 1/Myr 用 | **仍存在** | [main.c](mcluster/main.c)：rtide 用 G 的 pc/Myr 制与 VG/RG 的 km/s/pc 混合；与 MC-01 同属单位体系问题，数值影响范围依参数。 |
| `ASTRA-MC-24` | D4 | `--blackhole` 的 `R/L` 与 SSE 的 `kw=14` 约定不同 | **仍存在** | [main.c](mcluster/main.c)：--blackhole 仍设置 R=0.01、L=0，与 BSE 紧致天体约定不同；后续 BSE 刷新可覆盖，属于输出初始恒星字段不一致。 |
| `ASTRA-MC-25` | D4 | DM 与恒星共用 `drand48` 流，`seed` 只写 `.dm.info` | **非缺陷** | [darkmatter.h](mcluster/darkmatter.h)：单一随机流在改变恒星数量后令后续 DM 序列变化是正常 PRNG 行为；未见“DM 样本必须独立于 Nstar”的承诺，不按 bug 计。 |
| `ASTRA-MC-26` | D4 | DM 的 `step>0.065` 分支不可达（最大 0.0586） | **非缺陷** | [darkmatter.h](mcluster/darkmatter.h)：前置参数限制使 grid step 守卫不可达，不等价于错误结果；保留为冗余检查信息。 |
| `ASTRA-MC-27` | D4 | 极端但合法 DM 输入触发 SIGFPE 而非报错 | **仍存在** | [main.c](mcluster/main.c)：正数解析不排除极小尺度，后续平方/比例下溢与启用浮点陷阱组合仍可异常退出；极端输入健壮性项，未复跑原命令。 |
| `ASTRA-MC-28` | D4 | 自带测试的编译命令在 GCC 15 下失效 | **已修改，待专项验证** | [main.c](mcluster/main.c)：源文件的 MC-16 声明修复同时覆盖测试编译入口；现代 C 严格声明检查通过，自带 DM 测试结果另记。 |

## 5. 本次验证与限制

本次没有修改正式源码、已有测试、构建配置或原登记册。测试产物/日志写入 `/tmp`，没有启动生产模拟或覆盖用户快照。

| 验证 | 本次结果 | 证明范围 |
|---|---|---|
| 登记册解析与逐项映射 | 192 个唯一编号全部覆盖，无漏项/重复 | 编号及模块分类完整性，不是192次动态复现 |
| `bash PeTar_new_version/test/bugfix/run.sh` 的 core 阶段 | 通过，ASan+UBSan | 3–6粒子树边界、6叶拓扑复制、ARI2/SYS2往返、List赋值、块步长0位图、changeover、EnergyAndMomentum序列化 |
| 同脚本 integration 阶段 | **中止，退出码134**，异常 `Coincident escape candidates without softening` | 此前顺序执行的 ASCII deadline/打印、邻居平移/ID表、SIMD堆缓冲检查已通过；不能声称整套 integration 通过，BSE单测在该脚本中未到达 |
| 单独编译/运行 `test/bugfix/python_bridge.cxx` | 通过 | remnant异常/坏数组、有向角、近共振恢复及引用管理；使用注入模块，不等于真实完整进动物理验收 |
| `python3 PeTar_new_version/test/bugfix/precession_sign.py` | 通过 | 本机实际 Python 库 cyclesign 回归 |

集成失败原因已核对：[sampleParticle](PeTar_new_version/test/bugfix/integration.cxx:5) 把所有测试粒子位置设为 `(2,0,0)`，旧测试依赖缓存零势；B45 后 [escape_removal.hpp](PeTar_new_version/src/escape_removal.hpp) 使用当前直接势能并拒绝无软化重合点。因此该旧夹具与新删除算法不兼容，不能仅据本次中止断言 EN-02/EN-05 又回归，也不能忽略测试已失败的事实。

本次日志：`/tmp/astra_bugfix_validation.log`、`/tmp/astra_python_bridge.log`、`/tmp/astra_precession_sign.log`。这些是临时文件，关键结果已写入本文。B45 的1/2/4 rank历史验证见 [专项调试日志](diagnostics/escape_removal_fix/debug.log) 和 [专项测试说明](PeTar_new_version/test/escape_removal.md)；本文核对了新记账公式，**未重新执行多rank专项测试或完整星团长程运行**。旧基线探针、Helgrind、全部非默认宏、GCC15、未来NumPy版本亦未全量重跑。

## 6. 未编号候选和撤回项

原登记册第11节已撤回的 N7、`vkick[4..7]` 未初始化等不重新加入192项。第15节尚未亲核候选不自动升级为已确认 bug：已编号的 `DG-09/10/11..14`、`BL-07/10/11`、`GEN-02..05` 等已在主表覆盖；剩余如伪粒子 `e>=1` 输入守卫、cluster 混合CM/成员判据、MPI计时归约开销等，需要独立调用路径和反例，本文不另造编号或沿用旧文中的执行建议。


## 7. 本轮修复安排（2026-11-07，进行中）

用户要求：修复已确认项，未确认项保留；可能改变并合逻辑的修复先征求意见。B 类已获整体授权；A 类逐项批准 IO-01、BSE-06/07、BSI-01/03/13/14。其余 A 类继续保留待确认，未把等待当作同意。下列是原设计清单，落实情况见第 8 节。

### 待决定 A：直接影响并合路径

- IO-01：将 ASCII 后缀升级为有版本的完整预测状态，保存 deadline、轨道检查点时间/a/e、累积动力学量；旧格式继续兼容并明确无法逐事件恢复。新格式重启首次检查不重建已有效的截止时刻。
- BSE-06/07：CE/MIX 后用新的恒星类型、质量和年龄重新计算结构参数，再回填半径/光度/包层/自旋相关量。可能改变下一次接触判定，必须以角动量守恒方式更新自旋而非任意重设。
- BSI-01/03/13/14：分别修正接触潮汐平方位置、保留 CE 超新星 kick 后轨道、限制 fallback 到物理范围、处理零 kick。会改变演化/反冲，逐项回归后才能接受。
- BH-05/08/09/10/11：接受圆轨道合法输出形状；仅裁剪舍入量级的负 e²；部署或定位 precession 包；对端点奇点做有限性处理；修复宿主 Python 已初始化的接口。可能将旧回退 BSE/直接拟合切换为成功的进动拟合，故需确认后实施。
- AR-08/09：修正并合后 ARRAY 漂移尺度和零质量节点轨道状态；可能改变并合后的组拆分/成员轨迹。
- MC-03/04/05/14：修正初态双星预演化的单位、并合/解绑输出、kick保留及碰撞类型初始化，可能改变交给 PeTar 的并合产物/双星名单。

### 待决定 B：间接改变并合历史的动力学修正

TR-11 将邻居判据的 G=1 改为运行时 G；TR-10 修正 Kepler Newton 的 cos 符号；AR-10/12、H4-06、TR-17/19 修正时间同步、扰动尺度及非默认力核；MC-01/02/06/07/08/23、GEN-01/04 修正初态单位/归一化与抽样。它们不增加新的并合判据，但可能改变相遇和最终并合历史。只有用户明确同意这一类后才实施；不把初态改变包装成纯诊断修复。

## 8. 第一批实施快照（2026-11-07；最新增量见第 9 节）

本轮涉及 **55 个原编号**（含重复根因），并非 55 个独立补丁：9 项已有针对性回归，45 项已修改但仍需各自触发条件下的专项验证，FD-17 为部分修复。默认构建和通用回归不能替代每个非默认宏、MPI 分支及生产路径的验证。逐项表中的“已修改，待专项验证”不计为完全关闭。

未确认的 18 项全部保留，未据猜测修改。尚未处理的已确认项也仍在表中，不代表本轮已完成全表修复。前面第 2 节统计、第 3 节优先级、第 5 节旧测试记录是修复前审计快照。

授权记录：B 类可实施并记录数值影响；A 类已批准 IO-01、BSE-06/07、BSI-01、BSI-03、BSI-13、BSI-14。BH、AR 并合后状态及 mcluster 预演化等其余 A 类未获逐项批准前不实施。具体补丁、验证范围和数值变化见 [BUG_FIX_SUMMARY.md](BUG_FIX_SUMMARY.md#astra-登记册修复进展2026-11-07)。

### 下一批 A 类具体待确认方案

- BH-10（分阶段）：先处理已明确的角度端点舍入越界：余弦只在距离 [-1,1] 为舍入误差时裁剪，明显越界保留失败；不把所有非有限量强行置零。真正零自旋及未证明为重根的退化解继续保留未关闭，后续另列方案。
- BH-11：分开记录“宿主 Python 已存在”和“本桥接 NumPy API/模块路径已初始化”；宿主解释器下在持有 GIL 时完成本模块初始化，退出时不 finalize 宿主解释器。测试独立启动与宿主预初始化两种生命周期。
- AR-08：仅在 AR_SLOWDOWN_ARRAY 的单存活粒子分支，将漂移时间乘当前 slowdown 的逆因子，与同宏下的一般漂移公式一致；TREE 分支另按现有约定运行。以 κ=1 和 κ>1 的位移对照验证。

- MC-03：一般双星分支计算交给 BSE 的周期时，将归一质量乘回总质量 M，与 OB 分支及以 AU/太阳质量/年表示的 Kepler 关系保持一致。预演化分离和后续并合历史会改变。
- MC-14：在 mcluster 第一次调用恒星/双星演化前执行一次 instar，初始化 COMMON /TYPES/ 的碰撞产物类型表；不修改表中处方。可能将原来错误的 kw=0 产物改为表规定的类型。

### 剩余 A 类待决定方案

- AR-09：追加零质量成员时，先清空人工节点的轨道缓存并重新计算成员质量、质心和当前分离；只在轨道定义良好时重算轨道。退化节点不再沿用上一棵树的 semi/ecc/r/stab，并明确作为非正常物理双星处理，保留内部正质量子树。需覆盖 4 正质量→2 正质量+2 零质量，以及单存活、全零质量、重合成员场景。
- MC-05：逃逸速度与 BSE kick 统一为 km/s，补上估算式缺失的 G；直接比较 kick 大小而非再次开方，任一组分超限就拒绝本次“双方均保留”的抽样。它会改变保留比例及后续双星名单。
- MC-04：现有 mcluster BSE 接口仅返回 kick 大小，不能可靠还原解体后两星的完整相空间。可先在重建轨道前明确拒绝 P<=0、非束缚偏心率及双零质量等不支持返回，避免生成假双星/除零；完整产物转换需要扩展演化接口并明确存活星/解绑星输出契约，不能用随机方向冒充原演化结果。

## 9. 后续授权与增量修复（2026-11-07）

A 类新增批准 BH-05/08/09、BH-10 的端点舍入阶段、BH-11、AR-08/09、MC-03/05/14；MC-04 用户明确选择暂不修改。BH-10 的零自旋/重根问题仍保留；没有以本阶段端点回归替代全部退化自旋验证。

当前实施清单共覆盖 **102 个原编号**，既有“专项回归”也有“已修改待专项验证”和部分修复，不等于同数量独立补丁或全表已关闭。逐项表为当前状态；第 2、3、5 节和第 8 节为较早审计/实施快照。未确认的 18 项保持原分类，剩余未处理的已确认项继续保留。

新测试入口包括 `bash mcluster/tests/run_astra.sh`；主回归新增真实 ARRAY integrator 的并合注入、宿主 Python 生命周期及端点舍入回归。模块路径选择顺序为显式 PETAR_PRECESSION_PATH、构建时绝对源码路径、编译时安装路径；显式路径无效会失败，不默默换包。安装规则可用 DESTDIR 暂存，当前仅在 /tmp 验证，没有部署到用户安装目录。

### 进动库辅助接口的下一批逐项方案

这些条目多在非默认辅助接口中；修复可能让原先报错的路径返回有效拟合输入，按“可能改变并合逻辑先询问”的要求单独确认。

- BH-15：两处预计算根解包直接接收三行 precomputedroots，去掉会删掉第三行的 `[:-1]`；检查形状并对照不传预计算根的结果。
- BH-16：widenutation_condition 的返回改为调用 `np.stack((which,kappa,chieff))`；保持各分量的顺序和定义。
- BH-17：为 eval_costheta12 的守恒量分支补充显式 r 参数，角度分支保留原调用兼容性；eccentricity 包装器相应转换 a/e，并验证两种输入形式一致。不凭空选择轨道分离。

### FDPS 与输出修复增量

FD-03/14/15 的部分修复和剩余限制已写回逐项表；FD-04/22/23 专项回归通过。IO-08 增加输入文件保护，DG-04 延后表头，正在运行输出专项回归。前一轮默认构建及完整 bugfix 回归已通过，包括 ARRAY 漂移、恒星结构、fallback 和宿主 Python；未据此声称全部 MPI/非默认配置通过。

### 后续 A 类具体方案（尚未批准，未实施）

- BH-21：eval_S 的互斥输入分支改为 if/elif，非法输入显式 raise TypeError；eval_kappa 同样显式抛错；eval_q 接受合法等质量 q=1。仅修复这些辅助接口，不声称消除其他演化式的等质量奇点。
- BH-22：无固定 chieff 的几何 deltachi 边界改为 ±(chi1+q*chi2)/(1+q)，用两自旋反向极点和角度采样检查包络。会放宽原来错误缩窄的允许区间。
- BH-23：deltachitildeav2 在小 m 区间使用由原公式展开的级数（常数项 3/8、一次项 m/16），避免椭圆积分相减再除 m² 的消去误差；对照高精度计算验证切换点，保留其他端点策略。可能改变平均自旋演化。

本轮输出专项与完整 `test/bugfix/run.sh` 已通过：硬链接输入快照保留原内容，下一编号照常输出，`-w 2` 表头包含粒子列。默认 PeTar 构建通过。mcluster 普通构建在固定种子 L3 测试中间歇 SIGSEGV，ASan/UBSan 构建连续 20 次未复现；仍在定位，不能宣称 mcluster 全部回归通过。

### 最新范围调整与验证

用户明确要求停止处理 mcluster；自此不再排查、修复或运行其测试。此前修改和失败记录保留，不据此追加修复，也不将尚未通过的结果改写为通过。工作范围继续为 SDAR、PeTar、FDPS。

BH-21/22/23 已逐项批准并落实。BH-15/16/17 尚无批准回复，保持待确认。GEN-03 已用 NDEBUG、单进程探针验证明确中止。新增 KDK/DKD 测试覆盖默认与 KDKDK_4TH 系数。

最新增量包含日志故障处理、GEN-01/04 生成器及 BL-08 参数检查。BL-13 的 configure 生成链已重新核对并澄清，未把合法配置期生成流程误写成完全没有生成规则。

本批收尾验证：`make -C PeTar_new_version -j2` 和 `bash PeTar_new_version/test/bugfix/run.sh` 均退出 0；额外纯 z 方向 CM 预测探针通过，并已加入测试入口。日志中的 “Failed to write merger event log” 是故障注入测试的预期报错，测试已验证半行回滚。`git diff --check`（兼容原 CRLF 文件设置）通过。

本轮累计涉及 102 个原编号，包含部分修复、文档澄清及仍待专项验证项，**不代表 102 项全部关闭**。未确认项、未获批的 A 类及剩余已确认项仍在分析表；mcluster 已按最新指示停止处理。MPI 多 rank 日志故障、全部非默认构建及长期星团演化未全量验证。
