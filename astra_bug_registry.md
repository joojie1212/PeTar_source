# astra 缺陷登记册（唯一编号来源；按 D 等级排序）

署名：**astra** 审阅体系 + 主代理复核
日期：2026-10-27（2026-10-27 按 D 等级重排）
基线：`PeTar_source@fbebd9c`（`robert/perf-main`）；正式源码未改。

本文件是**所有缺陷条目的唯一编号与等级来源**。此前五轮审阅使用了四套互相冲突的编号（主报告 `A1–A6/B1–B10/C1–C4/D1–D2/E1–E4`、[第一轮 BSE 报告](astra_physics_state_audit.md)的 `B1–B5`、[第一轮 SDAR 报告](astra_sdar_silent_error_audit.md)的 `A/B/C`、[第三轮](astra_round3_sdar_state.md)的 `R3-S1/S3-x/E3-x`），同一符号在不同文件指向不同缺陷。自本版起统一改为**按机制分区的稳定编号** `ASTRA-<AREA>-<NN>`，旧编号仅作别名保留。

**排序约定（本版更改）**：正文按**破坏性等级 D1 → 条件 D1 → D2 → 条件 D2 → D3 → D4 → 契约/非缺陷**排列，同一等级内按分区 `AR → H4 → BSE → BH → TR → FD → IO → EN → BL → DG → GEN → MC → BSI`、再按编号数字升序。**审阅轮次不再是组织顺序**，只保留在每行末尾的“来源/旧编号”列与[§14](#14-来源轮次索引不参与排序)里，供追溯用。编号仍按发布顺序分配，因此编号大小与严重程度无关。

阅读入口：[目录总览](README.md) → [索引](astra_review_index.md) → [评级与本册](astra_bug_severity.md) → 各轮专项正文（正文仍用旧编号，请按[对照表](#12-旧编号对照)换算）。

## 1. 编号规则

```
ASTRA-<AREA>-<NN>
```

| AREA | 机制范围 | 主要文件 |
|---|---|---|
| `AR` | AR 力核、递归事件状态、AR 组调度与 AR 序列化 | `src/ar_interaction.hpp`、`SDAR/src/AR/*` |
| `H4` | Hermite 积分器、块步长、组调度、H4 粒子 | `SDAR/src/Hermite/*` |
| `BSE` | BSE 边界时间步、单位、并合互斥、Fortran 输出 | `bse-interface/*`、`bse/*.f` |
| `BH` | BH 并合 remnant 拟合与 Python 接口 | `bhmerger/*` |
| `TR` | 树/邻居/FDPS/稳定性/changeover/SIMD 目标筛选 | `FDPS/src/*`、`src/cluster_list.hpp`、`src/stability.hpp`、`src/changeover.hpp`、`src/soft_force.hpp` |
| `FD` | **FDPS 框架本体**（树/力/LET/域分解/容器），第八轮只覆盖默认构建编译面 | `FDPS/src/tree*.hpp`、`particle_system.hpp`、`domain_info.hpp`、`ps_defs.hpp`、`memory_pool.hpp`、`reallocatable_array.hpp`、`sort.hpp` |
| `IO` | 检查点、序列化、重启契约 | `particle_base.hpp`、`AR/information.h`、`petar.hpp`(restart) |
| `EN` | 总能量账本与删除能量 | `src/petar.hpp`、`src/energy.hpp`、`src/hard.hpp`(能量) |
| `BL` | 构建/宏/CLI 接口可用性 | `Makefile`、`Common/list.h`、`src/petar.hpp`(参数) |
| `DG` | 计时、诊断、打印 | `time_profile_*`、`*::print`、并合日志 |
| `GEN` | 初态生成/域分解/伪粒子多极 | `src/particle_distribution_generator.hpp`、`src/domain.hpp`、`src/pseudoparticle_multipole.hpp` |
| `MC` | **初态生成器**：剖面采样、IMF/双星、维里归一化、速度单位、输出列布局 | `mcluster/main.c`、`mcluster/darkmatter.h`、`mcluster/mlwind.f` |

编号一经发布**不复用、不回收**：新缺陷在所属 AREA 分区内取下一个编号，但**正文位置由等级决定**（插入对应 D 级小节），不追加到文末。撤回条目移入[第 11 节](#11-已撤回条目不占编号)，编号同时作废。

## 2. 等级规则（按代码事实判定）

| 等级 | 判定条件（三条同时满足） |
|---|---|
| **D1** | ①机制上直接改变粒子位置/速度/质量/自旋、事件时刻、相互作用分组或存活成员；②路径在**默认构建 + 当前生产运行时配置**下可达；③有数值反例或真实消费端证据 |
| **条件 D1** | 满足 ①③，但路径依赖特定事件组合或运行时配置；**必须在条目里写明条件**，不得读作"当前生产已发生" |
| **D2** | 未初始化值、数据竞争、单位/参考系不一致使**数值控制或总能量账本不可靠**；路径在默认构建或受支持的运行时配置下可达；尚无轨迹污染实证 |
| **D3** | 仅在**非默认宏**或**特定构建目标**下可达的 UB/编译失败/内存安全，或仅影响辅助接口、输入解析、运行可用性；若存在"启用即崩溃/错值"的消费路径，标注**潜在 D2** |
| **D4** | 已核对仅影响计时、分类统计或打印，物理值与账本不受该缺陷影响 |

配套判定约定（本轮新增，消除此前的字面冲突）：

1. **力核越界不按 D3 处理。** 力/势计算路径上的越界或未定义值，即使只在非默认构建可达，也记 **D2**（数值本身不可信），并注明构建条件。
2. **"需非默认宏才可达"的 UB 默认记 D3（潜伏）**，不记 D2；一旦该宏在生产构建中被启用，按 ①–③ 重新评。
3. **限制在特定运行时选项的未初始化读**（不需要任何编译宏）：若该选项是受支持的生产配置，记 **D2**；若只在实验/调试选项下可达，记 D3。
4. **契约类缺口**（如重启随机状态）不强行套 D1–D4，单列 `契约`，写明在哪种契约下等价于 D1。
5. **非缺陷但需记录**的项记 `非缺陷`，附代码事实，避免反复被重新报告。

## 3. 等级分布（共 192 项）

| 等级 | 项数 | 排序口径 |
|---|---|---|
| **D1** | 2 | 默认构建 + 当前生产配置可达，直接改物理状态 |
| **条件 D1** | 9 | 同一机制，但需特定事件组合/运行时配置 |
| **D2** | 24 | 数值控制或总能量账本不可靠（含“潜在 D1”） |
| **条件 D2** | 3 | D2 机制，仅特定配置/启动方式可达 |
| **D3** | 98 | 非默认宏/特定构建目标；辅助接口/可用性 |
| **D4** | 54 | 仅计时、统计或打印 |
| **契约** | 1 | 不套 D1–D4，注明等价条件 |
| **非缺陷** | 1 | 记录在案，避免重复上报 |
| **合计** | **192** | |

同一等级内不再按“第几轮”排列；需要按轮次追溯时用[§14](#14-来源轮次索引不参与排序)。

## 4. D1（2 项）

判定：机制上直接改变粒子位置/速度/质量/自旋、事件时刻、相互作用分组或存活成员，且路径在**默认构建 + 当前生产运行时配置**下可达，并有数值反例或真实消费端证据。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-AR-01` | AR 组并行下，他组质量事件对扰动源的可见性依赖调度次序 | `ar_interaction.hpp:473,491,546-549,778-790`；`hermite_integrator.h:2478-2537` | 邻居表存 `&group.particles` 真实地址（`neighbor.h:181`），phase 1 即改 `cm.mass`；真实 SDAR 两双星回归 `Δvx≈1.00038e-7` | **D1** | A1 |
| `ASTRA-MC-01` | 输出速度实为 pc/Myr，表头/配套工具按 km/s 消费，差 1.0227 | `main.h:6`；`main.c:1191,1204,1214,5114` | `G=4.49850214e-3` 是 `pc³/(Msun Myr²)`（等价 `4.30092e-3×1.0227²`，逐位吻合）。端到端 `mcluster -u 1 -C 5` → `petar.init -v kms2pcmyr` → **Q=0.5230**；去掉 `-v` → **Q=0.5000**；生产脚本原参数（`test6.sh:149-150`，N=10000、30 M☉、Rh=0.0013 pc）同样 **Q=0.5230** | **D1** | 第十轮 |

## 5. 条件 D1（9 项）

判定：满足 D1 的机制与证据，但路径依赖特定事件组合或运行时配置。**每条都写明条件**，不得读作“当前生产已发生”。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-AR-02` | 递归 `modify_return` 取子分支最大值，局部 destroy 被提升为整子树 destroy | `ar_interaction.hpp:944-980`、`:1052-1057` | 五体注入 `return=3 status=3 live_mass=1`；真实消费端清零 CM 并跳过存活者漂移（应漂移 `0.000125`） | **条件 D1**（≥5 成员嵌套树 + 同次多分支毁灭） | A5 |
| `ASTRA-AR-03` | 束缚↔非束缚切换时复用旧 `ecca`，`calcParticles` 只按 `semi` 符号选公式 | `ar_interaction.hpp:1145-1151,1190-1200,1960-1990`；`two_body_tide.hpp:168-178` | canonical `COMM::Binary`：束缚→双曲 `r` `0.617579→0.882754`；反向 `1.674870→0.372964`；按当前 `r` 反解 `ecca` 后恢复 | **条件 D1/D2**（BSE 解体无 kick；`stellar-evolution=2` 潮汐交接） | N1 |
| `ASTRA-BH-01` | `BSEManager::merge` 把 Msun 质量与内部单位 G/dr/dv 混用 | `bse_interface.h:1739-1743`；`inspiralZJ.cpp:309-315` | 同物理事件、同 RNG，`mscale=1/2` 的 remnant 质量/自旋/kick 全部改变 | **条件 D1**（`BHMERGER_INTERFACE_LINKED` + BH-BH + `mscale!=1`；当前 `mscale=1`） | A6 |
| `ASTRA-BH-02` | remnant 方位角用无向 `acos`，丢绕角动量方向 | `inspiralZJ.cpp:341-345` | C++→Python 真实调用：第二颗自旋 y 取反，返回方向完全相同；有向角原型恢复符号 | **条件 D1**（BH 并合直接拟合路径） | A2 |
| `ASTRA-TR-01` | F32 邻居计数不具平移稳健性，真实连接被算成 `n_ngb==1` 而误归 single | `phantomquad_for_p3t_x86.hpp:28-36,135-142,165,550-595`；`cluster_list.hpp:454` | 平移 1e8 后 SIMD 计数 1、double scalar 与真实 FDPS 列表 2；double 对照恢复 | **条件 D1**（大绝对坐标 + 小间距；`use_simd_64=no` 下路径生效，四快照 39076 粒子未触发） | A4 |
| `ASTRA-TR-11` | 簇搜索两体判据把 `G` 写死 `1.0`；`-u 1` 生产配置下运行时 `G=4.4983e-3` → `peri` 小 222–250 倍，判据放宽 | `petar.hpp:892-893`；`cluster_list.hpp:400` | 真实函数 A/B（`probe_clv01.cpp`）：可达点 3/22 翻转、全网格 29/63，方向全为 G=1 多收；**真实快照 A/B**（`bench/th_0.7/input`，10,000 体、`-u 1`、`G=4.4983e-3`、`r_out=1.518e-5`）：`in_range=135`、判据通过 104→115、连通分量 **9896→9885**（11 个粒子被错误并组）；条件：v≳3σ 的退行邻居、`r∈(1.5r_out, r_search]` | **条件 D1** | 第十二轮 |
| `ASTRA-IO-01` | ASCII 检查点不保存 `time_merger/merger_semi/merger_ecc`，重启按 `now+lifetime` 重算 | `particle_base.hpp:168-208`；`ar_interaction.hpp:1870-1881` | 真实 write/read + 预测函数：deadline `976.5625→1708.984375`；binary roundtrip 保留原值 | **条件 D1**（ASCII 续跑且存在已安排的并合） | A3 |
| `ASTRA-MC-02` | `-Q` 不参与速度标定，只写进 Nbody6 的 ZETA | `main.c:287-289`、`938-949`、`1003-1013`、`4885` | `-Q 0.1` 与 `-Q 2.0` 的 `.txt` 逐字节相同；独立复算 Q 恒 0.5（code 单位）。`-C 3/5` 无消费者 | **条件 D1**（`-Q≠0.5`） | 第十轮 |
| `ASTRA-MC-03` | BSE 周期缺总质量因子：`P=sqrt(a³/(m1+m2))`，`m` 按 `M` 归一 | `main.c:3863`（对照 `3716`）；`evolv2.f:269` | 插桩：同函数 OB 分支写 `(m1+m2)*M`；`M=500` 时 BSE 实际间距 = 输入 × **7.9370 = 500^{1/3}**（3 例一致） | **条件 D1**（`-e>0` 且 `-b>0`） | 第十轮 |

## 6. D2（24 项）

判定：未初始化值、数据竞争、单位/参考系不一致，使**数值控制或总能量账本不可靠**；路径在默认构建或受支持的运行时配置下可达，但尚无轨迹污染实证。标注“潜在 D1”的条目仍归此级，是否升级取决于后续触发证据。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-H4-01` | `initial_step_flag` 在力循环内部被写，跨 item 读取无保护 | `neighbor.h:173,197`；`hermite_integrator.h:735,991-994` | Helgrind 对真实 `Neighbor` 方法报 2 个读写竞争；标志参与 `correctAndCalcDt4thList` 的步长重初始化（`:554-563`） | **D2**（未观测最小回归的数值差异） | B3 |
| `ASTRA-H4-03` | `calcNextDtLimit` 的 `bitmap==0` 导致 `while` 不终止 | `block_time_step.h:160-180` | `t=nextafter(dt_min,0)`，`t/dt_min=0.99999999999999989`，3 s `timeout` 返回 124 | **D2**（`time<dt_min` 时 H4 组步长挂起） | S3 |
| `ASTRA-BSE-01` | 紧致天体分支读取调用方未初始化的 `dtm` | `bse_interface.h:1841-1859`；`bse/trdot.f:77-81` | 槽值控制：0→0.1、0.2→2、7→70、100→500 Myr；Valgrind origin 指向 `getTimeStepStar` 栈 | **D2**（潜在 D1；tree 上界可能掩盖；NaN 时该粒子可能永不再触发 BSE，未证实） | B1 |
| `ASTRA-BSE-02` | 双星步长候选混用 IN 与 Myr 并重复除 `tscale` | `bse_interface.h:1874-1905` | 同物理状态 `tscale=1/2/0.5`：`308.88 / 154.44 / 617.77` Myr | **D2**（`tscale!=1` 才触发；本次历史全量 `tscale=1`） | B2 |
| `ASTRA-BSE-03` | 动力学 `merge` 直接调 `merge_`，未持 Fortran 共同锁 | `ar_interaction.hpp:1629/1633,1754`；`comenv.f:148,248,259`；`ran3.f:26-40` | `--wrap=merge_` 跟踪：`named_critical_depth=0`；可达共享 `/RAND3/` 写 | **D2（潜在 D1）**（需并发另组合并；未制造真实并发数值污染） | B4 |
| `ASTRA-BSE-04` | runtime 关闭恒星演化时 `modify_branch[2]` 未初始化即参与判定 | `ar_interaction.hpp:950,973` | 三体树（左叶子）、`stellar_evolution_option=0`，Memcheck 定位该栈变量；scratch 改 `{0,0}` 后 0 错误 | **D2**（不需任何宏即可达；正文曾建议 D3，按规则 3 定为 D2；未观测返回状态/轨迹差异） | B5 |
| `ASTRA-BSE-05` | 异步 `evolveBinary` 复用 `_out.dm` 作阶段性缓存，输出错的 `dm` | `bse_interface.h:1521-1522,1477-1486,1616/1622` | 真实接口 `actual_dm≈−2.631e-7`、`reported_dm≈+0.999999997`；最终质量正确 | **D2（潜在 D1）**（常规 AR 入口先同步 `time_record`，可达性未证） | B10 |
| `ASTRA-BSE-07` | `merge_` 回写的 `r/mc/rc/ospin` 是合并前遗留值，`lum` 完全不更新；`r` 被 `mergerCheckRadius` 消费 | `bse_interface.h:1754,1771-1784`；`merge.f:45-62`；`comenv.f:350-352`；`ar_interaction.hpp:890,1045-1046,70-74` | 真实 `libbse` 探针：CE 弹出后 `r=36.1177` vs 刷新 `0.102384`（352.8×）、`lum` 625.4 vs 13.34（46.9×）；MIX 直并 `r` 偏小 20.2%、`lum` 偏小 17.6× | **D2（潜在 D1）**（`merge_` 路径可达；CE 弹出/MIX 恒星并合影响显著，BH-BH 的 `r≈0` 影响小；需后续近遇才改变存活成员） | 第六轮 |
| `ASTRA-TR-03` | `generateBinaryTree` 的 `bin_host[ki]` 先比较后判界 | `binary_tree.h:699-704,716-721` | ASan `dynamic-stack-buffer-overflow` 于 `:703`（位置 `0,1,10`）与 `:720`（`0,9,10`），3 粒子链即触发 | **D2（潜在 D1）**（默认组树热路径 UB） | N2 |
| `ASTRA-TR-04` | `copyDataIter` 每子树只把 `inow` 减 1，右子树覆盖左子树 | `binary_tree.h:744-746,909-930`；`hard.hpp:2803` | 手工 6 粒子树 `0..5` 复制成 `4 5 2 3 4 5`，根 `getMemberN()` 仍 6；随机 6 粒子复制后 7 叶 | **D2（潜在 D1）**（`binary_table` 当前只填不读） | S1 |
| `ASTRA-TR-05` | quad + `AVX_PRELOAD` 下 `spjbuf[j+1]` 越过 `NJMAX` | `phantomquad_for_p3t_x86.hpp:1058` | 真实 PG 类 UBSan：`index 131072 out of bounds for type 'float [131072][3][4]'`；EP 路径不报但指针式预读未证安全 | **D2**（按规则 1；当前 `use_quad=no` 且不带 `AVX_PRELOAD`） | B6 |
| `ASTRA-TR-12` | `id_to_adr_pcluster_[...]` 用 `operator[]`，缺键静默插入 0，把边接到 `ptcl_cluster_[0][0]` | `cluster_list.hpp:592-596` | 相邻同类查表（`:819-821,826-838`）均有 assert；此处无 | **D2** | 第十二轮 |
| `ASTRA-IO-02` | `AR::Information` 把 8 字节 `ds` 按 `sizeof(int)` 写读 | `AR/information.h:91,368,378` | `ds_out=0.001234567890123`、`ds_in=3.7812373849316075e-315`，相对误差 1 | **D2**（检查点/续跑；`ar.cxx load_flag` 读后不重算） | S2 |
| `ASTRA-EN-01` | 孤立粒子演化阶段对四个共享能量字段裸 `+=` | `hard.hpp:2178-2200` | 源码确定；P1 两阶段设计可结构性消除 | **D2**（总账本不可靠；活跃演化事件动态对照仍缺） | B7 |
| `ASTRA-EN-02` | 逃逸删除在 OpenMP 循环内对四个共享能量字段裸 `-=` | `petar.hpp:2593-2620` | **已复核**：1 线程 `8`；8 线程三次 `99520 / 113824 / 111720`（应为 `8`）；scratch 线程归约后三次均 `8` | **D2**（单节点即可触发） | B8 |
| `ASTRA-EN-05` | 同一步删除多个逃逸粒子时逐粒子扣除 `m_i Φ_i`，把被删集合内部互势扣两次 | `petar.hpp:2583-2597,2605-2620`；`energy.hpp:238-241` | 两粒子互势 −1、同时逃逸：`ref_after=1`、真实末态 0、`ledger_error=-1`；单粒子对照 0 | **D2**（能量账本；单粒子路径不受影响） | 第六轮 |
| `ASTRA-MC-04` | BSE 并合/解绑（`P=-1`、两子星无质量）未处理：`0/0` 崩溃或假轨道 | `main.c:3870,3918,3933` | gdb：`i=1, m1=m2=0, P=-0.99997` → SIGFPE @ `3918`；`-N 20 -f 0 -m 20 -b 1 -e 10` 稳定复现（exit 136）。`P=-1` 而质量非零 → `a≈0.23 AU`、速度 ~1e6 的静默坏状态 | **D2** | 第十轮 |
| `ASTRA-MC-05` | 逃逸速度缺 `G`；双星判据用 `sqrt(vkick)` 且 `&&` | `main.c:1350,1489,1739,3669,3927` | `-N 1000 -R 1` 打印 28.2843 km/s，真值 1.855 km/s（15.2×）；100 km/s 的踢被当成 10 与 28 比较 → 几乎不拒绝 | **D2**（`-e>0`） | 第十轮 |
| `ASTRA-MC-06` | `-f 2` 在 `epoch=0` 时 `R=L=0`（`evolv1_` 被 `if (epoch)` 包住） | `main.c:1541`（对照 `1429/1773`） | `mcluster_sse -f 2 -e 0`：每行 `kw=1, R=0, L=0`；`-f 1` 同参数给出 `R=0.269, L=0.0093` | **D2** | 第十轮 |
| `ASTRA-MC-08` | 气体势重标：epot 仍 `i<j-1`，外势只算 0.5 倍；仅 `-C 3` | `main.c:5081,5091` | 子审：`-x 1000 -x 2 -x 1 -x 0`、N=1000 → `Q_self=1.017, Q_tot=0.332`（应 0.5）；`vscale` 1.42597 vs 1.75043 | **D2** | 第十轮 |
| `ASTRA-MC-09` | 内置 `mlwind.f` 与 `bse/mlwind.f` 不同（SSE 原始风 vs Vink+Belczynski） | `mcluster/mlwind.f` vs `bse-interface/bse/mlwind.f` | 子审：60 M☉→BH **7.683 vs 11.071 M☉（−30.6%）**；80: 10.195/13.538 | **D2**（`-e>0`） | 第十轮 |
| `ASTRA-BSI-01` | RLOF/接触分支对流阻尼因子写成 `ttid/(2tc)**2`（量纲 1/时间），detached 分支在下一行改正、接触分支未改 | `bse/evolv2.f:2009`（对照 `:600-601`） | 源码对照 + 量纲；`f_wrong/f_correct=1/ttid`，10 天轨道 `f` 偏大 36.5×、1 天 365×、100 天 3.65× | **D2** | 第十一轮 |
| `ASTRA-BSI-02` | `gntage` 对 `mc<mcagbf(Mhef)` 强制 `mmin=Mhef`，非单调 `mcheif` 的低质量根被排除，且无 `m0<=mt` 校验 | `bse/gntage.f:113-118` | 探针 `kw=4, mt=0.8, mc=0.477, fage=0.5` → `m0=1.844`（`m0/mt=2.31`）、`aj=1107 Myr`；4 组输入均返回 `m0≈1.8–2.8` | **D2** | 第十一轮 |
| `ASTRA-BSI-03` | CE 非并合分支在 `kick()` 后用 kick 前 `EORBF/ECIRC` 重算 `ECC`；MOBSE 有 `SNexplosion` 守卫，`bse`/`bseEmp` 没有；源码注释自认 “should not be done in case of CE SN ... fix” | `bse/comenv.f:355-365`；`kick` 调用 `:148,248,259` | 代码 + 变体对照；CE 中形成 NS/BH 且仍束缚时 `ECC` 与 `SEPF` 不一致，`JORB` 用错 `ECC` | **D2** | 第十一轮 |

## 7. 条件 D2（3 项）

判定：D2 机制，但只在特定运行时配置、构建选项或启动方式下可达。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-BH-09` | 安装树不含 `precession` 包；`sys.path` 用相对 CWD 的 `./bhmerger/precession-master` | `inspiralZJ.cpp:65-74`；`Makefile:296-315` | `petar-install/{bin,include}` 无 `bhmerger/`；从 `/tmp` 启动 `import precession` → `ModuleNotFoundError`；安装版在数据目录运行 → 每次并合 `remnant_zj_success=false` → 回退 BSE `merge_`（质量/kick 全变），进程继续 | **条件 D2**（取决于启动目录；源码树内运行不触发） | 第九轮 |
| `ASTRA-BH-25` | 并合事件日志写积分帧 `pos/vel`，快照写 `pos+pcm.pos`（仅 `ext_mode!=off`）；`d8fc857` 又删掉 `center.pos/vel` 六列，文件内无列可反推 | `bhmerger_output.hpp:114-119`；`particle_base.hpp:385-390`；`Makefile:95-98` | 已核对 `capture()` 与 `printColumnWithOffset`；质量/自旋/轨道/相对反冲不受影响 | **条件 D2** | 第十二轮 |
| `ASTRA-FD-01` | `TreeForForce` 默认构造不初始化裸指针，析构却无条件 `delete[]` | `tree_for_force.hpp:383-395`；`tree_for_force_impl.hpp:140-170` | placement-new 脏内存探针：构造后不 `initialize()` 即析构 → ASan `SEGV @ tree_for_force.hpp:385`。`petar.hpp:2960-2966` 已用注释 *"avoid segmentation fault due to FDPS clear function bug"* 打补丁绕开 | **条件 D2**（构造后未完成 initialize 即析构；任何新增早退/异常路径复发） | 第八轮 |

## 8. D3（98 项）

判定：仅在**非默认宏**或**特定构建目标**下可达的 UB/编译失败/内存安全，或仅影响辅助接口、输入解析、运行可用性；存在“启用即崩溃/错值”消费路径的标为**潜在 D2**，启用后按规则 1 重新定级。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-AR-05` | `AR_CHANGEOVER` 分支引用未声明的 `invr` | `ar_interaction.hpp:306-312` | `-DAR_CHANGEOVER -fsyntax-only`：`‘invr’ was not declared; did you mean ‘inv_r’?` | **D3**（该分支从未编译通过） | N4 |
| `ASTRA-AR-06` | `AR_TIME_FUNCTION_MULTI_R` / `AR_SLOWDOWN_MASSRATIO` 分支不可编译 | `symplectic_integrator.h:570-600,1031` | `:585/:595/:596` 报 `_bini`、`gt_kick_inv` 未声明；`:1031` 对 `BinaryTree*` 用 `.bin` | **D3** | S6 |
| `ASTRA-AR-07` | `SymplecticStep` 写裸指针，读回前先清零 size 跳过 `clear()` | `symplectic_step.h:337-357` | 文件含进程地址；重复读泄漏旧 `cd_pair_`；`ar.cxx load_flag` 路径可达 | **D3**（序列化卫生 + 泄漏；无物理值错误） | S7 |
| `ASTRA-AR-08` | `AR_SLOWDOWN_ARRAY` 下合并后单存活粒子直接漂移缺少 `1/κ` | `symplectic_integrator.h:2268-2285` vs `:1093-1106` | TREE 下 `Δτ·vel` 与积分算子一致（默认不受影响）；ARRAY 下正常漂移带 `kappa_inv`，该分支不带，位置多走 κ 倍，随后速度仍按 `1/κ` 写回 | **D3**（仅 ARRAY 构型；默认 TREE 已排除，见[第七轮](astra_sdar_line_audit.md)） | 第七轮 |
| `ASTRA-AR-09` | 零质量成员追加为人工树根后轨道量/`r`/`stab` 未计算，复用上一棵树残值 | `AR/information.h:284-300`；`hermite_integrator.h:1707-1757`；`symplectic_integrator.h:1364-1385` | 探针：4 体树→2 正+2 零后根 `semi/r/ecc/ecca` 与旧树逐位相同；`stableCheckIter` 返回 1000；`checkBreak` 消费这些字段 | **D3（潜在 D2）**（需组内保留零质量成员；未做 BSE 端到端） | 第七轮 |
| `ASTRA-AR-10` | 时间同步首个子步用绝对时间 `_time_end/time_table[k]` 缩放步长 | `symplectic_integrator.h:2767-2775` vs `:2790` | 组形成时刻非 0 时正确因子应为 `(_time_end-time_0)/(time_table[k]-time_0)`；现因子偏小，同步步数膨胀，极端触发 `step_count_max` | **D3**（步长控制/可用性） | 第七轮 |
| `ASTRA-AR-11` | team 并行 AR 组循环的中断记录按临界区到达顺序选取 | `hermite_integrator.h:2631-2637,2689` vs 串行 `:2715`、task `:2543` | 多组同一步中断时返回哪个二元组不确定；当前 `hard.hpp:727-747` option 1 读后即 `clear()`，不可观测 | **D3 潜伏**（需 `PETAR_INNER_THREADS>1`+team） | 第七轮 |
| `ASTRA-AR-12` | `AR_SLOWDOWN_PERT_R4` 下 `calcPertFromForce` 返回 `G·(m_p m_pert/r⁴)`，与 `calcPertFromMR`/`calcPertFromBinary` 的 `r⁻⁴` 归一化不一致 | `ar_interaction.hpp:644-649` | SDAR 用 `calcPertFromBinary` 作 pert_in、`calcPertFromForce` 作 pert_out（`hermite_integrator.h:1819-1822`）→ 比值差 G≈3.95e13 | **D3**（非默认宏） | 第十二轮 |
| `ASTRA-H4-02` | `table_check_mask` 按活组数分配、按全局槽位索引 | `hermite_integrator.h:2842-2847` | 槽位 `{0,2}` + dt-changing 中断：ASan `dynamic-stack-buffer-overflow` 于 `:2845` | **D3**（生产构建不编译；**`Makefile:254` 的 `HARD_DEBFLAGS` 启用 `HERMITE_DEBUG`，`petar.hard.debug` 受影响**） | C1 |
| `ASTRA-H4-05` | `checkNewGroup` 的 `used_mask` 按 `n_particle+n_group` 定尺却用 `index_offset_group_` 索引 | `hermite_integrator.h:1933-1944` | 定尺依赖 `index_offset_group_==particles.getSize()`，仅由默认构建的空 `ASSERT` 保护；`reserveMem(N)` 后填 `n<N` 即写越界 | **D3 潜伏**（生产 `linkMemberArray` 满容量，不触发） | 第七轮 |
| `ASTRA-H4-06` | `ForceSoft::clear()` 与 `FPSoft::clearForce()` 不清 `acorr`，4 阶修正跨步累积 | `soft_ptcl.hpp:16-23,268-279`；`soft_force.hpp:74-76,118` | 探针（`-DKDKDK_4TH`）：`clear()` 后 `acorr=(4,5,6)`、`clearForce()` 后 `acorr=(9,9,9)`；SIMD 主核不重置，`GradientKick` 的 `CalcCorrect~` 是 `+=2*acorr` | **D3**（非默认宏；力路径数值错误） | 第十三轮 |
| `ASTRA-BSE-06` | MIX 分支不写 `menv/renv`；COEL 分支写入的是合并前包层值，C++ 无条件拷贝 | `bse_interface.h:1684,1786-1793`；`merge.f:45-62` | Valgrind definedness 四字段失败，origin 为 `BSEManager::merge` 栈；第六轮补充：CE 弹出后 `menv` 仍是合并前巨星值；“无物理消费者”只对 `menv/renv` 成立，同一次回写的 `r/mc/ospin` 见 `BSE-07` | **D3**（仅输出/诊断） | C3 |
| `ASTRA-BSE-08` | `evolveBinary` 双 `kw==15` 早退不写事件终止符，调用方读未初始化 `record[9][*]` | `bse_interface.h:1527-1532,541-545`；`ar_interaction.hpp:1386-1394` | `BinaryEvent.record[20][9]` 无初始化；早退路径不 `recordEvent`/不置 `-1`；对照单 `15` 分支会写终止符 | **D3（潜伏）**（需两成员同时 `kw=15`；常规 `postProcess` 会先 `destroy`，未构造出可达场景） | 第六轮 |
| `ASTRA-BH-03` | Python 调用返回 NULL 后直接 `PyArray_Check` | `inspiralZJ.cpp:397-425,477-488,510-529,553-573,593-605` | 注入 `remnantmass` 异常 → SIGSEGV，gdb 定位 `inspiralZJ.cpp:481` | **D3（潜在 D2）**（并合中崩溃即丢事件） | N5 |
| `ASTRA-BH-05` | `e==0` 时 `inspiral_precav` 的 `e` 返回 `(1,1)`，C++ 按 `[0,1]` 越界读 | `inspiralZJ.cpp:411-419`；`eccentricity.py:556-558,582-590` | 真实调用：`e=0` → shape `(1,1)`，`e=0.5/0.9` → `(1,2)`；越界读 8 字节，垃圾 `e` 可致 NaN 方向或错误方向并入回退 | **D3（潜在 D2）**（需 `_bin.ecc` 精确为 0） | 第六轮 |
| `ASTRA-BH-06` | `remnantkick` 返回非数组时静默零 kick 且不 abort；无形状/有限性检查 | `inspiralZJ.cpp:577-607,338-339`；`bse_interface.h:1817-1827` | 失败分支无 `CAL_REMNANT_RETURN`，`vkick` 保持 `{0,0,0}`；stderr 误写 `remnantspindirection`；`full_output=False` 时读 `data[1..3]` 越界 | **D3（潜在 D2）**（`full_output=True` 正常路径不受影响） | 第六轮 |
| `ASTRA-BH-08` | 束缚轨道 `e = sqrt(1+2EL²/(μG²m₁²m₂²))` 未裁剪，`e≲1e-8` 时舍入成负数 → `e=NaN` | `inspiralZJ.cpp:313-315,446-458` | 忠实复刻探针：20,000 条严格圆轨道 **8,065（40.3%）** radicand<0 → `e=NaN`；`e=1e-8` 2.6%、`e=1e-6` 0。`is_bound_inspiral` 变 false，**无诊断**走 direct 分支（注释假定 `e>=1`），remnant 方向/自旋/kick 改用遭遇时角度 | **D3**（潜在 D2：若生产并合中 `e≈0` 普遍） | 第九轮 |
| `ASTRA-BH-10` | 退化/近对齐自旋：`χ=0` 的 0/0、`θ=0/π` 的 `arccos(1±ulp)`、`deltachiroots` 近重根被 `wraproots` 置 NaN、随后 `np.random.uniform(0,τ=∞)` | `precession.py:943,1009,974,1040,1322,2981-2982,3624,3573` | `χ=0`：400/400 NaN（被 C++ 1e-6 扰动屏蔽）；`θ₁=θ₂=0`：`deltachiroots`→NaN→`OverflowError`；`θ₁=0`（θ₂ 任意）：246/2,000 `theta1=NaN`；κ 恰为 `kappalimits` 端点（文档的 `(kappa,chieff)` 初值）也 NaN；`tofdeltachi` 对 χ± 的 1-ulp 越界给 NaN → C++ `isfinite` 失败 → 静默用遭遇时角度 | **D3**（精确对齐/端点零测度；回退对对齐情形物理接近，对 `θ₁=0,θ₂` 任意则真实丢步） | 第九轮 |
| `ASTRA-BH-11` | `initializePython` 在宿主已初始化 Python 时早退，跳过 `_import_array()` 与 `sys.path` 设置 | `inspiralZJ.cpp:48` | 同进程嵌入（宿主已 `Py_Initialize`）→ `PyArray_API` 为空 → 首次并合 `PyArray_*` 崩溃；当前 AMUSE 接口为子进程，不触发 | **D3 潜伏** | 第九轮 |
| `ASTRA-BH-15` | `precomputedroots[:-1]` 把 3 行根解包给 3 个名字 | `precession.py:4376,4418` | 文档快路径传 3×N 根 → `ValueError: not enough values to unpack (expected 3, got 2)` | **D3**（库内部，无调用点） | 第九轮 |
| `ASTRA-BH-16` | `return np.stack[(which, kappa, chieff)]` 索引 ufunc 而非调用 | `precession.py:4745` | 任何 `widenutation_condition` 调用 → `TypeError: 'numpy._ArrayFunctionDispatcher' object is not subscriptable` | **D3**（库内部） | 第九轮 |
| `ASTRA-BH-17` | `eval_costheta12` 使用未定义 `r`（签名无此参数） | `precession.py:1045,1099` | 守恒量分支 → `NameError`；角度分支正常，`remnantspin` 只用角度分支 | **D3**（库内部） | 第九轮 |
| `ASTRA-BH-18` | `inspiral_hybrid` 把 `deltachi`/`r` 传给不接受它们的 `inspiral_precav`；残留 `print("uswitch, eswitch",…)` | `precession.py:5585`；`eccentricity.py:1320,1347` | 递减 r（正常吸积）→ `LinAlgError`；递增 r → 成功但 stdout 打印调试行 + LSODA 噪声 | **D3**（库内部，PeTar 不调用） | 第九轮 |
| `ASTRA-BH-20` | `rhs_precav` 的 `q==1` 分支丢弃第二项；近 q=1 分支本身病态 | `precession.py:4783-4786` | `u<=0`（r=∞ 起始）+ q==1：RHS 0.10125 vs q=0.999 的 0.0772（+31%），q=0.9999/0.99999 为 0.064/−0.219；C++ 用 `q=1−1e-10` 避开 | **D3**（有限起始 `u>0`） | 第九轮 |
| `ASTRA-BH-22` | `deltachilimits_definition` 用 `\|χ₁−qχ₂\|/(1+q)` 作 Δχ 上下限，与 `eval_deltachi` 的实际跨度 `±(χ₁+qχ₂)/(1+q)` 不符 | `precession.py:2756,3129-3130` | q=0.7、χ=(0.8,0.6)：定义给 ±0.2235，实际可达 ±0.7176（3.2×）；任何用 `deltachilimits(q,χ₁,χ₂)` 做合法性检查的调用都会拒绝约一半物理构型 | **D3**（库内部；`deltachilimits` 的几何分支） | 第九轮 |
| `ASTRA-BH-23` | `deltachitildeav2` 在 `m→0` 用 `tol=1e-7` 硬夹，但闭式在该区间灾难性相消 | `precession.py:3329-3331` | `m=1e-7`（及更小）返回 0.35527，精确极限 3/8=0.375（**低 5.3%**）；`m=1e-6` 仍正确；注释"m→0 极限精确到 1e-9"只对 `deltachitildeav` 成立。`eval_chip_rms` 在近对齐/共振处因此偏低 | **D3**（库内部；`eval_chip_rms` 路径） | 第九轮 |
| `ASTRA-BH-24` | `kapparesonances` 的 up-down/down-up 递归分支不可达；数组输入被 `map` 截断到最短维 | `precession.py:2558-2586,2594` | 近 `chieff=±(χ₁±qχ₂)/(1+q)` 的 20,000 组中 **3,458 组**抛 `ValueError: Input values are not compatible`（本应可算）；`kapparesonances([20,8],0.2,0.7,0.8,0.6)` 只返回 r=20 的 (2,1) 结果，其余元素被静默丢弃 | **D3**（库内部；`kappalimits`/`kapparescaling` 会套错二元组的限） | 第九轮 |
| `ASTRA-BH-26` | 并合日志只记检测时刻 `time_now`，未抓 `time_merger/merger_semi/merger_ecc`（随后被清零） | `bhmerger_output.hpp:159-181`；`ar_interaction.hpp:1745-1746` | `ParticleState` 无这 4 字段；偏晚量上界 `dt_drift` 但不可度量 | **D3** | 第十二轮 |
| `ASTRA-BH-27` | `write()` 忽略 `writeLocked` 的返回值；写循环中途失败会留下**半行**，事件被静默丢弃 | `bhmerger_output.hpp:50-66,179` | ENOSPC/权限错误时无任何提示；半行使整份列解析失效 | **D3** | 第十二轮 |
| `ASTRA-BH-28` | 表头只在 `!_append && _initialize_file` 时写，且 `_initialize_file` 仅 rank0 为真；rank0 表头写失败只把它自己的 `enabled_` 置 false | `bhmerger_output.hpp:137-147`；`petar.hpp:3056-3061` | 其余 rank 仍以 `enabled_=true` 写数据 → 有数据无表头；`append` 到不存在的文件同样无表头 | **D3** | 第十二轮 |
| `ASTRA-BH-29` | 表头/数据无单位与参考系标注；time/semi/pos 用代码单位，姊妹 `.bse` 日志同一事件乘 `tscale·year_to_day` 与 `rscale` | `bhmerger_output.hpp:139-145,171-177`；`ar_interaction.hpp:1688-1694` | 跨文件比对差 ~365（天 vs Myr）与 ~4.4e7（Rsun vs pc）；`tools/analysis/` 无 `.bhmerger` 读取器，错不会被工具拦住 | **D3** | 第十二轮 |
| `ASTRA-TR-02` | `getEpjFromId` 用混合 `tp_glb_[i]` 下标过滤压缩 `epj_sorted_` | `tree_for_force_impl.hpp:18-29`；`make_tree.hpp:475-555` | 合法重排注入：`[EP0,SP0,EP1]` 查 id22 返回 null，`[SP0,EP0,EP1]` 查 id11 返回 null | **D3**（需多 rank LET-SP；消费点 `cluster_list.hpp:483` 在 `SAVE_NEIGHBOR_ID_IN_FORCE_KERNEL` 内且**直接解引用**返回值 → 启用即空指针） | E1 |
| `ASTRA-TR-06` | `INTEGRATED_CUTOFF_FUNCTION` 下 `q_` 从不赋值，`readBinary` 允许半读 | `changeover.hpp:14-16,44-59,148-155,202-229` | Valgrind：`log` 依赖未初始化值，栈来源 `main`，读点 `changeover.hpp:226`；`rcount<1` 只要求读到 1 个 | **D3（潜伏）** | N3 |
| `ASTRA-TR-07` | 稳定性 `acos` 只夹上界，`ratio<-1` 得 NaN 后被判"稳定" | `stability.hpp:170,265-266`；`binary_tree.h:589` | `ratio=-1.0000000000000004` → `acos=nan` → `stable3check` 返回 1（=stable） | **D3（潜在 D1）**（极端反平行角动量） | N10 |
| `ASTRA-TR-08` | CM 预测位移写成 `dx²+2·dy·dz`，漏 `dz²` | `status.hpp:317` | 源码确认为 `dpos_pred[1]*dpos_pred[2]` 加两次；仅 `SMOOTH_CM_USING_VEL_PRED` 下编译（默认 `SMOOTH_CM_USING_RECORES`） | **D3（潜伏）**（已确认公式错误，动态未测） | E4 |
| `ASTRA-TR-10` | Kepler 方程 Newton 用 `sqrt(1-sin²)` 代替 `cos`，丢符号 | `binary_tree.h:253-277` | 探针：e=0.3/0.5/0.9 分别 9/18/32（共 39）个 M 不收敛（10⁵ 次后 abort）；调用点 `particle_distribution_generator.hpp:156` 的 e≪0.3 | **D3**（潜伏；`solveKepler` 无调用点） | 第七轮 |
| `ASTRA-TR-13` | 判据对 NaN 无防护（`mass=0`/`r=0`/`e²<0`）；NaN 使 `peri<r_crit` 恒 false，静默丢边 | `cluster_list.hpp:399-405` | `particleToSemiEcc` 的 `sqrt` 可为 NaN；`_index` 不收录 | **D3** | 第十二轮 |
| `ASTRA-TR-15` | `FPSoft(){}` 不初始化 `acc/pot/pos/vel`；转换构造器不补 `acorr` | `soft_ptcl.hpp:50,76-87,96-98` | 唯一默认构造路径是 `CORRECT_FORCE_DEBUG`（`petar.hpp:1245-1249,1266`） | **D3** | 第十二轮 |
| `ASTRA-TR-17` | `const PS::F64vec da = _pi.acc - _pi.acc;` 恒为 0；KDKDK4 的 `acorr` 丢掉全部 `da` 项 | `hard.hpp:1689,1750` | 应为 `_pi.acc - _pj.acc`；仅 `step_mode=kdkdk4` | **D3** | 第十二轮 |
| `ASTRA-TR-18` | `setDKDMode()` 把 `mode_` 置 1 后，`getDtStartContinue()` 返回 `coff_one_step_[0][1]`（drift 系数 `ds`），而调用方按 kick 使用；`next_is_kick_flag_` 也被翻反 | `kickdriftstep.hpp:27,56-73,150-158`；`petar.hpp:4282,4333` | 全树无 `setDKDMode` 调用点；启用 DKD 即静默互换 kick/drift | **D3**（潜伏，无调用者） | 第十二轮 |
| `ASTRA-TR-19` | `CalcForceEpSpQuadSimd` 用 j 全局下标 `i` 而不是槽位 `i_tmp` 调 `set_spj_one` | `soft_force.hpp:528-533`；对照 `phantomquad_for_p3t_x86.hpp:1305,1426` | 第 0 块正确，`ih>0` 时写错槽位、`run_spj` 读到旧值；需 `n_jp>NJMAX=131072` 才有第 2 块 → 当前不可达；`USE_QUAD` 亦非默认 | **D3 潜伏**（构建宏 + 输入规模；可达时按规则 1 升 D2） | 第十三轮 |
| `ASTRA-FD-02` | `MomentMonopole::set()` 对零质量 cell 做 `0/0`，NaN 沿 `accumulate` 爬到根 | `tree.hpp:564-570` | 零质量 EP 累积 → `pos=(-nan,-nan,-nan)`，父 cell 同样 NaN；但 `tree_walk.hpp:884` 的 `r_crit_sq < dis_sq` 对 NaN 为假 → 走树安全打开，`zeromass2.cpp` 在 K=4..64 零质量粒子下与直接求和**误差 0**；多 rank LET 判据用 `<=`（`tree_for_force_utils.hpp:721` 等）会把该 cell 静默不发送 | **D3（潜在多 rank D2）** | 第八轮 |
| `ASTRA-FD-03` | `REUSE_LIST` 复用上一轮排序排列与 `tc_*`，却按新粒子数 resize，新增下标未初始化 | `tree_for_force_impl.hpp:278-314, 383-431, 485-509`；`tree_for_force.hpp:846-861` | `adr_org_from_adr_sorted_*` 不重算 → `epj_org_[adr]` 越界读、力写错粒子 | **D3**（仅 `REUSE_LIST`；PeTar 用 `MAKE_LIST`） | 第八轮 |
| `ASTRA-FD-04` | θ=0 交互表分支只填 EP、不填 SP；且用 `tp_glb_` 位置当 `epj_sorted_` 下标 | `tree_for_force_impl.hpp:693-706, 710-735` | SP 相互作用被静默丢弃；多 rank 下下标错位 | **D3**（θ=0 且多 rank） | 第八轮 |
| `ASTRA-FD-05` | `calcForceDirectParallel` 镜像展开用 `epj[j*n_image+k]`，步长应为 `n_recv` | `tree_for_force_impl.hpp:911-917`；`particle_system.hpp:1390-1397` | `n_image>n_recv` 时堆越界写；大小关系相反时镜像错乱 | **D3**（直接力 API + 周期边界；全树无调用点） | 第八轮 |
| `ASTRA-FD-06` | `-DLOOP_TREE` 下 per-thread 数组维度硬编码 12，初始化上界却是 `omp_get_max_threads()` | `tree_for_force_impl_force.hpp:1421-1434`；`ps_defs.hpp:2667` | >12 线程在并行区外越界写栈 | **D3**（`-DLOOP_TREE` 未定义） | 第八轮 |
| `ASTRA-FD-07` | 以"缓存的 `omp_get_max_threads()`"当实际 team 大小做切分/索引 | `ps_defs.hpp:2334-2354`；`particle_system.hpp:879-967`；`utils2.hpp:164-171`；`sort.hpp:445-453` | `OMP_DYNAMIC=true`/`OMP_THREAD_LIMIT` 下缺失分片不处理 → 未初始化下标打包、树序错、粒子静默丢失 | **D3**（非默认 OpenMP 运行时） | 第八轮 |
| `ASTRA-FD-09` | `collectSampleParticle` 首次分支 `freeMem()` 把 `size_` 归零，二次采样写回旧样本 | `domain_info.hpp:215-257` | 第二个 `ParticleSystem` 的样本被丢弃/读陈旧内存（`petar.hpp:2074-2079` 的 DM 双系统路径） | **D3**（`DARKMATTER` 构建） | 第八轮 |
| `ASTRA-FD-10` | `reserve()` 先更新 `capacity_` 再 `ReallocInner`，后者按新容量 memcpy 旧段 | `reallocatable_array.hpp:339`（配合 `:207-208, 236-245`） | 从旧 `data_` 越界读 `(new-old)*sizeof(T)` 字节；池内分配使 ASan 不可见 | **D3**（现有调用点都安全，契约缺陷） | 第八轮 |
| `ASTRA-FD-11` | `MemoryPool::clearAll()` 置 `top_=NULL` 但保留 `cap_`，之后 `alloc` 直接返回 NULL | `memory_pool.hpp:330-351` | `flag_realloc` 恒假 → `ret=top_prev=NULL` 且记账成功 | **D3**（`clearAll` 无调用点） | 第八轮 |
| `ASTRA-FD-12` | `CalcRank1D` 的 `while` 对 `domain[]` 无界搜索 | `util.hpp:598-603`；`particle_system.hpp:872, 911` | 粒子恰好落在 root 域上界面（`contained` 为 high-exclusive）或越界时读 `domain[n_proc]`，`srank==n_proc` 可致堆越界写 | **D3**（root 域为 ±2.1e37 时不可达；防御性缺口） | 第八轮 |
| `ASTRA-FD-13` | `LinkCell` 在 `n_cnt==0` 时用 `adr[-1]` 作 cell 下标 | `tree_for_force_utils.hpp:335-336`（同类 `:478`） | 仅当 Morton 序在父区间被破坏时出现；`LinkCell` 有 assert（abort），`LinkCellST` 无 | **D3**（潜伏；不变式成立时不可达） | 第八轮 |
| `ASTRA-FD-14` | 旧 `MakeInteractionList*`/`MakeIPGroup*` 变体混用两套下标空间（SP 未 `ClearMSB`、EP 用 tp 下标） | `tree_for_force_utils.hpp:1268-1276, 1482-1490, 1556-1563, 547-548` | `sp_first[SetMSB(i)]` 负下标；EP 取到错粒子。全仓无调用点、未实例化 | **D3**（非默认；与 `TR-02` 同族但不同行） | 第八轮 |
| `ASTRA-FD-15` | `FindScatterParticleLongCutoffP2P` 函数体被掏空，所有出参原样返回 | `tree_for_force_utils2.hpp:331-583` | 存活段（349-379）零 `push_back`；启用后全局树只含本地粒子，远程引力静默丢失 | **D3（启用即 D1）**（`LONG_CUTOFF`+P2P+多 rank） | 第八轮 |
| `ASTRA-FD-21` | `Matrix2::operator*` 的 xy 分量写成 `xx + rhs.xy`，应为 `xx * rhs.xy` | `matrix2.hpp:59-61`（另 `:95-97` 的 second invariant 连带） | A=[[1,2],[3,4]]、B=I 时 `(A*B).xy=3`，应为 2；`Matrix2` 全仓仅 typedef、未实例化 | **D3**（non-default） | 第八轮 |
| `ASTRA-FD-22` | `MortonKey::getKey` 把负坐标 `(U64)` 转换（UB），越界点静默折叠到 key 0/nmax；`half_len==0` 时 `0*inf=NaN` | `key.hpp:161-178` | `getKey((-1e12,0,0))` 与 `getKey((-1,0,0))` 在 box=[-1,1]³ 下逐级 `getCellID` 相同；当前 root cell 由 `pos±r_search*1.000001` 构造，不变量成立 | **D3**（潜伏；一旦先更新位置再建树即静默错树） | 第八轮 |
| `ASTRA-FD-23` | 矩检查函数用 `tp` 位置直取 `ep[i]`/`sp[i]`，未按 `adr_ptcl_`/`ClearMSB` 映射 | `tree_for_force_utils_moment.hpp:574-583` | 混合全局树下参考 COM 错、`i>=n_sp` 时越界读 `sp`；仅 debug check，`check_impl` 未 include | **D3**（non-default） | 第八轮 |
| `ASTRA-FD-24` | LET-P2P 的 `pos_domain_allgather_[adr_sp_tmp]` 只有下界检查、无上界 | `tree_for_force_utils_moment.hpp:331-343` | 335 行留有注释掉的 assert；块次序假设不成立即越界读并污染 `geo_.setSize()` | **D3**（需 `EXCHANGE_LET_P2P_FAST`，默认 A2A） | 第八轮 |
| `ASTRA-IO-04` | `format_transfer` 的 `FileHeaderNoOffset` 用 `sizeof(FileHeader)` 读写 | `format_transfer.cxx:45-62,74-80,84-86`；`io.hpp:210-268` | `ext_mode=off`：24/24 正常；`RECORD_CM_IN_HEADER`：72 vs 24；ASan `WRITE of size 72` 于 `format_transfer.cxx:75` | **D3**（`ext_mode!=off` + `-b -c`；当前生产 `ext_mode=off` 不触发） | 第六轮 |
| `ASTRA-IO-05` | 组积分器 `writeBinary` 漏写 `gt_kick_inv_`、sd 能量与 `binary_slowdown` | `symplectic_integrator.h:3505-3521,3526-3564` | 读回后 TTL `getH()` 除零、sd 账本全 0、ARRAY `binary_slowdown[0]` 越界；当前 `src/` 无调用点 | **D3 潜伏**（与 AR-07 机制不同） | 第七轮 |
| `ASTRA-IO-06` | 二进制快照 `fwrite(&acc, 6)` 依赖字段布局；`KDKDK_4TH` 在 `acc` 后插入 `acorr[3]` → 写出 acc+acorr，`pot_tot/pot_soft/n_ngb` 从未落盘 | `soft_ptcl.hpp:122-129,157-174` | 字段顺序与 fwrite 计数对照；ascii 路径仍写 5+1 项，两格式不一致 | **D3**（启用即错值） | 第十二轮 |
| `ASTRA-IO-07` | 快照头与粒子均为裸 `sizeof` 内存转储，**无 magic/version/构建签名**；`FileHeader` 在 `RECORD_CM_IN_HEADER` 下多 6 个 double，粒子结构随 `BHMERGER`/`BSE` 宏改变 | `io.hpp:210-285`；`particle_base.hpp:243-245`；`ptcl.hpp:158-163` | 读端只检查 `rcount<1`；用不同宏集构建读旧快照时头/字段静默错位（时间、质量、位置全变垃圾），不 abort | **D3（潜在 D2）** | 第十二轮 |
| `ASTRA-EN-03` | 删除动能只减 `0.5 m v_rel²`，未加回 `pcm.vel` | `petar.hpp:2140-2141,2162-2163,2576,2603`；`energy.hpp:234-235` | `RECORD_CM_IN_HEADER` 探针：同一物理状态 `V_cm=0` 误差 0、`V_cm=3` 误差 −7.5 | **D3（潜伏）**；`ext_mode!=off` 时 **D2**。**当前 `Makefile:15` 为 `ext_mode=off`，第 97 行只在 `ifneq($(ext_mode),off)` 定义该宏 → 缺陷代码未被编译** | B9 |
| `ASTRA-EN-04` | `writeBinary` 从 `&ekin` 起写整个对象 | `energy.hpp:202-203` | ASan `stack-buffer-overflow`；`src/` 内未找到主程序调用点 | **D3** | C4 |
| `ASTRA-EN-06` | `finishInterruptDrift` 用局部 `mass_modify_list` 遮蔽成员，中断收尾的 `dm` 不进能量账本也不清零 | `petar.hpp:1911-1913`；对照 `:1827,2003,2536-2555` | `correctSoftPotMassChange` 消费成员表并清 `dm`；局部表随作用域销毁。单 rank 该段不编译 | **D3**（MPI 多 rank；若 MPI 计入受支持配置则为 D2 账本项） | 第十三轮 |
| `ASTRA-BL-01` | `--id-offset` 用 32-bit `atoi` 解析 `PS::S64` | `petar.hpp:455` | `--id-offset 3000000000 -h` 输出 `−1294967296` | **D3**（`-h` 即退出，未证明静默 ID 碰撞） | C2 |
| `ASTRA-BL-02` | `COMM::List::operator=` 调用不存在的 `reserveMem(nmax,mode)` | `list.h:48,121,127`；`particle_group.h:37-43` | `List<int> a,b; a=b;` 编译报 `no matching function for call to reserveMem(const int&, ListMode&)` | **D3**（任何 List/组赋值实例化即失败） | S4 |
| `ASTRA-BL-03` | `--disable-print-info` 早于 `-p` 时整段参数文件读取被跳过 | `petar.hpp:585-599,415-418`；`bse_interface.h:1042-1060` | 三次对照：`--disable-print-info -p missing.par` 不报 `.par`，只报 `.par.bse`；反序则报 `.par` | **D3**（输入解析；需要非推荐顺序） | 第六轮 |
| `ASTRA-BL-04` | `checkParams` 的 `data_format` 条件用逻辑或连接、恒真；NDEBUG 下 `-i 4` 静默不写快照 | `petar.hpp:673,511-516,2250-2255,3155` | 非 NDEBUG 构建由 `:514` 的 `&&` assert 拦截（实测 abort）；`output()` 无 4 分支但 `nfile++` 已执行 | **D3（潜伏）**（需 `-DNDEBUG` 或程序化赋值） | 第六轮 |
| `ASTRA-BL-05` | `ForceSoft::grav_const` 缺省 `1.0` 与 par 缺省同值；AU 数据若既不写 G 也不带 `-u 1`，引力静默弱 39.5× | `static_variables.hpp:11`；`petar.hpp:187,3457` | `-u 1` 会覆写 G（已核对）；故属输入解析/毒默认值 | **D3** | 第十二轮 |
| `ASTRA-BL-06` | `HERMITE_PERT_FORCE` 分支调用的 `calcAccJerkPerturber` 在 `HermiteInteraction` 中不存在 | `hermite_integrator.h:769-772,835-838` | 全树无该成员定义；该宏从未编译过 | **D3** | 第十二轮 |
| `ASTRA-BL-07` | `static HardDumpList hard_dump;` 写在头文件里（内部链接），每个包含它的 TU 各有一份；唯一 `initial()` 在 `petar.hpp:3084` | `hard_assert.hpp:222`；`petar.hpp:3084` | 当前单 TU 构建（`Makefile:238` 只编 `main.cc`）无影响；将来新增一个使用 `DATADUMP` 的 TU 会拿到 `size=0` 的列表 → assert/空指针 | **D3**（构建结构依赖） | 第十二轮 |
| `ASTRA-BL-08` | `checkParams` 的软化上界 assert 写成 `eps <= ratio_r_cut`（无量纲），意图是 `eps <= r_in = ratio_r_cut*r_out` | `petar.hpp:667`；注释同处；自势修正 `hard.hpp:1802,3368,3469` | `r_out=0.01 pc` 时 `eps=0.05` 仍通过检查，但树里的自势是 `Gm/eps` 而修正用 `Gm/r_out` → 残留假势 `Gm(1/r_out-1/eps)`；默认 `eps=0` 不触发 | **D3**（配置错误被静默放行） | 第十二轮 |
| `ASTRA-BL-11` | `static_variables.hpp` 头内定义 `Ptcl::`/`EPISoft::`/`ForceSoft::` 静态成员；`r_out=0`、`search_factor=0` 毒默认 | `static_variables.hpp:3-11`；包含图仅 `petar.hpp:86` ← `main.cc` | 单 TU 链接通过；第二 TU 包含即重复定义；`r_out=0` 在 `initialParameters` 前被读则搜索半径为 0 | **D3 潜伏**（ODR）；默认值 D4 | 第十三轮 |
| `ASTRA-BL-12` | `regularTimeStep` 对 `_dt=0`/`+inf` 死循环；`dt_soft` 的零/无穷回退无保护 | `petar.hpp:850-859,3573,3586,3598-3612`；`checkParams:675` 只断言 `>=0` | 零速度散度 → `0.1*r_out/0=inf` → `while(dt<=inf)` 恒真；`n_glb==1` 且显式 `-t 0`（默认 `time_end=10`）→ `dt_soft=0` → `setStep(0)` 全零系数、时间不前进；均无报错 | **D3**（输入退化 → 静默挂起） | 第十三轮 |
| `ASTRA-DG-07` | `dump_number++` 在并行路径非原子；`hard_large_energy` 的 DATADUMP 位于并行 cluster 循环 | `hard_assert.hpp:196,209`；`hard.hpp:1260,2543` | 竞争 + 文件名重复覆盖；错误路径可达 | **D3** | 第十二轮 |
| `ASTRA-DG-12` | `readOneCluster` 先 `resizeNoInitialize(n_group)` 再检查 `n_group<0` | `hard_assert.hpp:139-152` | 负值先转成巨大 size_t → `bad_alloc`，负值检查成死代码；仅损坏/手工 dump | **D3** | 第十二轮 |
| `ASTRA-DG-14` | `ASSERT` 宏展开为裸 `if(!(expr)){...}`，无 `do{}while(0)`、无尾分号 | `hard_assert.hpp:231-236` | `if(c) ASSERT(x); else` 会绑错 else；仅 `HARD_DEBUG`（默认关） | **D3**（潜伏） | 第十二轮 |
| `ASTRA-GEN-01` | Kepler 盘三个分布函数同用 `init_genrand(getRank())`，`a` 与 `e/i` 取自同一 MT 流的相邻位置；`seed` 只改角度/相位 | `particle_distribution_generator.hpp:15,36,61,130-132` | 三处种子逐字相同；`main.cc:31` 已注释 `generateKeplerDisk()` → 潜伏 | **D3**（启用即 D2） | 第十二轮 |
| `ASTRA-GEN-03` | `assert(n_proc%4==0)` + 同名 `if(n_proc%4==0){…}` 包住整个函数体 | `domain.hpp:29-152` | NDEBUG 下 assert 消失、函数体被静默跳过；`DomainDecision` 无调用点 → 不可达 | **D3**（潜伏） | 第十二轮 |
| `ASTRA-GEN-04` | Plummer 归一化 `coef=1/sqrt(r_scale)` 隐含 G=M=1，缺 `sqrt(G·mass_glb)` | `particle_distribution_generator.hpp:237-242` | `mass_glb≠1` 时速度小 √M 倍 → Q=0.5/M；调用点硬编码 `m_tot=1.0`，当前不触发 | **D3**（潜伏） | 第十二轮 |
| `ASTRA-GEN-05` | Kepler 盘生成器默认 `gravitational_constant=1.0`（AU/Msun/5e6s 体系）而 peTar 用 pc/Myr；且不生成中心恒星 | `particle_distribution_generator.hpp:127`；`petar.hpp:3387` | `main.cc:31` 已注释 `generateKeplerDisk()` → 不可达；启用即坐标/速度/力场不自洽 | **D3**（潜伏，启用即 D1） | 第十二轮 |
| `ASTRA-MC-07` | `-f 2` + 单个 `-m`：质量写成 `1/N`，`Rh3D` 累计不到 0.5M → 除零 | `main.c:595-601`、`962` | `-N 100 -f 2 -m 1` → exit 136，gdb 停在 `main.c:962`；`--dm-nfw` 下不崩溃但 `Mstar=1.0`（应 100） | **D3** | 第十轮 |
| `ASTRA-MC-10` | `-f 3` / `-f 4` 被解析器拒绝（只收 0..2），实现与帮助都声称支持 | `main.c:317-325`；`1571/1652`；`5383-5385` | `-f 3` → `Error: Mass function (f) needs to be between 0 and 2` | **D3** | 第十轮 |
| `ASTRA-MC-11` | DM 输出缺 `--with-interrupt=sse` 的 24 列布局 | `darkmatter.h:317-334` | 用 PeTar 真实 `readAscii` 交叉读：`(sse, 20 列)` 第 1 行把 id 读进 `dm`、第 3 行 abort；(off, 34 列) 同样 abort(134) | **D3** | 第十轮 |
| `ASTRA-MC-12` | DM 输出在质心系而非密度中心；`\|cm\|=1.7–2.2 pc > Rh` | `darkmatter.h:308-315,356-357` | 以原点量 `Rh=1.96/2.40/2.11/2.32`，加回 `subtracted_common_position` 后 1.012/0.942/1.012/1.013 | **D3** | 第十轮 |
| `ASTRA-MC-13` | `epochstar` 恒 0，与 `-e` 演化后的质量/`kw` 矛盾 | `main.c:647,650` | 25.93 M☉ @5 Myr 的 `epoch1=-0.0869`；`petar.init` 一律按 t=0 起算 → 二次演化 | **D3** | 第十轮 |
| `ASTRA-MC-14` | 缺少 `CALL instar`，`/TYPES/ ktype` 全 0 → 并合产物 `kw=0` | `evolv2.f:1193`、`mix.f:35`、`comenv.f:113` | 子审；仅 `-B` 且发生并合 | **D3** | 第十轮 |
| `ASTRA-MC-15` | >100 M☉ 的 100 M☉ 夹紧被启用 | `star.f:56`；`hrdiag.f:76/78/780-785` | 150 M☉ ZAMS：`R=17.226 vs 22.861`（−24.6%）、`L=1.271e6 vs 2.480e6`（−48.8%） | **D3** | 第十轮 |
| `ASTRA-MC-16` | GCC 15 下 Makefile 构建失败（`feenableexcept` 缺 `_GNU_SOURCE`；`output5` 无原型） | `main.c:166,1244,5201` | `make mcluster` → 2 error；可用 `-D_GNU_SOURCE -std=gnu17 -Wno-error=implicit-function-declaration` 绕过 | **D3** | 第十轮 |
| `ASTRA-MC-17` | King/EFF 二分查找在端点可越界 | `main.c:2344,3018,3040` | 静态；当前采样不触发 | **D3** | 第十轮 |
| `ASTRA-BSI-04` | 默认 `bse/instar.f` 缺 `ktype(11,11)=12`；`bseEmp`/`mobse` 都有 | `bse/instar.f:31-36`；`bseEmp/instar.f:37` | 变体 `diff`；CO+CO 双 WD 并合在默认 BSE 被当 CO WD（总量超 `Mch` 时 `kw=15` 无残留），不是 ONe | **D3** | 第十一轮 |
| `ASTRA-BSI-05` | BSE 解析器复用 PeTar `-p`，无条件打开 `<prefix>.bse`，缺失即 `abort()`；`stellar-evolution=0` 时也执行 | `bse_interface.h:1042-1060`；`petar.hpp:2953` | `./petar.bse -p /tmp/nonexistent 1.0` 打印 “Cannot open file ... .bse” 后终止；`petar -p input.par` 无 `input.par.bse` 即启动失败 | **D3** | 第十一轮 |
| `ASTRA-BSI-06` | `BSEEMP` 帮助允许 `--bse-metallicity 0`，但 `checkParams()` 断言 `z>0`，`zcnsts` 对 `log10(0)` 无保护 | `bse_interface.h:1139,752`；`bseEmp/zcnsts.f:50` | debug 断言失败；NDEBUG 下 `nint(-inf)`、z 被改成 1.22e-4 | **D3** | 第十一轮 |
| `ASTRA-BSI-07` | `bse_test` 错误路径备份写成 `bin[1].star[1]`，应为 `bin[i]`；`nbin==1` 时越界读 | `bse_test.cxx:400` | `-D_GLIBCXX_ASSERTIONS` 单双星运行在 `stl_vector.h:1263` 断言 `__n < size()` | **D3**（测试工具） | 第十一轮 |
| `ASTRA-BSI-08` | `bse_test` 在 `omp parallel for` 内对共享 `std::vector` 做 `push_back`，并并发调 Fortran BSE 无 `critical` 锁 | `bse_test.cxx:374,494-508` | 代码；vector 重分配竞争 + `RAND3`/COMMON 竞争；生产 `ar_interaction.hpp` 有 `critical(petar_bse_fortran)` | **D3**（测试工具） | 第十一轮 |
| `ASTRA-BSI-12` | `dgcore` 的 `temp` 计算后未用：`if(temp.gt.the)`/`if(temp.gt.tc)` 被注释，`x=1`、`y=0` 写死；He→C 无条件、C/O→Ne 永不发生 | `bse/dgcore.f:77-86` | 代码；三变体一致；`kw3=12` 通道不可达，`de`/`ebinde` 偏高 | **D3**（模型） | 第十一轮 |
| `ASTRA-BSI-13` | `FBFAC=FBTOT/(mt-mcx)` 无上限；`nsflag=3` 薄包层 `mt-mcx<0.2` 时 `FBFAC>1`，`kick` 的 `vk*(1-FBFAC)` 反向放大 | `bse/hrdiag.f:647-648`；`bse/kick.f:155` | 代码；`mt==mcx` 除零；触发率依赖是否存在 `mt-mcx<0.2` 的 CCSN 前身 | **D3** | 第十一轮 |

## 9. D4（54 项）

判定：已核对仅影响计时、分类统计或打印，物理值与账本不受该缺陷影响。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-AR-04` | destroy 分支漏记 slowdown 事件分类项 | `symplectic_integrator.h:2130-2136` vs `:2225-2228` | 真实两体回归：分类项 `0→0.702481…`，总修正不变 | **D4**（总参考仍正确） | D2 |
| `ASTRA-AR-13` | `findCloseSoftPert` 在 `_n_tt==0` 时读 `_tt[-1].group_id` | `ar_perturber.hpp:46-76`；调用点 `hard.hpp:637,812` 均在 `/* */` 块内 | 无活调用者；`SOFT_PERT` 默认定义但 `soft_pert` 由 `hard.hpp:430/547/796` 直接赋值 | D4（死代码） | 第十三轮 |
| `ASTRA-H4-04` | `ParticleH4::print` 把 `acc1` 打成 `acc0`，数组按地址输出 | `hermite_particle.h:104-109` | 探针 `acc0=0x… acc1=0x…` 两地址相同 | **D4** | S5 |
| `ASTRA-H4-07` | `HermiteInformation::time_origin` 无初值、`clear()` 不置位 | `hermite_information.hpp:8,36`；赋值点 `hard.hpp:658,848`；唯一读取 `printColumn` | 默认构建不打印（`HARD_DEBUG_PRINT` 未定义）；启用打印且在 `initialIntegration` 前调用即读未初始化值 | D4 | 第十三轮 |
| `ASTRA-BH-04` | 每次 BH 并合泄漏 kwargs float、空 tuple 与一个 array 引用 | `inspiralZJ.cpp:368-386,374-383` | 20000 次调用 RSS `27736→38408 KB`，约 `546 B/次`；40000 次线性增长 | **D4** | N6 |
| `ASTRA-BH-07` | `PyArray_SimpleNewFromData` 包装栈数组并额外 `INCREF`，每次 inspiral 泄漏且留下悬垂指针 | `inspiralZJ.cpp:368-371,377,386`；`eccentricity.py:571` | 净引用 1；`a_arr_data` 为栈数组；`outcome['a']` 是同一对象 | **D4**（当前读取路径不触发 UAF；与 BH-04 同类但为独立站点） | 第六轮 |
| `ASTRA-BH-12` | `odeint` 容差守卫拼写为 `'rol'`/`'aol'`，用户 `rtol/atol` 被静默覆盖 | `precession.py:4852-4854,5265-5267` | 传 `rtol=1e-13` 实际仍用 1e-10；默认调用不变 | **D4**（潜伏） | 第九轮 |
| `ASTRA-BH-13` | `rhs_precav` 返回 `float(1-element ndarray)`，NumPy≥2 弃用 | `precession.py:4808` | 一次 inspiral 43 条 `Conversion of an array with ndim > 0 to a scalar is deprecated`；NumPy 移除该转换后变 `TypeError` → 生产积分器停止 | **D4 现在 / 前瞻 D2** | 第九轮 |
| `ASTRA-BH-14` | `np.in1d`（NumPy≥2 弃用） | `precession.py:5018,5463` | 每次 inspiral 1 条 `in1d is deprecated` warning；移除后 `TypeError`；上游为 `np.isin` | **D4** | 第九轮 |
| `ASTRA-BH-19` | `remnantspindirection` 文档说返回 `thetaL` 实返 `Ĵ`；`remnantkick` 文档默认 `kms=False/full_output=False` 实为 True；`deltachiroots`/`kappalimits`/`deltachilimits` 的 Returns 顺序与代码相反 | `precession.py:5950,5987,5991,6021-6022,6031`；`2416-2419,2660-2663,2967-2972,3115-3120` | 按文档调用得 `(N,4)` 与 km/s；把 χ₃ 当 χ₋ 会静默错位；C++ 恰好依赖 `True/True`，把默认值"改回文档"会破坏接口 | **D4**（契约，只改文档） | 第九轮 |
| `ASTRA-BH-21` | `eval_S`/`eval_kappa` 的 `TypeError` 未 `raise`（第二分支写成 `if`）；`eval_q` 断言 `q<1` | `precession.py:538,1440-1441,1479,1489-1490` | `eval_S()`→`UnboundLocalError`；`eval_q(2,2)`→`AssertionError` | **D4**（库内部） | 第九轮 |
| `ASTRA-BH-30` | `capture()` 不记录 `binary_state`（pair id/中断态） | `bhmerger_output.hpp:110-124` | 事件无法与 `setBinaryPairID`/中断类型对齐；两侧粒子 id 已在，合并树仍可重建 | **D4** | 第十二轮 |
| `ASTRA-TR-14` | 覆写 `sys[i].n_ngb` 为「判据通过数、不含自身」，与 `:454` 单粒子分支保留的 1 语义漂移；只影响输出列 `N_b` | `cluster_list.hpp:493,525` | 两分支对「无邻居」编码不同（1 vs 0）；`driftSingle` 死代码不消费 | **D4** | 第十二轮 |
| `ASTRA-TR-16` | `group_data` 是 `artificial` 与 `cm` 的 16 B 完全重叠 union；`EPISoft::copyFromFP` 无条件按 `.artificial` 判 `type`，而建树时 mode 可能仍是 `cm` | `ptcl.hpp:20-38`；`soft_ptcl.hpp:298` | 默认路径确有该窗口（`initial()`→`setParticleGroupDataToCMData`→`treeNeighborSearch`）；当前 `type` 只被力核使用且力核在切回 `artificial` 之后 | **D4（条件 D1 潜伏）** | 第十二轮 |
| `ASTRA-TR-20` | `KickDriftStep` one-step API 无调用者，`%=coff_one_step_.size()` 在 `setStep` 前除零；`coff_[0]` 空表下标 | `kickdriftstep.hpp:28-31,112-119` | 生产只用 continue API（`petar.hpp:3891,3922,3943,4053-4417`），其顺序本轮核对自洽 | D4（潜伏 UB） | 第十三轮 |
| `ASTRA-FD-08` | `RandomSampling` 的 `j` 抽自 `[0,n_src)` 而非 `[i,n_src)`，子集分布偏斜 | `ps_defs.hpp:2377-2382`（同型 2405-2410、2439-2444） | 4×10⁵ 次探针：`{0,1}` 0.2496（应 0.1667）、`{2,3}` 0.1250；只影响域分解样本/排序枢轴质量 | **D4**（默认可达，采样质量） | 第八轮 |
| `ASTRA-FD-16` | `~TreeForForce` 漏删 `adr_epj/adr_spj/adr_ipg_for_force_` 与 MPI `status_` | `tree_for_force.hpp:384-395`；`tree_for_force_impl.hpp:157-159, 170` | 串行 ASan：7 次构造 → 3 组 `Direct leak` 共 21 个分配、2856 B | **D4**（资源泄漏） | 第八轮 |
| `ASTRA-FD-17` | `getMemSizeUsed()` 循环体用 `tmp =` 而非 `tmp +=`，且漏计多组数组 | `tree_for_force_impl.hpp:43-48` | 只统计最后一个线程的 5 个数组，8 线程时低报约 8× | **D4** | 第八轮 |
| `ASTRA-FD-18` | `ParticleSystem::n_proc_exch_` 未初始化即 `+=` | `particle_system.hpp:23, 146-151, 1064-1066` | 唯一置零处 `clearCounterAll()` 只在 `profile.hpp:168` 的注释块内被调用 | **D4**（只污染交换计数） | 第八轮 |
| `ASTRA-FD-19` | `TimeProfile` 五个成员构造/`clear()` 前未初始化，却每步参与 `operator+` | `ps_defs.hpp:2035-2039, 2117-2138, 2167-2171` | `walK_tree/keep_list/copy_ep/dispatch/retrieve` 不在 dump/total 中，故无可见输出 | **D4** | 第八轮 |
| `ASTRA-FD-20` | 周期边界闸门只判 `SEARCH_MODE_LONG`，漏掉 `LONG_SYMMETRY`/`LONG_SCATTER` | `tree_for_force_impl_exlet.hpp:96-100`；`utils2.hpp:64-89` | `TagSearchBoundaryConditionOpenOnly` 丢弃 shift，周期轴下发送未平移镜像 → 静默错力 | **D4**（PeTar 不开周期边界） | 第八轮 |
| `ASTRA-IO-08` | 重启后 `initialStep` 无条件 `file_header.nfile--; output()`，使初始快照写回**重启输入文件本身**（data.N 被覆盖） | `petar.hpp:3925-3926,3994-3995,2235-2254` | 意图是避免同时间重复文件，但原快照被销毁；若写盘中途失败，重启文件损坏 | **D4** | 第十二轮 |
| `ASTRA-BL-09` | `checkParams` 对 `data_format` 的 assert 写成 `>=0 \|\| <=3`（恒真）；命令行路径 `case 'i'` 的 assert 才是正确的 `&&` | `petar.hpp:673`；对照 `:514`；`petar.hpp:3155-3158,2251-2254` | 从 par 文件读入非法值（如 7）时无校验：读取走 binary 分支，写出不匹配任何 `if` → **静默无快照输出**，运行继续 | **D4**（par 文件输入） | 第十二轮 |
| `ASTRA-BL-10` | `HardDumpList`：`dumpThread` 直接索引不校验；`initial()` 重复调用泄漏；`readOneCluster` 先按未校验 `n_group` resize 再查 `<0` | `hard_assert.hpp:173-177,204-214,139-152` | `HARD_DUMP` **默认开**（见下"宏更正"），但 dump 只在调试分支触发；正常路径 `ith<size` | D4 | 第十三轮 |
| `ASTRA-BL-13` | `get_version.hpp` 是 `@VERSION@` 模板产物但 `Makefile` 无生成规则；banner 版本恒为写死的 `1268_298`；`printFeatures` 不打印 `HARD_DUMP`/`AR_DEBUG_DUMP` | `get_version.hpp:5`；`get_version.hpp.in`；`Makefile`（无 `get_version` 目标）；`petar.hpp:2820-2900` | `make` 构建的实际提交与 banner 版本无绑定；默认构建的 dump 宏在 banner 上不可见 | D4 | 第十三轮 |
| `ASTRA-DG-01` | Tree_NB 计时偏移与 `time_profile_` 累加竞争；AR 计步与 `noteArRegion` 裸累加 | `tree_for_force_impl_force.hpp:1284-1303`（第六轮补记 `:1297,1303`）；`hermite_integrator.h:2626-2627`；`inner_parallel.hpp:392-393` | 源码；`:1297`/`:1303` 的 `time_profile_.* += ...` 位于 `#pragma omp parallel for`（`:1289`）内，读改写竞争；`tree_profile.patch` 原型 32t 有收益迹象但不稳定 | **D4**（不参与力计算） | D1 |
| `ASTRA-DG-02` | `PtclHard::print` 的 `id_cluster`/`adr_org` 写到 `std::cerr` | `hard_ptcl.hpp:88-92` | 静态即可确认：目标流缺列、stderr 多出半行 | **D4** | N9 |
| `ASTRA-DG-03` | 并合路径无条件输出 `DEBUG Dynamic_merge invoked` 到 stderr 与 `fout_bse` | `ar_interaction.hpp:1672-1682` | 无 `#ifdef`、无 verbose 判断；`fout_bse` 中该行不符合列格式 | **D4**（会破坏按列解析） | N8 |
| `ASTRA-DG-04` | `-w 2` 的 `.status` 表头在读数据前写出，粒子列标题全缺 | `petar.hpp:3016-3018`；`main.cc:19-34`；`status.hpp:100` | 6 粒子运行：表头 38 列 vs 数据 260 列（缺 222 个粒子列名） | **D4**（输出格式；不影响物理） | 第六轮 |
| `ASTRA-DG-05` | `_n_ptcl>1000` 的“Large cluster”无条件写 stderr | `hard.hpp:320-333` | `#ifdef HARD_DEBUG` 已被注释，`std::cerr` 无任何条件；每个大组每次漂移一次 | **D4**（日志噪声） | 第六轮 |
| `ASTRA-DG-06` | `checkNeighborExist` 把平方距离 `r_min_sq` 与半径 `r_out` 直接比较 | `hard.hpp:924,937` | 量纲不一致；`n_neighbor_zero` 系统性偏低（只进统计） | **D4** | 第十二轮 |
| `ASTRA-DG-08` | `calcEnergyPertOneGroup` 不使用传入的 `_perturber`，实取 `_group.perturber.soft_pert`；NULL 时静默 0 | `hermite_interaction.hpp:448-460` | 力与势口径当前一致（同组扰动也为 0） | **D4** | 第十二轮 |
| `ASTRA-DG-09` | `Tprofile::end()` 连续两次 `PS::GetWtime()`，`tbar` 与 `time` 在不同时刻采样 | `profile.hpp:191-194` | 仅影响诊断计时（亚微秒级偏差） | **D4** | 第十二轮 |
| `ASTRA-DG-10` | `NumCounter* iptr = (NumCounter*)this+i` 跨成员指针算术（技术上 UB） | `profile.hpp:528-551` | 当前 16/19/14 个成员与类布局逐一相符、列头顺序一致，不会错列 | **D4** | 第十二轮 |
| `ASTRA-DG-11` | dump 只存 4 个 `Ptcl::` 静态量，不存 `EPISoft::eps/r_out` 与 `ForceSoft::grav_const` | `hard_assert.hpp:70-86,121-130` | 重放时若 par 的 eps/G 与 dump 不一致，软化尺度与截断半径不自洽；格式缺项，当前 hard_debug 读取同一 par 故无现状后果 | **D4** | 第十二轮 |
| `ASTRA-DG-13` | `readOneCluster` 覆写进程全局的 `Ptcl::search_factor/r_search_min/mean_mass_inv/r_group_crit_ratio` | `hard_assert.hpp:127-130` | 重放结束或再读第二个 dump 时静态量是陈旧的；单次重放是有意为之 | **D4** | 第十二轮 |
| `ASTRA-DG-15` | profile Min/Max 行取 `time−tbar`，普通行取 `time`；`getMax/getMin` 只回填 `time`，`tbar` 留 0 | `profile.hpp:208-214,403-421`；`petar.hpp:2442-2457,2498-2499` | 探针（真实 `profile.hpp`）：`time=10,tbar=4` → `getMin/getMax=6`、`dump=10`、`dumpBarrier=4`；min 的 barrier 列恒 0 | D4（诊断口径） | 第十三轮 |
| `ASTRA-DG-16` | 错误信息中的读取数与实际不一致 | `soft_ptcl.hpp:137`（读 5 报 6）；`hard_assert.hpp:124`（读 4 报 3） | 文本；诊断 | D4 | 第十三轮 |
| `ASTRA-GEN-02` | `DomainDecision` 全树无调用点（实际走 `dinfo.decomposeDomain()`） | `domain.hpp:26-152`；`petar.hpp:2079-2083` | `grep -rn DomainDecision` 仅命中定义；覆盖清单里“域分解已审”是假象 | **D4** | 第十二轮 |
| `ASTRA-MC-18` | 双星在 Q/Rh 归一化之后插入，此后不重标 | `main.c:1164` | `-b 0.5`：`Q=0.4913`、`Rh=0.8054`；`-p 0`：0.4993/0.7982。默认周期分布下 ≲2% | **D4** | 第十轮 |
| `ASTRA-MC-19` | `fractalize` 用陈旧 `i==N` 做速度归一化，实际未归一化任何星 | `main.c:3629-3634` | 静态；后续剖面都会重设速度大小 | **D4** | 第十轮 |
| `ASTRA-MC-20` | `strcat` 拼到 argv 字符串；`output4/5` 用 50 字节栈缓冲 | `main.c:4861-4869,5129,5204` | 长 `-o` 溢出/截断 | **D4** | 第十轮 |
| `ASTRA-MC-21` | "final energy check" 仍 `i<j-1`；`check` 无开关 | `main.c:1267` | N=100：`Epot=-4.35905` vs `-4.41478`（低 1.26%） | **D4** | 第十轮 |
| `ASTRA-MC-22` | `energy_order` 动能项带 `1/Nstars`，排序几乎只按 `-1/r` | `main.c:4430` | 静态；径向单调性不变 | **D4** | 第十轮 |
| `ASTRA-MC-23` | 潮汐半径的 `omega`（km/s/pc）当 1/Myr 用 | `main.c:769` | `rtide=12.9306` vs 正确 12.74（−1.5%） | **D4** | 第十轮 |
| `ASTRA-MC-24` | `--blackhole` 的 `R/L` 与 SSE 的 `kw=14` 约定不同 | `main.c:1057-1066` | `R=0.01 Rsun, L=0` vs SSE `R=4.24e-6·M, L=1e-10`；PeTar `-s bse` 会重算 | **D4** | 第十轮 |
| `ASTRA-MC-25` | DM 与恒星共用 `drand48` 流，`seed` 只写 `.dm.info` | `darkmatter.h:302-306` | 同 `-s 42`、`-N 256` vs `-N 512` 的首个 DM 粒子不同 | **D4** | 第十轮 |
| `ASTRA-MC-26` | DM 的 `step>0.065` 分支不可达（最大 0.0586） | `darkmatter.h:104-106` | 子审 | **D4** | 第十轮 |
| `ASTRA-MC-27` | 极端但合法 DM 输入触发 SIGFPE 而非报错 | `main.c:166`；`darkmatter.h:36-43` | `-R 1e-160 --dm-rs 1e-160 …` → exit 136，无输出 | **D4** | 第十轮 |
| `ASTRA-MC-28` | 自带测试的编译命令在 GCC 15 下失效 | `tests/test_darkmatter.py:37-41` | 改 `CC` 后全部通过（8.3 s） | **D4** | 第十轮 |
| `ASTRA-BSI-09` | `StarParameter::initial()` 不初始化 `lum`；`recordEvent()` 在首次 `evolveBinary` 前读它写初始事件 | `bse_interface.h:291-302,564-565` | 代码；初始事件记录 `lum1/lum2` 为栈上垃圾，仅诊断 | **D4** | 第十一轮 |
| `ASTRA-BSI-10` | `mlwind.f` 的 `/VALUE1/` 声明 4 个 double（`neta,bwind,hewind,mxns`），其余为 3 个；`hewind` 从未使用 | `bse/mlwind.f:10-11`；`bse/const_bse.h` | 代码；多出的 `mxns` 越界读，未使用故默认无数值影响；`--bse-hewind` 空操作 | **D4** | 第十一轮 |
| `ASTRA-BSI-11` | `zcnsts` 的 `msp(81)` 三次项写成 `lzs*lzs*xr(83)`，相邻多项式均为 `c0+lzs*(c1+lzs*c2)` | `bse/zcnsts.f:177` | 代码；Z=0.001 时 R_TMS(1.6–2 M☉) 偏 ~0.7%，太阳金属度无影响 | **D4** | 第十一轮 |
| `ASTRA-BSI-14` | `kick.f` 在 `DISP0=0`（`--bse-sigma 0`）时 `v(4)=0`，`v(k)*vk/v(4)` 产生 `0/0=NaN` | `bse/kick.f:178-181` | 代码；`vkick[3]=0` 使 `postProcess` 不写粒子速度，但污染诊断/后续 `merge` 初值 | **D4** | 第十一轮 |

## 10. 契约与非缺陷（2 项）

契约类缺口不强行套 D1–D4，写明在哪种契约下等价于 D1；“非缺陷”条目附代码事实，避免反复被重新报告。

| 编号 | 机制 | 文件/行 | 代码事实与证据 | 等级 | 来源/旧编号 |
|---|---|---|---|---|---|
| `ASTRA-IO-03` | 常规重启未串接 Fortran `/RAND3/`+`/VALUE3/` 与 NumPy RNG 状态；`inspiralZJ.cpp` 连 seed 都不设 | `ran3.f:10-12`；`bse_interface.h:1151-1180`；`ar_interaction.hpp:1370`；`hard_debug.cxx:181`；`precession.py:6081`；`eccentricity.py:569` | 显式 dump/read 可逐位恢复 8 个 float；仅重设同 seed 时第 2 次 kick z `+326.71` vs 连续 `−654.67`；dump 只在错误路径、read 只在 hard debug | **契约**（承诺逐事件连续重启时等价 D1；否则为可复现性限制） | E2 |
| `ASTRA-TR-09` | wrapper 只装载 `n_ip_local` 个目标却给 kernel 传 `n_ip` | `soft_force.hpp:330,404,470,535` | 三处 `set_*_one(n_ip_local)` 配 `run_*(n_ip)`；写回只到 `n_ip_local`，故不污染有效力 | **非缺陷**（纯浪费；data.70 过滤 0%、data.52 0.03%，A/B 无收益） | E3 |

## 11. 已撤回条目（不占编号）

| 旧编号 | 原主张 | 撤回理由 |
|---|---|---|
| N7 | 单星质量变化时的势能修正错误 | `hermite_interaction.hpp:76` 的 `gmor` 只含 `m_j`，`force_[k].pot` 已是单位质量势；`de_pot = pot*dm` 与解析式一致。组 CM 路径定义相同 |
| —（历史） | `StarParameterOut::vkick[4..7]` 未初始化 | `vkick` 已零初始化，`merge()` 先把两个输出的四分量各自复制到 `vkick[8]`；[PeTar_merger_bug_analysis.md](../02_baseline_and_defects/PeTar_merger_bug_analysis.md) 的该项在当前源码不成立 |
| —（B4 段内） | BSE 并合锁定修复会引入重入死锁 | capped/uncapped 两个接触案例均未重入死锁；**但"不能把整个 lambda 套同名锁"仍是有效设计约束**（`postProcess` 会调已加锁的 `modifyOneParticle`，OpenMP critical 不可重入） |

## 12. 旧编号对照

| 新编号 | 旧编号 | 新编号 | 旧编号 |
|---|---|---|---|
| `ASTRA-AR-01` | A1 | `ASTRA-TR-01` | A4 |
| `ASTRA-AR-02` | A5 | `ASTRA-TR-02` | E1 |
| `ASTRA-AR-03` | N1 | `ASTRA-TR-03` | N2 |
| `ASTRA-AR-04` | D2 | `ASTRA-TR-04` | S1 |
| `ASTRA-AR-05` | N4 | `ASTRA-TR-05` | B6 |
| `ASTRA-AR-06` | S6 | `ASTRA-TR-06` | N3 |
| `ASTRA-AR-07` | S7 | `ASTRA-TR-07` | N10 |
| `ASTRA-H4-01` | B3 | `ASTRA-TR-08` | E4 |
| `ASTRA-H4-02` | C1 | `ASTRA-TR-09` | E3 |
| `ASTRA-H4-03` | S3 | `ASTRA-IO-01` | A3 |
| `ASTRA-H4-04` | S5 | `ASTRA-IO-02` | S2 |
| `ASTRA-BSE-01` | B1 | `ASTRA-IO-03` | E2 |
| `ASTRA-BSE-02` | B2 | `ASTRA-EN-01` | B7 |
| `ASTRA-BSE-03` | B4 | `ASTRA-EN-02` | B8 |
| `ASTRA-BSE-04` | B5 | `ASTRA-EN-03` | B9 |
| `ASTRA-BSE-05` | B10 | `ASTRA-EN-04` | C4 |
| `ASTRA-BSE-06` | C3 | `ASTRA-BL-01` | C2 |
| `ASTRA-BH-01` | A6 | `ASTRA-BL-02` | S4 |
| `ASTRA-BH-02` | A2 | `ASTRA-DG-01` | D1 |
| `ASTRA-BH-03` | N5 | `ASTRA-DG-02` | N9 |
| `ASTRA-BH-04` | N6 | `ASTRA-DG-03` | N8 |

## 13. 其它四套旧编号的换算

| 文件 | 其编号 | 实际机制 | 新编号 |
|---|---|---|---|
| [第一轮 BSE](astra_physics_state_audit.md) | B1 / B2 / B3 / B4 / B5 | `dtm` / `tscale` / 并合漏锁 / `menv` / 孤立能量 | `BSE-01` / `BSE-02` / `BSE-03` / **`BSE-06`** / **`EN-01`** |
| [第一轮 SDAR](astra_sdar_silent_error_audit.md) | A / B / C | 组并行事件 / 邻居标志 / 稀疏 mask | `AR-01` / `H4-01` / `H4-02` |
| [第三轮 SDAR](astra_round3_sdar_state.md) | R3-S1 / R3-S2 | 递归 destroy / 分类漏记 | `AR-02` / `AR-04` |
| [第三轮 恒星](astra_round3_stellar_state.md) | S3-1 / S3-2 / S3-3 | BH 单位 / 异步 dm / RNG | `BH-01` / `BSE-05` / `IO-03` |
| [第三轮 能量](astra_round3_energy_state.md) | E3-1 / E3-2 / E3-3 | CM 参考系 / 逃逸能量 / 越界读 | `EN-03` / `EN-02` / `EN-04` |

> 注意：`S3-x`、`E3-x` 中的数字是**轮次编号**，不是新体系的分区号，勿与 `ASTRA-TR-*`/`ASTRA-EN-*` 混读。

## 14. 来源轮次索引（不参与排序）

同一批条目的审阅出处。**它只回答“哪一轮发现的”，不决定条目在正文中的位置**；正文位置只由等级与分区决定。

| 来源标注 | 项数 | 范围 |
|---|---|---|
| 旧别名 | 42 | 第一至五轮（`A/B/C/D/E/N/S` 旧编号） |
| 第六轮 | 11 | 主循环检查点、工具、能量账本 |
| 第七轮 | 7 | `SDAR/src` 11,867 行逐行补审 |
| 第八轮 | 24 | FDPS 默认构建编译面 34 文件 / 19,146 行 |
| 第九轮 | 17 | `bhmerger` 9,034 行逐行 |
| 第十轮 | 28 | `mcluster` 13,307 行逐行 |
| 第十一轮 | 14 | `bse-interface` 默认构建 11,070 行 |
| 第十二轮 | 37 | PeTar 主程序编译面 34 文件 / 11,771 存活行 |
| 第十三轮 | 12 | 默认构建剩余存活行收口 |
| **合计** | **192** | |

> 旧别名（`A1`–`E4`、`B1`–`B10`、`N1`–`N10`、`S1`–`S7`）来自第一至第五轮，逐条换算见[§12](#12-旧编号对照)与[§13](#13-其它四套旧编号的换算)；别名中的数字**不是**新体系的等级或分区号。

## 15. 尚未亲核候选（不占编号）

以下条目来自本轮三个子代理的逐行报告，已有逐字代码引用，但主代理尚未独立复核，列此备查；确认后再按分区入册。原始证据见 [`.scratch/astra_round12_20261027/registry_candidates.md`](../../.scratch/astra_round12_20261027/registry_candidates.md) 与各子代理报告。

- `profile.hpp`：`[Min/Max]` 两列口径不同（D4）；`n_loop==1` 整数除法（D4，除数为 1 故无实际影响）；`print()` 死代码除零（D4）；计时器内 16 次 MPI allreduce（D3，多 rank + `n_glb==1`）。（`end()` 双取时钟与跨成员指针算术已入册 `DG-09`/`DG-10`。）
- `hard_assert.hpp`：四条已入册 `DG-11..14`；头内 `static hard_dump` 见 `BL-07`。**原"`n_arti==0 && n_group>0` 时 dump 与 reader 计数不一致"经第十三轮逐字段核对推翻**：写/读顺序与元素个数对称，`n_arti==0` 时只是不读人工粒子数组；真实缺陷改记为 `BL-10`（索引无界、`initial()` 泄漏、先 resize 后校验）。
- `kickdriftstep.hpp`：单步 API 无调用者、空表时 `%= size()` 除零（D3）。（`setDKDMode` 互换已入册 `TR-18`。）
- `static_variables.hpp`：头内非 inline 静态定义（D3，单 TU 下无影响）；`group_data_mode=none` 哨兵（D3）；`eps=0/r_out=0` 毒默认（D3）。
- `hermite_information.hpp`/`hermite_perturber.hpp`：`time_origin` 无初始化（D4）；空 perturber 壳（已在 `BL-06`）。
- `domain.hpp`：`DomainDecision` 死代码与 `n_proc%4` 静默跳过已入册 `GEN-02/03`；4 重折叠 + 镜像域盒从不校验粒子是否在盒内（D3，需接入才可达）。
- `pseudoparticle_multipole.hpp`：`sqrt((1-e²)/4)` 对 e≥1 无守卫（D3）。
- `particle_distribution_generator.hpp`：`coef` 缺 `sqrt(GM)` 与 Kepler 盘 G=1 已入册 `GEN-04/05`；`r_cutoff` 用字面 -0.25（D4）；CM 归零只对本 rank（D4）；`h` 用互 Hill 半径（D4）；`mass_planet_glb` 未乘 `f_ice`（D4）。
- `bhmerger_output.hpp`：四条已入册 `BH-27..30`。
- `cluster_list.hpp`：组成员「成员位置 + CM 质量/速度」混合判据（潜在 D2/D1，疑点），第十三轮未定级。

---

## 16. 宏更正与阴性结论（第十三轮）

**宏更正（影响此前若干判断）**：`make -n -B build/petar.mpi.omp.avx2.bse.tt2nd` 实测默认宏含 `-D HARD_DUMP` 与 `-D AR_DEBUG_DUMP`。因此 `hard.hpp:2507` 的整簇 `backup`、`symplectic_integrator.h:2057` 的 `backupIntData` 都在默认路径上（仅开销、不改物理，未量化）。`ASSERT` 仍为空（`HARD_DEBUG` 未定义）。

**阴性/降级**：`n_arti` 计数不一致不成立（§15 已改）；`hermite_interaction.hpp` 六个 pair kernel 与 `calcEnergy` 逐式对照通过；`kickdriftstep` continue 状态机与主循环顺序自洽；`CalcForcePPNoSimd/CalcForcePPSimd` 无调用者（`eps2=0` 不入册）；`calcDtLimit` 的 `bitmap==0` 与第六轮同判（调用点已注释，不单列）。**已核对阴性（原文）**：`two_body_tide.hpp:137-144` 的 Mardling 修正指数 `0.5 + 0.25·(0.5|η−2|)^1.5` 与 Mardling & Aarseth 2001（MNRAS 321, 398, DOI `10.1046/j.1365-8711.2001.03974.x`）式 (A4) `ζ=η[2/(1+e)]^{α/2}` + 式 (A5) `α=1+½|(η−2)/2|^{3/2}` **逐位一致**（2,000 点最大差 `1.11e-16`，`lit/compare_ma2001.py`）；η 的定义与 `eccentricity_switch=0.9` 也对上原文适用范围。**疑点（不占编号）**：`evolveOrbitDynamicalTide` 的 `ASSERT(peri>R1+R2)` 默认空，近接触外推仅由 `e=0.9` 交换点截断。
