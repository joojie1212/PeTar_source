# PeTar：N × 初始半质量半径扫描

本目录可单独上传，不依赖被 Git 忽略的 `run_bh_cluster/`。依据当前实际
`run_n_series.sh`、单模型 runner 和 HARD_REOPT supervisor 制作；旧 README 的
100 trh/0.2 trh/72226 与实际入口不一致，这里遵循实际入口的 300/0.5/72222。

## 模型与资源

- N = 500, 1000, 2000, 4000, 10000。
- N_ref=10000 时 rh_ref = 0.0013047672, 0.013047672, 0.13047672 pc。
- rh = rh_ref (N/10000)^(1/3)，每档保持解析中心密度相同；随机有限粒子实现会有涨落。
- 单颗 30 Msun、Plummer、无初始双星、同一 seed=72222、`petar.init --blackhole`
  零自旋初始化，速度由 km/s 转换为 pc/Myr。
- G=0.00449830997959438，Coulomb gamma=0.11，运行 300 初始 trh，每 0.5 trh 输出。
- 逃逸半径=100 初始 rh；其正值仍要求满足程序的能量判据。
- 原积分参数保持：-a 1、-u 1、--write-group-info 0、--stellar-evolution 2、
  --energy-err-hard 1e-3；HARD_REOPT ratio/windows/cooldown=0.9/3/10。
- 一节点、16 CPU 核、最多四个模型并发；每模型一个 MPI rank、4 OpenMP 线程。
  每个独立 srun step 使用 --exclusive --exact -n1 -c4 --cpu-bind=cores --mpi=none。
  在每个 step 内用 OpenMPI 4 的 mpirun --mca plm isolated --mca ras ^slurm
  --bind-to none -np 1 启动本地单 rank，继承 step 的四核亲和性。
  主模拟和 find.dt 使用相同启动器，避免依赖服务器缺失的 Slurm PMI 直接启动支持，
  也避免 mpirun 再申请嵌套 Slurm step。服务器上仍需短作业验证。
- 不执行 petar.data.gether 或 petar.data.process，不生成自动后处理派生文件。
- 初始 input 仍由 petar.init 生成 ASCII。首轮用 -i 2（读 ASCII、写二进制）；
  HARD_REOPT 重启用 -i 0（读写二进制），petar.find.dt 同样使用 -i 0。
  -a 1 保留原含义，不是格式开关。状态、并合记录和运行日志仍为文本。
  二进制文件应使用匹配构建的 PeTar/格式转换工具读取。

## 提交

先在服务器确认 `~/petar_env.sh` 可加载 MPI、petar-runtime、PeTar 与 McLuster。
需要 petar.format.transfer，供 petar.find.dt 读取二进制快照头。

```bash
cd "$HOME/PeTar_source/server/n_rh_sweep"
source "$HOME/petar_env.sh"
bash run_n_rh_series.slurm --print-config
sinfo -p tiny -o '%P %l %c %m'
scontrol show partition tiny
```

脚本默认请求 `tiny`、16G、5 天。这是初始资源预算，不是实测完成时间。
如果 tiny 不允许 16 核/5 天，应改用实际允许的分区，例如：

```bash
# 将 YOUR_PARTITION 替换为服务器上真实的合适分区
sbatch --partition=YOUR_PARTITION --time=5-00:00:00 run_n_rh_series.slurm --run
```

如果 tiny 的限制满足需求，直接 `sbatch run_n_rh_series.slurm --run`。
请从本目录提交；从别处提交时设置 SWEEP_SCRIPT_DIR 为本目录绝对路径。
`--print-config` 不创建模型/目录，也不提交作业。没有 --run 不启动模拟。
可设置 OUTPUT_ROOT、MCLUSTER_SEED、PETAR_ENV_FILE。保持资源参数为一节点/4 tasks/4 CPUs。

默认输出 `$HOME/scratch/petar_n_rh_results/batch_<jobid>_seed72222_<随机后缀>/`。
各模型生成独立目录，manifest.tsv 给出半径、密度、trh、终止时间，status.tsv 记录
queued/launched/completed/failed（launched 包含 Slurm step 等待时间）。
每个任务退出后自动补位；单任务失败会记录并继续其余模型，最终批作业返回非零。
提交时保留配套脚本副本，并使用副本执行后续任务，防止运行中修改源脚本影响批次。
resources.txt 记录 GNU time 的墙钟时间及进程 RSS，workflow.log 记录全过程。
GNU time 的 MaxRSS 不能当作多个并发 Python worker 的同时总内存；需要结合 Slurm 监控。

```bash
squeue -u "$USER"
sacct -j JOBID --format=JobID,JobName%36,State,Elapsed,AllocCPUS,MaxRSS,ExitCode
```

调度器超时/取消不是 HARD_REOPT 的退出码 42，不会由此脚本自动恢复。
重新提交会创建新批次，不会自动续跑旧批次。异常强制终止时 status.tsv 可能停在 launched；
以 sacct 记录为准。普通 HARD_REOPT 的快照重启和 find.dt 流程保留。

## 内存、时间与磁盘估算

没有目标服务器上的 wall-time/MaxRSS 基准，不能承诺具体完成时间。
初始建议每个 step 4G、四个并发共 16G，包含 Python/SciPy、树结构、紧密组积分和
模拟运行的余量；不是测量值，也不是严格上限保证。每个模拟线程的 OMP_STACKSIZE=128M，
4 线程的栈配置量级约 512M（不代表实际常驻占用）。10000 粒子、37 列 double 的基础
快照数组只有约 2.82 MiB，不能据此推断程序总内存。若监测后需要每任务 8G，必须同时
将 Slurm 总申请改为 32G，并修改 sweep.py 中 step 的 --mem=8G，否则可能只能并发两个任务。

弛豫时间公式：
trh=0.138 sqrt(N) rh^(3/2) / [sqrt(G*m) ln(0.11N)]。
相同 N 下，三档物理 trh/终止时间比例为 1:31.6228:1000；这不是墙钟耗时比例。
牛顿同构模型按 trh 演化时步数可以相近，但这里保留引力波、并合和紧密双星演化，
不同密度下耗时可能显著不同，不能仅用物理终止时间外推。

若只为制定预算，假设每步代价约 O(N) 且每 trh 步数随 N/ln(0.11N) 增长，
则相同密度档的总代价约 N²/ln(0.11N)。相对 N=500 的倍数约：
1、3.41、11.89、42.14、228.89。这只是简化估算，忽略树的对数项、负载变化和紧密组。
建议依据各密度档 N=500/1000 的服务器实测调整时间；5 天是申请时限，不能保证
N=10000 在此时限内完成。完整 15 个任务总墙钟至少为 max(最长单任务时间, 各任务时间和/4)，
还要加作业内 step 等待、I/O、并发竞争及尾部空闲时间。tiny 示例的 30 分钟不能作为本批次的完成预算。

每模型约 601 个常规快照，所有模型按初始粒子数保持不变计共约 3155 万粒子记录。
二进制输出大小取决于本次构建的原生结构体布局，不能用 ASCII 的 867 字节/粒子直接估算。
取消后处理也避免生成 .single/.binary 等派生数据。可在服务器输出首个快照后，以其大小
乘 601 估算该模型的常规快照占用，再汇总 15 个模型；逃逸/并合会减少后期粒子数。
HARD_REOPT 额外快照和文本日志仍需另留空间。尤其 workflow.log 和各 cycle 日志仍重复
保存模拟输出，五天运行时应监测日志及存储配额。

## 本地验证

```bash
python3 test_sweep.py
```

使用模拟 srun（不运行真实模拟）验证：15 个唯一模型、同档密度不变、最高四任务并发、
4 CPU/4G step 参数、失败后继续补位和非零批次退出码。
增加模拟 HARD_REOPT 测试，检查首轮 ASCII→二进制、优化读取二进制及重启读写二进制。
未提交服务器作业，未进行 300 trh 性能测试。

## 可选逃逸半径（质量质心修复后）

需要重新编译安装修复后的 PeTar；仅更新脚本不会修复中心判据。
默认仍为 100 初始 rh，可用以下命令选择 200：

```bash
ESCAPE_RADIUS_RH=200 bash run_n_rh_series.slurm --print-config
ESCAPE_RADIUS_RH=200 sbatch run_n_rh_series.slurm --run
```

每个模型的 config.txt 保存实际半径。增大半径会推迟移除并增加计算负担，
不保证更高数值精度；建议先对比 100/200/300 的演化和并合统计。
