#pragma once

#include <cstdio>
#include <cstdlib>
#include <omp.h>

//! Intra-cluster (level-2) parallelism of the hard force loops -- v2, scoped ICV
/*!
  v1 的致命错误：在 main() 里调用 omp_set_max_active_levels(2)，这是**程序级**的开关。
  FDPS 的树代码里有 "parallel for 写在 parallel 区域内部" 的写法，长期依赖
  "内层被自动降级成 1 线程"。v1 打开 level 2 后，每个外层线程又开一队，最高 32x32=1024
  线程，导致 Tree_NB 劣化 215 倍、PP_cluster 劣化到 200 倍。

  v2 的做法：只在**一个 hard cluster 任务**的生命期内打开 level 2，退出立即恢复 1。
  这样 FDPS 的树代码永远看不到 level 2（已验证：cluster 任务内不可达任何树代码）。

  环境变量：
    PETAR_INNER_THREADS  k   每个 cluster 任务内部的线程数（1=关闭，默认 1）
    PETAR_OUTER_THREADS  m   "每集团一个任务"这一层的线程数（默认 = OMP_NUM_THREADS）
    PETAR_PAR_H4         0/1 是否并行 calcAccJerkNBList（默认 1）
    PETAR_PAR_AR         0/1 是否并行 AR 组循环（默认 1）
    PETAR_PARSTAT        0/1 退出时打印 [PARSTAT]（默认 1）

  这些开关只影响性能，不改变浮点运算顺序，因此结果应与串行逐位一致。
 */
namespace InnerPar {

    inline int n_var     = 0;
    inline int outer_var = 0;
    inline int use_h4    = 1;
    inline int use_ar    = 1;
    inline int stat_on   = 1;
    inline int use_tasks = 0;   // PETAR_PAR_MODE=task：用 taskloop 复用现有 team，不建嵌套 team
    inline int use_lpt    = 1;   // PETAR_AR_LPT=0：关闭 AR 组循环的"大组优先"(LPT) 调度

    //! 诊断计数器：region 数与被切分的条目数（用于回答"粒度够不够粗"）
    inline long long n_region = 0;
    inline long long n_item   = 0;

    //! AR 组循环 task 化 / BSE 串行化的诊断计数器（P1, 2026-10-18）
    inline long long n_ar_region = 0;
    inline long long n_ar_item   = 0;
    inline long long n_bse_call  = 0;
    inline double    t_bse       = 0.0;
    //! AR 组循环每次调用的组数分布（b0:<8, b1:8-31, b2:32-127, b3:>=128）
    inline long long n_ar_calls  = 0;
    inline long long n_ar_work   = 0;
    inline long long n_ar_b[4]   = {0,0,0,0};
    //! [P2] item 级代价分布（per-thread，避免原子）
    inline const int NP_IT_TH = 64;
    inline double    t_item[NP_IT_TH][3] = {};
    inline long long n_item_th[NP_IT_TH][3] = {};
    inline double    t_item_max[NP_IT_TH] = {};
    inline long long item_max_mem[NP_IT_TH]  = {};
    inline long long item_max_kind[NP_IT_TH] = {};
    inline long long item_max_tot[NP_IT_TH]  = {};
    inline long long n_member_hist[8] = {};   // resolved group 的成员数直方图
    inline long long sum_member_res = 0;
    inline long long n_mem_region = 0;        // member 循环被 task 化的 region 数/成员数
    inline long long n_mem_item   = 0;
    //! region 级效率：sum(ideal)/sum(wall)，ideal = Σitem_time / n_threads
    inline double    sum_reg_wall  = 0.0;
    inline double    sum_reg_ideal = 0.0;
    inline double    sum_reg_item  = 0.0;
    //! item 代价 vs 邻居数（回答"item 代价是否 ∝ 邻居数"）
    inline long long n_nb_th[NP_IT_TH][10] = {};
    inline double    t_nb_th[NP_IT_TH][10] = {};
    inline double    t_regwall_th[NP_IT_TH] = {};
    inline void noteItemNeighbor(const long long _n_nb, const double _dt) {
        if (!stat_on) return;
        const int ith = omp_get_thread_num();
        const int it = (ith<NP_IT_TH) ? ith : NP_IT_TH-1;
        const int ib = (_n_nb<=2)?0:(_n_nb<=4)?1:(_n_nb<=8)?2:(_n_nb<=16)?3:(_n_nb<=32)?4:
                       (_n_nb<=64)?5:(_n_nb<=128)?6:(_n_nb<=256)?7:(_n_nb<=512)?8:9;
        n_nb_th[it][ib] += 1;
        t_nb_th[it][ib] += _dt;
    }
    //! region 墙钟（per-thread，无原子）
    inline void noteRegionWall(const double _wall) {
        if (!stat_on) return;
        const int ith = omp_get_thread_num();
        t_regwall_th[(ith<NP_IT_TH)?ith:NP_IT_TH-1] += _wall;
    }

    //! region 尺寸直方图（n_tot 分桶）与每桶的 item CPU
    inline long long n_reg_size[6]  = {0,0,0,0,0,0};   // <32,32-63,64-127,128-255,256-511,>=512
    inline double    sum_item_size[6] = {0,0,0,0,0,0};
    inline double    sum_wall_size[6] = {0,0,0,0,0,0};
    //! 注意：一个 region 只有一次调用，所以这里的原子不构成热点（与 per-item 原子不同）
    inline void noteRegionSize(const int _n_tot, const double _wall, const double _sum_item) {
        if (!stat_on) return;
        const int ib = (_n_tot<32)?0:(_n_tot<64)?1:(_n_tot<128)?2:(_n_tot<256)?3:(_n_tot<512)?4:5;
#pragma omp atomic
        n_reg_size[ib]++;
#pragma omp atomic
        sum_item_size[ib] += _sum_item;
#pragma omp atomic
        sum_wall_size[ib] += _wall;
    }

    //! 每个 region 实际参与的线程数（按 n_tot 分桶：<64 与 >=64）
    inline double    sum_thr_big = 0.0;   // >=64 个 item 的 region：Σ参与线程数
    inline long long n_reg_big   = 0;
    inline double    sum_thr_small = 0.0;
    inline long long n_reg_small   = 0;
    inline double    sum_wall_big = 0.0;
    inline double    sum_item_big = 0.0;
    inline double    sum_wall_small = 0.0;
    inline double    sum_item_small = 0.0;

    inline void noteRegionThreads(const int _n_threads, const int _n_tot,
                                  const double _wall, const double _sum_item) {
        if (!stat_on) return;
        if (_n_tot>=64) {
#pragma omp atomic
            sum_thr_big += _n_threads;
#pragma omp atomic
            n_reg_big++;
#pragma omp atomic
            sum_wall_big += _wall;
#pragma omp atomic
            sum_item_big += _sum_item;
        }
        else {
#pragma omp atomic
            sum_thr_small += _n_threads;
#pragma omp atomic
            n_reg_small++;
#pragma omp atomic
            sum_wall_small += _wall;
#pragma omp atomic
            sum_item_small += _sum_item;
        }
    }

    inline void noteRegionStat(const double _sum_item, const double _wall, const int _nth) {
        if (!stat_on) return;
        const double ideal = _sum_item / (double)(_nth>0?_nth:1);
#pragma omp atomic
        sum_reg_wall  += _wall;
#pragma omp atomic
        sum_reg_ideal += ideal;
#pragma omp atomic
        sum_reg_item  += _sum_item;
    }

    //! n_group_tot 的精确分布（0..7）
    inline long long n_ar_hist[8] = {0,0,0,0,0,0,0,0};
    //! 空组循环 vs 有组循环的耗时
    inline double    t_ar_empty = 0.0;
    inline double    t_ar_work  = 0.0;
    inline int       ar_min_work = 2;   // PETAR_AR_MINWORK（实测每组循环只有 1-4 个活跃 AR 组）
    inline int       use_mem    = 0;   // PETAR_PAR_MEM：resolved group 的 member 循环 task 化
    inline int       mem_min    = 16;  // PETAR_MEM_MIN：成员数下限
    inline int       mem_per_task = 8; // PETAR_MEM_PER_TASK：每个 task 的成员数
    inline int       h4_task_mult = 1; // PETAR_H4_TASKS：H4 item 循环的 num_tasks 倍数
    inline int       h4_ipt      = 0;  // PETAR_H4_IPT：每个 task 的目标 item 数（>0 启用尺寸自适应）
    inline int       item_stat   = 0; // PETAR_ITEM_STAT：item 级计时（会扰动调度，默认关）
    inline int       cluster_lpt  = 0; // PETAR_CLUSTER_LPT：大 cluster 优先分发（实测无收益，默认关）

    //! 串行段计时（PETAR_PHASE_TIMER=1）：报告最慢 cluster 的 all/AR/single/force
    inline int phase_on = 0;
    inline double g_phase_all = 0, g_phase_ar = 0, g_phase_single = 0, g_phase_force = 0;
    inline void phaseReport() {
        if (!phase_on) return;
        fprintf(stderr,
                "[PHASE] all=%.6f ar=%.6f single=%.6f force=%.6f "
                "other=%.6f predict_correct=%.6f\n",
                g_phase_all, g_phase_ar, g_phase_single, g_phase_force,
                g_phase_all - g_phase_ar - g_phase_single,
                g_phase_single - g_phase_force);
    }

    inline int envInt(const char* _name, const int _dflt) {
        const char* e = std::getenv(_name);
        if (e == nullptr) return _dflt;
        return std::atoi(e);
    }

    inline void report() {
        if (!stat_on) return;
        if (n_region > 0)
            fprintf(stderr, "[PARSTAT] regions=%lld items=%lld\n", n_region, n_item);
        {
            long long cnt[3] = {0,0,0};
            double    tim[3] = {0.0,0.0,0.0};
            double tmax = 0.0; long long mmax=0, kmax=-1, rmax=0;
            for (int it=0; it<NP_IT_TH; it++) {
                for (int k=0;k<3;k++) { cnt[k] += n_item_th[it][k]; tim[k] += t_item[it][k]; }
                if (t_item_max[it] > tmax) {
                    tmax = t_item_max[it]; mmax = item_max_mem[it];
                    kmax = item_max_kind[it]; rmax = item_max_tot[it];
                }
            }
            if (cnt[0]+cnt[1]+cnt[2] > 0) {
                fprintf(stderr, "[PARITEM] items single=%lld res=%lld cm=%lld | "
                        "time single=%.3f res=%.3f cm=%.3f s | mean_members_res=%.1f\n",
                        cnt[0], cnt[1], cnt[2], tim[0], tim[1], tim[2],
                        cnt[1]>0 ? double(sum_member_res)/double(cnt[1]) : 0.0);
                fprintf(stderr, "[PARITEM] h4_task_mult=%d (num_tasks = %d x nthreads) cluster_lpt=%d\n",
                        h4_task_mult, h4_task_mult, cluster_lpt);
                {
                    static const long long ub[10] = {2,4,8,16,32,64,128,256,512,1000000000};
                    long long nn[10] = {0,0,0,0,0,0,0,0,0,0};
                    double    tt[10] = {0,0,0,0,0,0,0,0,0,0};
                    double    wall_all = 0.0;
                    for (int it=0; it<NP_IT_TH; it++) {
                        wall_all += t_regwall_th[it];
                        for (int ib=0; ib<10; ib++) { nn[ib]+=n_nb_th[it][ib]; tt[ib]+=t_nb_th[it][ib]; }
                    }
                    fprintf(stderr, "[PARITEM] item cost vs neighbor count:\n");
                    for (int ib=0; ib<10; ib++)
                        fprintf(stderr, "[PARITEM]   nb<=%-10lld n=%-8lld mean=%.2f us  total=%.3f s\n",
                                ub[ib], nn[ib], nn[ib]>0 ? tt[ib]/nn[ib]*1e6 : 0.0, tt[ib]);
                    fprintf(stderr, "[PARITEM] region wall (sum over calls) = %.3f s ; "
                            "effective concurrency = item_cpu / region_wall\n", wall_all);
                }
                fprintf(stderr, "[PARITEM] region size histogram: buckets(<32,32-63,64-127,128-255,256-511,>=512)\n");
                for (int ib=0; ib<6; ib++)
                    fprintf(stderr, "[PARITEM]   size_bucket[%d] n=%-6lld item_cpu=%7.3f s wall=%7.3f s "
                            "concurrency=%.1f\n", ib, n_reg_size[ib], sum_item_size[ib], sum_wall_size[ib],
                            sum_wall_size[ib]>0.0 ? sum_item_size[ib]/sum_wall_size[ib] : 0.0);
                fprintf(stderr, "[PARITEM] threads/region: big(>=64 items) n=%lld avg_thr=%.2f "
                        "avg_wall=%.3f ms avg_item=%.3f ms | small n=%lld avg_thr=%.2f\n",
                        n_reg_big, n_reg_big>0? sum_thr_big/n_reg_big : 0.0,
                        n_reg_big>0? sum_wall_big/n_reg_big*1e3 : 0.0,
                        n_reg_big>0? sum_item_big/n_reg_big*1e3 : 0.0,
                        n_reg_small, n_reg_small>0? sum_thr_small/n_reg_small : 0.0);
                fprintf(stderr, "[PARITEM] region_eff: sum_wall=%.3f s sum_ideal=%.3f s "
                        "sum_item_cpu=%.3f s -> eff=%.3f\n",
                        sum_reg_wall, sum_reg_ideal, sum_reg_item,
                        sum_reg_wall>0.0 ? sum_reg_ideal/sum_reg_wall : 0.0);
                fprintf(stderr, "[PARITEM] member_hist(<=2,<=4,<=8,<=16,<=32,<=64,<=128,>128)="
                        "%lld,%lld,%lld,%lld,%lld,%lld,%lld,%lld | "
                        "max_item=%.3f ms kind=%lld members=%lld region_items=%lld mem_par=%d "
                        "mem_regions=%lld mem_items=%lld\n",
                        n_member_hist[0], n_member_hist[1], n_member_hist[2], n_member_hist[3],
                        n_member_hist[4], n_member_hist[5], n_member_hist[6], n_member_hist[7],
                        tmax*1e3, kmax, mmax, rmax, use_mem, n_mem_region, n_mem_item);
            }
        }
        if (n_ar_calls > 0) {
            fprintf(stderr, "[PARSTAT-AR] calls=%lld work_sum=%lld "
                    "buckets(<8,8-31,32-127,>=128)=%lld,%lld,%lld,%lld "
                    "ar_regions=%lld ar_items=%lld bse_calls=%lld bse_time=%.6f\n",
                    n_ar_calls, n_ar_work, n_ar_b[0], n_ar_b[1], n_ar_b[2], n_ar_b[3],
                    n_ar_region, n_ar_item, n_bse_call, t_bse);
            fprintf(stderr, "[PARSTAT-AR] hist(n_group:count)=0:%lld 1:%lld 2:%lld 3:%lld "
                    "4:%lld 5:%lld 6:%lld 7:%lld | t_empty=%.6f t_work=%.6f min_work=%d\n",
                    n_ar_hist[0], n_ar_hist[1], n_ar_hist[2], n_ar_hist[3],
                    n_ar_hist[4], n_ar_hist[5], n_ar_hist[6], n_ar_hist[7],
                    t_ar_empty, t_ar_work, ar_min_work);
        }
    }

    inline void init() {
        if (n_var != 0) return;
        int n = envInt("PETAR_INNER_THREADS", 1);
        if (n < 1) n = 1;
        n_var   = n;
        use_h4  = envInt("PETAR_PAR_H4", 1);
        use_ar  = envInt("PETAR_PAR_AR", 1);
        stat_on = envInt("PETAR_PARSTAT", 1);
        use_lpt = envInt("PETAR_AR_LPT", 1);
        ar_min_work = envInt("PETAR_AR_MINWORK", 2);
        use_mem      = envInt("PETAR_PAR_MEM", 0);
        mem_min      = envInt("PETAR_MEM_MIN", 16);
        mem_per_task = envInt("PETAR_MEM_PER_TASK", 8);
        h4_task_mult = envInt("PETAR_H4_TASKS", 1);
        h4_ipt       = envInt("PETAR_H4_IPT", 0);
        item_stat    = envInt("PETAR_ITEM_STAT", 0);
        cluster_lpt  = envInt("PETAR_CLUSTER_LPT", 0);
        phase_on = envInt("PETAR_PHASE_TIMER", 0);
        const char* mode = std::getenv("PETAR_PAR_MODE");
        use_tasks = (mode != nullptr && mode[0]=='t' && mode[1]=='a'
                     && mode[2]=='s' && mode[3]=='k');
        // 关键：这里**不**设 max_active_levels。level 2 只在 enterCluster() 内打开。
        fprintf(stderr,
                "[PAR] v5 scoped=1 mode=%s inner=%d omp_max=%d h4=%d ar=%d "
                "max_active_levels_at_start=%d\n",
                use_tasks ? "task" : "team", n_var, omp_get_max_threads(),
                use_h4, use_ar, omp_get_max_active_levels());
        atexit(report);
        atexit(phaseReport);
    }

    inline int n() { if (n_var == 0) init(); return n_var; }

    //! 进入一个 hard cluster 任务：仅在该任务作用域内允许 level 2
    inline void enterCluster() { if (n_var > 1 && !use_tasks) omp_set_max_active_levels(2); }

    //! 退出该任务：**必须**恢复，否则跑过 cluster 的线程会让下一轮树相位也看到 level 2
    inline void leaveCluster() { if (n_var > 1 && !use_tasks) omp_set_max_active_levels(1); }

    //! 树力入口处兜底：保证 FDPS 永远只看到 level 1
    inline void forceLevel1() { if (n_var > 0) omp_set_max_active_levels(1); }

    //! "每集团一个任务"那一层的队大小（PETAR_OUTER_THREADS，默认 = 默认队）
    inline int outerFor(const int _n_default) {
        if (outer_var == 0) {
            const char* e = std::getenv("PETAR_OUTER_THREADS");
            const int m = (e != nullptr) ? std::atoi(e) : 0;
            outer_var = (m >= 1) ? m : _n_default;
        }
        return (outer_var < _n_default) ? outer_var : _n_default;
    }

    //! 是否够粗到值得 fork/join，并受分项开关控制
    inline int nFor(const int _n_work, const int _min_work, const int _enable) {
        if (!_enable) return 1;
        const int nth = n();
        return (_n_work < _min_work) ? 1 : nth;
    }

    //! AR 组循环的 region 计数（只在遇到 region 的线程里调用，无需原子）
    //! 记录一次 AR 组循环调用的组数（跨线程，用原子）
    inline void noteArCall(const int _n_work) {
        const int ib = (_n_work<8) ? 0 : (_n_work<32) ? 1 : (_n_work<128) ? 2 : 3;
#pragma omp atomic
        n_ar_calls++;
#pragma omp atomic
        n_ar_work += _n_work;
#pragma omp atomic
        n_ar_b[ib]++;
        if (_n_work<8) {
#pragma omp atomic
            n_ar_hist[_n_work]++;
        }
    }

    //! [P2] 记录一个 item 的代价（kind: 0=single, 1=resolved group, 2=cm group）
    inline void noteItem(const int _kind, const int _members, const double _dt, const int _n_tot) {
        const int ith = omp_get_thread_num();
        const int it = (ith<NP_IT_TH) ? ith : NP_IT_TH-1;
        t_item[it][_kind] += _dt;
        n_item_th[it][_kind] += 1;
        if (_dt > t_item_max[it]) {
            t_item_max[it]  = _dt;
            item_max_mem[it]  = _members;
            item_max_kind[it] = _kind;
            item_max_tot[it]  = _n_tot;
        }
        if (_kind==1) {
            const int ib = (_members<=2)?0:(_members<=4)?1:(_members<=8)?2:(_members<=16)?3:
                           (_members<=32)?4:(_members<=64)?5:(_members<=128)?6:7;
#pragma omp atomic
            n_member_hist[ib]++;
#pragma omp atomic
            sum_member_res += _members;
        }
    }

    inline void noteMemberRegion(const int _n) {
        if (!stat_on) return;
#pragma omp atomic
        n_mem_region++;
#pragma omp atomic
        n_mem_item += _n;
    }

    //! H4 item 循环的 task 数
    /*!  两种口径：
         * PETAR_H4_IPT>0（推荐）：按"每个 task 多少 item"定，nt = ceil(n_item/IPT)，
           clamp 到 [1, 线程数] —— 小 region 自然只有 1-2 个 task，大 region 才吃满 32；
         * PETAR_H4_IPT=0：旧行为 num_tasks = 线程数 × PETAR_H4_TASKS。 */
    inline int h4TaskCount(const int _n_item) {
        const int t = n();
        if (h4_ipt > 0) {
            int nt = (_n_item + h4_ipt - 1) / h4_ipt;
            if (nt > t) nt = t;
            return (nt<1) ? 1 : nt;
        }
        int nt = t * (h4_task_mult>0 ? h4_task_mult : 1);
        if (nt > _n_item) nt = _n_item;
        return (nt<1) ? 1 : nt;
    }

    //! [P2] member 循环的 task 数；<=1 走串行
    inline int memberTaskCount(const int _n_member) {
        const int t = n();                      // 触发 init()
        if (!use_tasks || !use_mem || t<=1) return 1;
        if (!omp_in_parallel()) return 1;
        if (_n_member < mem_min) return 1;
        int nt = (_n_member + mem_per_task - 1) / mem_per_task;
        if (nt > 2*t) nt = 2*t;                 // 允许 2x 超订，压尾部
        if (nt > _n_member) nt = _n_member;
        return (nt<1) ? 1 : nt;
    }

    inline void noteArRegion(const int _n_item) {
        if (!stat_on) return;
        n_ar_region++;
        n_ar_item += _n_item;
    }

    //! task 模式下 AR 组循环的 task 数；<=1 表示走原串行路径
    inline int arTaskCount(const int _n_work, const int _min_work) {
        const int t = n();                 // 触发 init()
        if (!use_tasks || !use_ar) return 1;
        if (t <= 1 || _n_work < _min_work || _n_work < ar_min_work) return 1;
        return (t > _n_work) ? _n_work : t;
    }

    inline void noteRegion(const int _n_item) {
        if (!stat_on) return;
#pragma omp atomic
        n_region++;
#pragma omp atomic
        n_item += _n_item;
    }

}
