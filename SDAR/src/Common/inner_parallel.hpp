#pragma once

#include <omp.h>

//! Helpers for the intra-cluster parallelism of the hard loops.
/*!
  PeTar already parallelises the hard clusters: every cluster is one iteration
  of an OpenMP parallel loop.  The loops *inside* one cluster (the H4 force
  loop) are therefore split with `taskloop` statements: the tasks are executed
  by the team that is already there, so no nested team is created.  That is
  important because the FDPS tree code relies on nested OpenMP regions being
  serialised.

  Every loop that uses these helpers must satisfy two conditions for its result
  to be independent of the task schedule:
    1. an iteration writes only objects owned by that iteration;
    2. every value written by one iteration and read by another has to be taken
       from a snapshot that is built before the loop starts.
*/
namespace InnerPar {

//! True when the calling thread is already inside an OpenMP parallel region.
inline bool inTeam() { return omp_in_parallel()!=0; }

//! Number of tasks to generate for a loop of _n independent items.
inline int taskCount(const int _n) {
    const int n_thread = omp_get_num_threads();
    return (_n<n_thread) ? _n : n_thread;
}

//! Whether a loop of _n independent items should be split into tasks.
/*! A task loop costs a fork/join even when there is one item per thread, and
    the hard-cluster loops are entered thousands of times per step, so the
    items are only split when there are at least _min_per_task of them per
    thread.  Below that the plain loop is used. */
inline bool useTask(const int _n, const int _min_per_task = 2) {
    if (_n<=1 || !inTeam()) return false;
    return _n >= _min_per_task*omp_get_num_threads();
}

} // namespace InnerPar
