#pragma once
#include <vector>
#include "hard_ptcl.hpp"

// Event-local bookkeeping. No particle/dump layout changes; callbacks are owned
// by the active hard integrator, not by the shared ARInteraction manager.
struct MergerPartitionEvent {
    PtclHard *first, *second;
    PtclHard before_first, before_second;
    double time;
    std::vector<PtclHard> state;
    std::vector<unsigned char> in_group;
    double origin[3] = {0,0,0};
    double group_potential = 0;
    int first_index = -1, second_index = -1;
    struct Context {
        void* owner = nullptr;
        void (*begin)(void*, MergerPartitionEvent&) = nullptr;
        void (*finish)(void*, MergerPartitionEvent&) = nullptr;
    };
    static Context& context() { static thread_local Context value; return value; }
    struct Binding {
        Context saved;
        Binding(Context next): saved(context()) { context()=next; }
        ~Binding() { context()=saved; }
        Binding(const Binding&)=delete;
        Binding& operator=(const Binding&)=delete;
    };
    Context bound;
    MergerPartitionEvent(PtclHard& a, PtclHard& b, double t, bool enabled):
        first(&a), second(&b), before_first(a), before_second(b), time(t),
        bound(enabled ? context() : Context{}) {
        if (bound.begin) bound.begin(bound.owner,*this);
    }
    ~MergerPartitionEvent() {
        if (bound.finish && before_first.mass>0 && before_second.mass>0
            && (first->mass==0 || second->mass==0)) bound.finish(bound.owner,*this);
    }
    MergerPartitionEvent(const MergerPartitionEvent&)=delete;
    MergerPartitionEvent& operator=(const MergerPartitionEvent&)=delete;
};
