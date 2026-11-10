#pragma once
#include <vector>
#include "hard_ptcl.hpp"

// Event-local bookkeeping for mergers and non-merging orbital updates.
// No particle/dump layout changes; callbacks belong to the active hard
// integrator, not to the shared ARInteraction manager.
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
        if (!bound.finish) return;
        // A velocity-only update is already covered by AR/CM kinetic accounting.
        auto changed=[](const PtclHard& a,const PtclHard& b) {
            return a.mass!=b.mass || a.pos.x!=b.pos.x
                || a.pos.y!=b.pos.y || a.pos.z!=b.pos.z;
        };
        if (changed(before_first,*first) || changed(before_second,*second))
            bound.finish(bound.owner,*this);
    }
    MergerPartitionEvent(const MergerPartitionEvent&)=delete;
    MergerPartitionEvent& operator=(const MergerPartitionEvent&)=delete;
};
