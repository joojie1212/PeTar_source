#pragma once
#include <utility>
#include <vector>
#include "ptcl.hpp"

// Tree forces represent resolved groups through artificial particles. A kick
// or output writeback may already have restored the members' physical masses.
// Temporarily restore the tree representation without changing group metadata.
// Arrays must not be resized while this guard is alive.
class TreeSourceMassGuard {
    std::vector<std::pair<Ptcl*, PS::F64>> saved_;
public:
    TreeSourceMassGuard() { assert(Ptcl::group_data_mode==GroupDataMode::artificial); }
    TreeSourceMassGuard(const TreeSourceMassGuard&)=delete;
    TreeSourceMassGuard& operator=(const TreeSourceMassGuard&)=delete;

    template<class Tarray>
    void mask(Tarray& particles, const PS::S32 count) {
        for(PS::S32 i=0;i<count;++i) {
            auto& p=particles[i];
            // Binaries without artificial particles remain real tree sources.
            if(p.group_data.artificial.isMember() && p.getParticleCMAddress()>0) {
                saved_.emplace_back(&p,p.mass);
                p.mass=0.0;
            }
        }
    }

    ~TreeSourceMassGuard() {
        for(const auto& entry:saved_) entry.first->mass=entry.second;
    }
};
