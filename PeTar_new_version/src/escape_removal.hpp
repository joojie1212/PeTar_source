#pragma once
#include <particle_simulator.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

// A synchronized real particle, never an artificial group/CM particle.
// Trivially copyable so FDPS can exchange candidates across ranks.
struct EscapeCandidate {
    PS::S64 id=0, partner=0;
    PS::F64 mass=0;
    PS::F64vec pos=0.0, vel=0.0;
    int forced=0;
};

namespace EscapeRemoval {
inline void sumPotentials(const std::vector<double>& local, std::vector<double>& global) {
    if(local.empty()) return;
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Allreduce(local.data(),global.data(),int(local.size()),MPI_DOUBLE,MPI_SUM,PS::Comm::getCommunicator());
#else
    global=local;
#endif
}

inline double pairPotential(const EscapeCandidate& a, const EscapeCandidate& b,
                            double G, double eps2) {
    const auto dr=a.pos-b.pos;
    const double r2=dr*dr+eps2;
    if (!(r2>0)) throw std::runtime_error("Coincident escape candidates without softening");
    return -G/std::sqrt(r2); // multiply by source mass for specific potential
}

inline std::vector<EscapeCandidate> gather(std::vector<EscapeCandidate>& local) {
    const int n=local.size(), nr=PS::Comm::getNumberOfProc();
    std::vector<int> counts(nr), offsets(nr+1,0);
    PS::Comm::allGather(&n,1,counts.data());
    for(int i=0;i<nr;++i) offsets[i+1]=offsets[i]+counts[i];
    std::vector<EscapeCandidate> global(offsets.back());
    // The non-MPI FDPS wrapper also writes offsets[1].
    EscapeCandidate dummy;
    PS::Comm::allGatherV(n?local.data():&dummy,n,
                        global.empty()?&dummy:global.data(),counts.data(),offsets.data());
    std::sort(global.begin(),global.end(),[](const EscapeCandidate& a,const EscapeCandidate& b){return a.id<b.id;});
    for(size_t i=1;i<global.size();++i)
        if(global[i-1].id==global[i].id) throw std::runtime_error("Duplicate real particle ID in escape candidates");
    return global;
}

struct Decision {
    std::vector<int> remove, forced, mate;
    double energy_change=0.0;
};

// phi includes current potentials from all real particles and any external
// field. It is NEVER read from the previous tree's per-particle potential cache.
inline Decision decide(const std::vector<EscapeCandidate>& p,
                       const std::vector<double>& phi, double radius2,
                       bool check_energy, double G, double eps2,
                       const PS::F64vec& center_pos, const PS::F64vec& center_vel) {
    Decision d;const size_t n=p.size();
    d.remove.assign(n,0);d.forced.assign(n,0);d.mate.assign(n,-1);
    std::unordered_map<PS::S64,size_t> ids;
    for(size_t i=0;i<n;++i) ids.emplace(p[i].id,i);
    // Only reciprocal, present, currently bound companions are protected.
    // Stale IDs after disruption must not glue unrelated single particles.
    for(size_t i=0;i<n;++i) {
        auto it=ids.find(p[i].partner);
        if(it==ids.end() || it->second==i) continue;
        size_t j=it->second;
        if(p[j].partner!=p[i].id) continue;
        auto dv=p[i].vel-p[j].vel;
        if(.5*(dv*dv)+(p[i].mass+p[j].mass)*pairPotential(p[i],p[j],G,eps2)<0)
            d.mate[i]=j;
    }
    for(size_t i=0;i<n;++i) {
        int j=d.mate[i];
        if(j>=0 && size_t(j)<i) continue;
        double mass=p[i].mass, pot=phi[i];
        auto pos=p[i].pos,vel=p[i].vel;
        int forced=p[i].forced;
        if(j>=0) {
            mass+=p[j].mass;
            pos=(p[i].mass*p[i].pos+p[j].mass*p[j].pos)/mass;
            vel=(p[i].mass*p[i].vel+p[j].mass*p[j].vel)/mass;
            // External potential of the pair: omit its own binding energy.
            // Weighted member evaluations avoid a point-CM approximation.
            double u=p[i].mass*p[j].mass*pairPotential(p[i],p[j],G,eps2);
            pot=(p[i].mass*phi[i]+p[j].mass*phi[j]-2*u)/mass;
            forced=forced||p[j].forced;
        }
        const auto dr=pos-center_pos,dv=vel-center_vel;
        const bool escape=forced || ((dr*dr)>=radius2 && (!check_energy || .5*(dv*dv)+pot>=0));
        d.remove[i]=escape;d.forced[i]=forced;
        if(j>=0) {d.remove[j]=escape;d.forced[j]=forced;}
    }
    return d;
}

// delta E = -(K_S + U_SR + U_SS + U_external,S).
// Sum(m_i phi_i) contains U_SS twice. Correct once, even across MPI ranks.
inline double energyChange(const std::vector<EscapeCandidate>& p,
                           const std::vector<double>& energy_phi,
                           const std::vector<int>& removed, double G, double eps2,
                           const PS::F64vec& inertial_velocity_offset) {
    double delta=0;
    for(size_t i=0;i<p.size();++i) if(removed[i]) {
        auto v=p[i].vel+inertial_velocity_offset;
        delta-=p[i].mass*(.5*(v*v)+energy_phi[i]);
        for(size_t j=0;j<i;++j) if(removed[j])
            delta+=p[i].mass*p[j].mass*pairPotential(p[i],p[j],G,eps2);
    }
    return delta;
}
} // namespace EscapeRemoval
