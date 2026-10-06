#pragma once
#include "bh_peters.hpp"

namespace BHGW {
// PeTar owns this orbit-averaged prescription. SDAR still advances the
// conservative orbit and perturbations between dissipative callbacks.
// For BH-BH only, time_interrupt=-(1+t_last_GW) marks an initialized GW
// clock. Legacy snapshots have nonnegative time_interrupt and initialize at
// their first callback, without replaying past BSE losses. SSE resets this
// marker when a component becomes a single. time_merger is a prediction only.
template<class Interaction, class Interrupt, class Binary>
bool evolve(Interaction& interaction, Interrupt& interrupt, Binary& bin,
            int& modified, const double cutoff_rg=10.0) {
    auto& p=*bin.getLeftMember();
    auto& q=*bin.getRightMember();
    const double now=interrupt.time_now;
    const double c=interaction.bse_manager.getSpeedOfLight();
    const double G=interaction.gravitational_constant;
    const double cutoff=cutoff_rg*G*(p.mass+q.mass)/(c*c);
    const bool same_pair=p.getBinaryPairID()==q.id && q.getBinaryPairID()==p.id;
    const bool clock_valid=same_pair && p.time_merger>=0 && q.time_merger>=0
        && p.time_interrupt<=-1 && q.time_interrupt<=-1;
    const double dt=clock_valid
        ? std::max(0.0,now-std::max(-1-p.time_interrupt,-1-q.time_interrupt)) : 0.0;
    p.setBinaryPairID(q.id);
    q.setBinaryPairID(p.id);
    bin.calcOrbit(G);
    bool reached=false;
    double predicted=-1;
    double lag=0;
    if (bin.semi>0 && bin.ecc>=0 && bin.ecc<1
        && std::isfinite(bin.semi) && std::isfinite(bin.ecc)) {
        const double b=BHPeters::beta(G,c,p.mass,q.mass);
        const auto result=BHPeters::advance({double(bin.semi),double(bin.ecc)},dt,b,cutoff);
        if (result.orbit.a!=bin.semi || result.orbit.e!=bin.ecc) {
            const double mass=p.mass+q.mass;
            double rcm[3],vcm[3];
            for (int k=0;k<3;++k) {
                rcm[k]=(p.mass*p.pos[k]+q.mass*q.pos[k])/mass;
                vcm[k]=(p.mass*p.vel[k]+q.mass*q.vel[k])/mass;
            }
            bin.semi=result.orbit.a;
            bin.ecc=result.orbit.e;
            // Retain the conservative eccentric anomaly and orbital plane.
            // This is operator splitting, not phase-resolved PN integration.
            bin.calcParticles(G);
            for (int k=0;k<3;++k) {
                p.pos[k]+=rcm[k];q.pos[k]+=rcm[k];
                p.vel[k]+=vcm[k];q.vel[k]+=vcm[k];
            }
            bin.calcOrbit(G);
            interrupt.adr=&bin;
            if (interrupt.status==decltype(interrupt.status)::none)
                interrupt.status=decltype(interrupt.status)::change;
            modified=2;
        }
        reached=result.reached_cutoff;
        if (reached) {
            predicted=now-dt+result.elapsed;
            lag=dt-result.elapsed;
        }
        else {
            const auto remaining=BHPeters::advance(result.orbit,
                std::numeric_limits<double>::infinity(),b,cutoff);
            predicted=now+remaining.elapsed;
        }
    }
    // An unbound orbit cancels its old bound-inspiral deadline. Its GW capture
    // is handled by the encounter prescription, not by Peters' bound equations.
    p.time_merger=q.time_merger=predicted;
    p.time_record=q.time_record=now;
    p.time_interrupt=q.time_interrupt=-(1+now);
    p.star.tphys=q.star.tphys=now*interaction.bse_manager.tscale;
    if (reached && interaction.stellar_evolution_write_flag) {
        interaction.fout_bse<<"BH_GW_terminal "<<std::setprecision(17)
            <<now<<" "<<p.id<<" "<<q.id
            <<" predicted="<<predicted<<" delay="<<lag
            <<" cutoff="<<cutoff<<" a="<<bin.semi<<" e="<<bin.ecc<<std::endl;
    }
    return reached;
}
} // namespace BHGW
