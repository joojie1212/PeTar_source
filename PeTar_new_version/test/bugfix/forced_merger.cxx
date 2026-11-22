#include <iostream>
#include <cassert>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024
#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"

int main() {
    ARInteraction interaction;
    interaction.gravitational_constant=0.00449830997959438;
    interaction.bse_manager.vscale=0.977813107686401;
    interaction.stellar_evolution_option=1;
    interaction.stellar_evolution_write_flag=false;
    PtclHard p1, p2;
    p1.id=1; p2.id=2;
    for (auto* p: {&p1, &p2}) {
        p->mass=30.; p->dm=0.; p->radius=0.;
        p->star.initial(30.,14);
        p->time_record=1.; p->time_interrupt=10.;
        p->binary_state=0;
    }
    AR::BinaryTree<PtclHard> bin;
    bin.setMembers(&p1,&p2,0,1);
    bin.mass=60.;
    const double tight=5e-9, wide=5e-5, ecc=0.2;
    auto setState = [&](double semi, double speed_factor) {
        const double r=semi*(1.-ecc);
        const double v=speed_factor*std::sqrt(
            interaction.gravitational_constant*60.*(1.+ecc)/r);
        p1.pos=PS::F64vec(-r/2.,0.,0.);
        p2.pos=PS::F64vec(r/2.,0.,0.);
        p1.vel=PS::F64vec(0.,-v/2.,0.);
        p2.vel=PS::F64vec(0.,v/2.,0.);
    };
    auto staleCache = [&](double semi) { bin.semi=semi; bin.ecc=ecc; };
    AR::InterruptBinary<PtclHard> event;

    // A cached compact orbit must not rescue a currently unbound/wide pair.
    for (bool unbound: {true, false}) {
        setState(unbound ? tight : wide, unbound ? 2. : 1.);
        staleCache(tight); event.clear();
        assert(!interaction.tryForcedMerger(p1,p2,bin,1e-6,event));
        assert(unbound ? bin.semi<0. : std::abs(bin.semi/wide-1.)<1e-12);
        assert(event.status==AR::InterruptStatus::none);
    }

    // A stale wide orbit must not hide a compact candidate. The updated cache
    // is also what the integrator uses to rank eligible leaves by semi.
    setState(tight,1.); staleCache(wide); event.clear();
    assert(interaction.tryForcedMerger(p1,p2,bin,1e-6,event));
    assert(std::abs(bin.semi/tight-1.)<1e-12);
    assert(std::abs(bin.ecc-ecc)<1e-12);
    assert(event.status==AR::InterruptStatus::forced_merge);
    assert(event.forced_id[0]==1 && event.forced_id[1]==2);
    assert(std::abs(event.theoretical_merger_time/
        interaction.petersMergerLifetime(tight,ecc,30.,30.)-1.)<1e-11);

    // An eligible request can become stale before consumption. Revalidate it
    // without losing the interrupt clock or merging either particle.
    for (bool unbound: {true, false}) {
        setState(tight,1.); event.clear();
        assert(interaction.tryForcedMerger(p1,p2,bin,1e-6,event));
        event.adr=&bin; event.time_now=1.; event.time_end=1.+1e-6;
        event.theoretical_merger_time+=event.time_now;
        setState(unbound ? tight : wide, unbound ? 2. : 1.);
        assert(interaction.modifyAndInterruptIter(event,bin)==0);
        assert(p1.mass==30. && p2.mass==30.);
        assert(event.status==AR::InterruptStatus::none && !event.forced);
        assert(event.time_now==1. && event.time_end==1.+1e-6);
        assert(event.forced_id[0]==-1 && event.forced_id[1]==-1);
        assert(event.theoretical_merger_time==-1.);
    }

    interaction.stellar_evolution_option=0;
    setState(tight,1.); event.clear();
    assert(!interaction.tryForcedMerger(p1,p2,bin,1e-6,event));

    // A still-eligible request must continue through the real BSE merge path.
    IOParamsBSE params;
    params.rscale.value=4.4334e7;
    params.vscale.value=0.977813107686401;
    interaction.bse_manager.initial(params,false);
    interaction.stellar_evolution_option=1;
    p1.star.tphys=p2.star.tphys=1.;
    event.clear();
    assert(interaction.tryForcedMerger(p1,p2,bin,1e-6,event));
    event.adr=&bin; event.time_now=1.; event.time_end=1.+1e-6;
    // Corrupt both cached orbit and diagnostic time after issuing the request.
    staleCache(wide); event.theoretical_merger_time=123.;
    assert(interaction.modifyAndInterruptIter(event,bin)==2);
    assert(event.status==AR::InterruptStatus::merge && event.forced);
    assert((p1.mass>0.) != (p2.mass>0.));
    assert(std::abs(bin.semi/tight-1.)<1e-12);
    assert(std::abs(event.theoretical_merger_time-
        (1.+interaction.petersMergerLifetime(tight,ecc,30.,30.)))<1e-12);
    std::cout<<"forced merger current-orbit regressions passed\n";
}
