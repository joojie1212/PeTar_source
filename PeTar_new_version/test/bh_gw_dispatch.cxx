#include <iostream>
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
    IOParamsBSE io;
    interaction.bse_manager.initial(io,false);
    interaction.gravitational_constant=1;
    interaction.eps_sq=0;
    interaction.stellar_evolution_option=2;
    interaction.stellar_evolution_write_flag=false;
    interaction.tide.gravitational_constant=1;
    interaction.tide.speed_of_light=interaction.bse_manager.getSpeedOfLight();
    PtclHard p[2];
    for (int i=0;i<2;++i) {
        p[i].id=i+1;p[i].mass=30;
        p[i].star.initial(30);
        p[i].star.kw=14;p[i].star.mt=p[i].star.mc=30;
        p[i].star.r=0.0001272;
        p[i].radius=interaction.bse_manager.getMergerRadius(p[i].star);
        p[i].time_record=0;p[i].time_interrupt=1e10;
        p[i].time_merger=0; // expired LEGACY deadline must not trigger merger
        p[i].spin=PS::F64vec(0.1,0.2,0.3);
        p[i].setBinaryPairID(2-i);
        p[i].setBinaryInterruptState(BinaryInterruptState::none);
    }
    AR::BinaryTree<PtclHard> bin;
    bin.setMembers(p,p+1,0,1);
    bin.m1=bin.m2=30;bin.mass=60;
    bin.incline=0.3;bin.rot_horizon=0.4;bin.rot_self=0.5;
    bin.slowdown.setSlowDownFactor(1);
    auto place=[&](double a,double e,double phase) {
        bin.semi=a;bin.ecc=e;bin.ecca=phase;
        bin.calcParticles(1);
        for(int i=0;i<2;++i) {
            p[i].pos+=PS::F64vec(0.02,-0.03,0.01);
            p[i].vel+=PS::F64vec(1,2,3);
        }
        bin.calcCenterOfMass();bin.calcOrbit(1);
    };
    AR::InterruptBinary<PtclHard> intr;
    auto update=[&](double t) {
        intr.clear();intr.time_now=t;intr.time_end=t+1;
        interaction.modifyAndInterruptIter(intr,bin);
    };
    place(1e-5,0.5,2.5);
    update(10);
    assert(p[0].mass==30 && p[1].mass==30);
    assert(std::abs(bin.semi/1e-5-1)<1e-8); // no replay of old SSE interval
    const double early=p[0].time_merger;
    place(2e-5,0.5,2.5); // external perturbation widens the orbit
    update(10);
    assert(p[0].time_merger>early);
    const double wide=p[0].time_merger;
    update(10+(wide-10)*0.01);
    assert(bin.semi<2e-5 && bin.ecc<0.5);
    for(int k=0;k<3;++k) {
        assert(std::abs((p[0].vel[k]+p[1].vel[k])/2-(k+1))<1e-8);
    }
    const double a=bin.semi;
    update(intr.time_now);
    assert(std::abs(bin.semi/a-1)<1e-8); // repeated callback is idempotent
    place(-2e-5,1.1,-0.5);
    update(intr.time_now+1e-6);
    assert(p[0].time_merger<0 && p[1].time_merger<0);
    assert(p[0].mass==30 && p[1].mass==30);
    place(1e-5,0.2,1);
    update(intr.time_now);
    // Finish the bound inspiral at the first callback after its analytic endpoint.
    const double endpoint=p[0].time_merger;
    update(endpoint+1e-5);
    assert(p[0].mass>0 && p[0].mass<60 && p[1].mass==0);
    assert(intr.status==AR::InterruptStatus::merge);
    const auto fit=BHMerger::getRemnant();
    for (int k=0;k<3;++k) {
        const double vcm=k+1;
        const double kick=fit.vkick_nor[k]/interaction.bse_manager.vscale;
        assert(std::abs(p[0].vel[k]-vcm-kick)<1e-6);
        const double change=p[0].mass*p[0].vel[k]-60*vcm;
        const double expected=-(60-p[0].mass)*vcm+p[0].mass*kick;
        assert(std::abs(change-expected)<1e-4);
    }
    assert(std::abs(p[0].pos.x-0.02)<1e-12);
    assert(std::abs(p[0].pos.y+0.03)<1e-12);
    assert(std::abs(p[0].pos.z-0.01)<1e-12);
    std::cout<<"BH GW dispatch: legacy clock, postponement, dissipation, CM, "
        "idempotence, unbinding, terminal merger and recoil momentum passed\n";
}
