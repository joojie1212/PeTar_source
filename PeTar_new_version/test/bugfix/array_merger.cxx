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

struct MergeInjection: ARInteraction {
    double kappa=1.;
    int modifyAndInterruptIter(AR::InterruptBinary<PtclHard>& event,
                               AR::BinaryTree<PtclHard>& bin) {
        auto* a=bin.getMember(0); auto* b=bin.getMember(1);
        a->mass+=b->mass; b->mass=0.;
        a->pos=PS::F64vec(0.,0.,0.); a->vel=PS::F64vec(1.,2.,3.);
        bin.slowdown.setSlowDownFactor(kappa);
        event.adr=&bin; event.status=AR::InterruptStatus::merge;
        return 2;
    }
};
int main() {
    for(double kappa: {1.,8.}) {
        AR::TimeTransformedSymplecticManager<MergeInjection> manager;
        manager.time_error_max=1e-10; manager.energy_error_relative_max=1e-8;
        manager.time_step_min=1e-12; manager.slowdown_pert_ratio_ref=1e-6;
        manager.slowdown_timescale_max=100.; manager.step_count_max=10000;
        manager.interrupt_detection_option=1; manager.step.initialSymplecticCofficients(2);
        manager.interaction.eps_sq=0.; manager.interaction.gravitational_constant=1.;
        manager.interaction.kappa=kappa;
        AR::TimeTransformedSymplecticIntegrator<PtclHard,PtclHard,ARPerturber,
            MergeInjection,AR::Information<PtclHard,PtclHard>> system;
        system.manager=&manager;
        system.particles.setMode(COMM::ListMode::local); system.particles.reserveMem(2);
        for(int i=0;i<2;++i) {
            PtclHard p; p.id=i+1; p.mass=1.;
            p.pos=PS::F64vec(i-0.5,0.,0.); p.vel=PS::F64vec(0.,i-0.5,0.);
            p.changeover.setR(1.,2.); system.particles.addMember(p);
        }
        system.info.reserveMem(2); system.info.r_break_crit=2.;
        system.info.generateBinaryTree(system.particles,1.);
        system.perturber.r_neighbor_crit_sq=4.;
        system.reserveIntegratorMem(); system.initialIntegration(0.);
        system.info.ds=0.01;
        const auto event=system.integrateToTime(2.);
        assert(event.status==AR::InterruptStatus::merge && system.getTime()==2.);
        bool checked=false;
        for(int i=0;i<2;++i) if(system.particles[i].mass>0.) {
            for(int k=0;k<3;++k) assert(std::abs(system.particles[i].pos[k]-2.*(k+1)/kappa)<1e-12);
            checked=true;
        }
        assert(checked);
    }
    std::cout<<"ARRAY merger drift regressions passed\n";
}
