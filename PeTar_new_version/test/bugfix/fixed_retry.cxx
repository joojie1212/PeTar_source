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
struct CountingInteraction: ARInteraction {
    double last=-1.; int calls=0;
    AR::FixStepOption* policy=nullptr; bool saw_adaptive=false;
    int modifyAndInterruptIter(AR::InterruptBinary<PtclHard>& event,
                               AR::BinaryTree<PtclHard>& bin) {
        // A rejected gravitational trial must not replay an accepted callback.
        assert(event.time_now>last); last=event.time_now; ++calls;
        saw_adaptive |= policy && *policy==AR::FixStepOption::none;
        return 0;
    }
};
int main() {
    for(double ds: {1e-5, 0.1}) {
        AR::TimeTransformedSymplecticManager<CountingInteraction> manager;
        manager.time_error_max=1e-10; manager.energy_error_relative_max=1e-9;
        manager.time_step_min=1e-14; manager.slowdown_pert_ratio_ref=1e-6;
        manager.slowdown_timescale_max=100.; manager.step_count_max=1000000;
        manager.interrupt_detection_option=1; manager.step.initialSymplecticCofficients(2);
        manager.interaction.eps_sq=0.; manager.interaction.gravitational_constant=1.;
        AR::TimeTransformedSymplecticIntegrator<PtclHard,PtclHard,ARPerturber,
            CountingInteraction,AR::Information<PtclHard,PtclHard>> system;
        system.manager=&manager;
        system.particles.setMode(COMM::ListMode::local); system.particles.reserveMem(3);
        for(int i=0;i<3;++i) {
            PtclHard p; p.id=i+1; p.mass=1.;
            p.pos=PS::F64vec(i==2?10.:i-0.5,0.,0.);
            p.vel=PS::F64vec(0.,i==2?0.3:(i?0.7071067811865475:-0.7071067811865475),0.);
            p.changeover.setR(100.,200.); system.particles.addMember(p);
        }
        system.info.reserveMem(3); system.info.r_break_crit=100.;
        system.info.generateBinaryTree(system.particles,1.);
        system.perturber.r_neighbor_crit_sq=40000.;
        system.reserveIntegratorMem(); system.initialIntegration(0.);
        system.info.ds=ds; system.info.fix_step_option=AR::FixStepOption::later;
        manager.interaction.policy=&system.info.fix_step_option;
        system.integrateToTime(.01);
#ifdef EXPECT_FIXED_RETRY
        if (ds==1e-5) assert(manager.interaction.saw_adaptive);
#else
        if (ds==1e-5) assert(!manager.interaction.saw_adaptive);
#endif
        assert(system.info.fix_step_option==AR::FixStepOption::later);
        assert(std::abs(system.getTime()-.01)<1e-9);
        std::cout << "ds="<<ds<<" calls="<<manager.interaction.calls
                  <<" trials="<<system.profile.step_count_sum<<"\n";
        if(ds==0.1) assert(system.profile.step_count_sum>static_cast<unsigned>(manager.interaction.calls));
    }
}
