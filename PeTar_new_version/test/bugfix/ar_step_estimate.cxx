#include <iostream>
#include <cassert>
#include <cfenv>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024
#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"

void testKeplerCandidates() {
    AR::Information<PtclHard,PtclHard> info;
    PtclHard p[3];
    AR::BinaryTree<PtclHard> inner, root;
    inner.setMembers(&p[0], &p[1], 0, 1);
    root.setMembers(static_cast<PtclHard*>(&inner), &p[2], -1, 2);
    inner.m1=inner.m2=1.; inner.semi=1.;
    root.m1=2.; root.m2=1.; root.semi=100.;
    const double scale=std::pow(0.1,1./6.);

    // The outer candidate has a much larger step. Moving the inner metric
    // across the minimum should select it throughout the near-tie window.
    for (double delta: {-0.005, 0.005, 0.02}) {
        inner.slowdown.setSlowDownFactor(std::sqrt(2./3.)*1000.*(1.+delta));
        const double got=info.calcDsKeplerBinaryTree(root,6,1.,1.);
        const double expected=info.calcDsElliptic(delta<=0.01 ? inner : root,1.)*scale;
        assert(std::abs(got-expected)<1e-12);
    }

    // A single elliptic binary retains the original estimate.
    assert(std::abs(info.calcDsKeplerBinaryTree(inner,6,1.,1.)
                    - info.calcDsElliptic(inner,1.)*scale)<1e-12);

    // Exercise both a hyperbolic pair and an entirely hyperbolic tree. No
    // elliptic candidate exists, so its maximum-value sentinel must not be scaled.
    inner.semi=-1.; root.semi=-100.;
    const double expected=info.calcDsHyperbolic(inner,1.);
    for (auto* bin: {&inner, &root}) {
        std::feclearexcept(FE_ALL_EXCEPT);
        const double got=info.calcDsKeplerBinaryTree(*bin,6,1.,1.);
        assert(std::fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)==0);
        assert(std::abs(got-expected)<1e-12);
    }

    // Even a finite elliptic metric can overflow when enlarged by 1%.
    inner.semi=1e100;
    const double elliptic_ds=info.calcDsElliptic(inner,1.)*scale;
    const double unscaled_metric=elliptic_ds*(2.*inner.semi);
    inner.slowdown.setSlowDownFactor((0.995*NUMERIC_FLOAT_MAX)/unscaled_metric);
    std::feclearexcept(FE_ALL_EXCEPT);
    const double got=info.calcDsKeplerBinaryTree(inner,6,1.,1.);
    assert(std::fetestexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO)==0);
    assert(std::abs(got/elliptic_ds-1.)<1e-12);
}

struct EntryStepInteraction: ARInteraction {
    Float* step=nullptr;
    Float first_step=-1.;
    int modifyAndInterruptIter(AR::InterruptBinary<PtclHard>&,
                               AR::BinaryTree<PtclHard>&) {
        if (first_step<0.) first_step=*step;
        return 0;
    }
};

void testEntryStepPolicy() {
    for (auto policy: {AR::FixStepOption::always, AR::FixStepOption::later,
                       AR::FixStepOption::none}) {
        for (double ds: {0.001, 1.}) {
            AR::TimeTransformedSymplecticManager<EntryStepInteraction> manager;
            manager.time_error_max=1e-10; manager.energy_error_relative_max=1e-9;
            manager.time_step_min=1e-14; manager.slowdown_pert_ratio_ref=1e-6;
            manager.slowdown_timescale_max=100.; manager.step_count_max=1000000;
            manager.interrupt_detection_option=1;
            manager.step.initialSymplecticCofficients(2);
            manager.interaction.eps_sq=0.; manager.interaction.gravitational_constant=1.;
            AR::TimeTransformedSymplecticIntegrator<PtclHard,PtclHard,ARPerturber,
                EntryStepInteraction,AR::Information<PtclHard,PtclHard>> system;
            system.manager=&manager;
            system.particles.setMode(COMM::ListMode::local); system.particles.reserveMem(2);
            for (int i=0;i<2;++i) {
                PtclHard p; p.id=i+1; p.mass=1.;
                p.pos=PS::F64vec(i-0.5,0.,0.);
                p.vel=PS::F64vec(0.,i?0.7071067811865475:-0.7071067811865475,0.);
                p.changeover.setR(100.,200.); system.particles.addMember(p);
            }
            system.info.reserveMem(2); system.info.r_break_crit=100.;
            system.info.generateBinaryTree(system.particles,1.);
            system.perturber.r_neighbor_crit_sq=40000.;
            system.reserveIntegratorMem(); system.initialIntegration(0.);
            system.info.ds=ds; system.info.fix_step_option=policy;
            manager.interaction.step=&system.info.ds;
            system.integrateToTime(.01);
            assert(std::abs(system.getTime()-.01)<1e-9);
            if (policy==AR::FixStepOption::always || ds==0.001)
                assert(manager.interaction.first_step==ds);
            else
                assert(manager.interaction.first_step>0. && manager.interaction.first_step<0.5*ds);
            if (policy==AR::FixStepOption::always) {
                assert(system.info.ds==ds);
                assert(system.info.fix_step_option==policy);
            }
        }
    }
}

int main() {
    testKeplerCandidates();
    testEntryStepPolicy();
    std::cout << "AR step estimate regressions passed\n";
}
