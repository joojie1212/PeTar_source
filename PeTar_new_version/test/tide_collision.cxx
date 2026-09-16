#include <iostream>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024
#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"

// Exercise the real BSE merger/tide dispatch using instantaneous contact.
// G=1 and the default BSE radius scale define consistent test input units.
int checkEncounter(double ecc, double factor, double phase) {
    ARInteraction interaction;
    IOParamsBSE io;
    interaction.bse_manager.initial(io, false);
    interaction.gravitational_constant=1;
    interaction.eps_sq=0;
    interaction.stellar_evolution_option=2;
    interaction.stellar_evolution_write_flag=false;
    interaction.tide.gravitational_constant=1;
    interaction.tide.speed_of_light=interaction.bse_manager.getSpeedOfLight();
    PtclHard p[2];
    for(int i=0;i<2;i++) {
        p[i].id=i+1;p[i].mass=1;p[i].star.initial(1);
        StarParameterOut out;interaction.bse_manager.evolveStar(p[i].star,out,0,true);
        p[i].radius=interaction.bse_manager.getMergerRadius(p[i].star);
        p[i].time_record=0;p[i].time_interrupt=1e10;
        p[i].setBinaryPairID(2-i);
        p[i].setBinaryInterruptState(BinaryInterruptState::none);
    }
    AR::BinaryTree<PtclHard> bin;
    bin.setMembers(p,p+1,0,1);
    bin.m1=bin.m2=1;bin.mass=2;
    bin.ecc=ecc;bin.semi=(p[0].radius+p[1].radius)*factor/(1-ecc);
    bin.incline=bin.rot_horizon=bin.rot_self=0;bin.ecca=phase;
    bin.slowdown.setSlowDownFactor(1);
    bin.calcParticles(1);bin.calcCenterOfMass();bin.calcOrbit(1);
    AR::InterruptBinary<PtclHard> intr;intr.time_now=0;intr.time_end=1e-6;
    const bool expected=(bin.r<p[0].radius+p[1].radius);
    interaction.modifyAndInterruptIter(intr,bin);
    bool merged=(p[0].mass==0 || p[1].mass==0);
    for (int i=0; i<2; i++) {
        if (!std::isfinite(p[i].mass)) return 4;
        for (int k=0; k<3; k++)
            if (!std::isfinite(p[i].pos[k]) || !std::isfinite(p[i].vel[k])) return 4;
    }
    if (expected && intr.status!=AR::InterruptStatus::merge) return 5;

    std::cout<<std::setprecision(16)<<"e="<<ecc<<" factor="<<factor<<" phase="<<phase<<" merged="<<merged<<" expected="<<expected<<"\n";
    return merged==expected?0:2;
}

int main() {
    int failures=0;
    for (double ecc: {0.5, 1.5}) {
        for (double factor: {0.5, 1.0-1e-12, 1.0, 1.0+1e-12, 10.0})
            failures += (checkEncounter(ecc, factor, 2.0)!=0);
        failures += (checkEncounter(ecc, 0.5, -2.0)!=0);
        // At peri-centre, factor is the instantaneous separation in units
        // of the summed stellar radius.
        failures += (checkEncounter(ecc, 0.5, 0.0)!=0);
        failures += (checkEncounter(ecc, 1.5, 0.0)!=0);
    }
    return failures ? 1 : 0;
}
