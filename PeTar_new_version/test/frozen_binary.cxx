#include <iostream>
#include <fstream>
#include <string>
#include <particle_simulator.hpp>
#define HARD_DEBUG_PRINT_FEQ 1024
#include "io.hpp"
#include "hard_assert.hpp"
#include "cluster_list.hpp"
#include "hard.hpp"
#include "soft_ptcl.hpp"
#include "static_variables.hpp"

#ifndef FROZEN_BINARY
#error "This regression must be compiled with FROZEN_BINARY"
#endif

int main() {
    ARInteraction interaction;
    interaction.gravitational_constant = 1.0;
    interaction.frozen_energy_factor = 10.0;
    interaction.frozen_kinetic_energy_ref = 1.0;
    interaction.frozen_min_interval = 10.0;
    interaction.frozen_perturbation_limit = 1.0e-6;
    interaction.frozen_radius_factor = 3.0;
    if (interaction.frozen_ecc_limit != 0.1) return 25;

    PtclHard p[2];
    for (int i=0; i<2; i++) {
        p[i].id = i+1;
        p[i].mass = 1.0;
        p[i].radius = 1.0e-4;
        p[i].star.initial(1.0);
        p[i].star.kw = 1;
        p[i].setBinaryInterruptState(BinaryInterruptState::none);
    }

    AR::BinaryTree<PtclHard> bin;
    bin.setMembers(p, p+1, 0, 1);
    bin.m1 = bin.m2 = 1.0;
    bin.mass = 2.0;
    bin.semi = 0.01;
    bin.ecc = 0.0;
    bin.ecca = 0.0;
    bin.incline = bin.rot_horizon = bin.rot_self = 0.0;
    bin.slowdown.pert_in = 1.0;
    bin.slowdown.pert_out = 1.0e-8;
    bin.calcParticles(1.0);
    bin.calcCenterOfMass();
    bin.calcOrbit(1.0);

    auto eligible = [&]() {
        return interaction.isBinaryFrozenEligible(bin, p[0], p[1], 0.0, 100.0);
    };
    if (!eligible()) return 1;

    p[0].time_record = p[1].time_record = 0.0;
    p[0].time_interrupt = p[1].time_interrupt = 100.0;
    AR::InterruptBinary<PtclHard> interrupt;
    interrupt.time_now = 1.0;
    interrupt.time_end = 2.0;
    if (interaction.modifyAndInterruptIter(interrupt, bin)!=0) return 2;
    if (!interaction.isFrozenPair(p[0], p[1])) return 2;
    if (p[0].getBinaryPairID()!=p[1].id || p[1].getBinaryPairID()!=p[0].id)
        return 3;
    p[0].setBinaryInterruptState(BinaryInterruptState::form);
    if (p[0].isBinaryFrozen() || p[0].getBinaryPairID()!=p[1].id) return 4;
    p[0].setBinaryInterruptState(BinaryInterruptState::none);
    interaction.setFrozenPair(p[0], p[1], true);

    const double semi_save = bin.semi;
    bin.semi = 1.0;
    if (eligible()) return 5;
    bin.semi = semi_save;

    const double pert_save = bin.slowdown.pert_out;
    bin.slowdown.pert_out = 1.0e-4;
    if (eligible()) return 6;
    bin.slowdown.pert_out = pert_save;

    const double ecc_save = bin.ecc;
    bin.ecc = 0.99;
    if (eligible()) return 7;
    bin.ecc = ecc_save;

    p[0].setBinaryInterruptState(BinaryInterruptState::tide);
    if (eligible()) return 8;
    p[0].setBinaryInterruptState(BinaryInterruptState::none);

    // A short BSE interval is unsafe for an eccentric pair, but an already
    // circular, detached pair may override BSE's repeated zero interval.
    bin.ecc = 0.21;
    if (interaction.isBinaryFrozenEligible(bin, p[0], p[1], 95.0, 100.0)) return 9;
    bin.ecc = 0.0;
    if (!interaction.isBinaryFrozenEligible(bin, p[0], p[1], 100.0, 100.0))
        return 10;

    interaction.setFrozenPair(p[0], p[1], false);
    if (p[0].isBinaryFrozen() || p[1].isBinaryFrozen()) return 11;

    p[0].time_record = p[1].time_record = 99.0;
    p[0].time_interrupt = p[1].time_interrupt = 100.0;
    interrupt.time_now = 100.0;
    interrupt.time_end = 101.0;
    if (interaction.modifyAndInterruptIter(interrupt, bin)!=0) return 12;
    if (!interaction.isFrozenPair(p[0], p[1])) return 13;
    if (p[0].time_interrupt<110.0 || p[1].time_interrupt<110.0) return 14;
    interaction.setFrozenPair(p[0], p[1], false);

    // Tight circular binaries bypass the configurable clearance, but contact,
    // active BSE events and strong perturbations must still prevent takeover.
    const double saved_r = bin.r;
    bin.semi = bin.r = 1.9*(p[0].radius+p[1].radius);
    bin.ecc = 0.0;
    if (!eligible()) return 17;
    bin.ecc = 0.09;
    if (!eligible()) return 26;
    bin.ecc = 0.099999;
    if (!eligible()) return 27;
    if (!interaction.isBinaryFrozenEligible(bin, p[0], p[1], 100.0, 100.0))
        return 28;
    bin.ecc = 0.1;
    if (eligible()) return 29;
    const char* reason = nullptr;
    interaction.frozen_ecc_limit = 0.2;
    if (interaction.isBinaryFrozenEligible(bin, p[0], p[1], 0.0, 100.0, &reason)
        || std::string(reason)!="eccentricity") return 30;
    interaction.frozen_ecc_limit = 0.1;
    bin.ecc = 2.0*interaction.frozen_ecc_limit;
    if (eligible()) return 18;
    bin.ecc = 0.0;
    bin.slowdown.pert_out = 1.0e-4;
    if (eligible()) return 19;
    bin.slowdown.pert_out = pert_save;
    p[0].setBinaryInterruptState(BinaryInterruptState::tide);
    if (eligible()) return 20;
    p[0].setBinaryInterruptState(BinaryInterruptState::none);
    bin.r = p[0].radius+p[1].radius;
    if (eligible()) return 21;
    bin.r = saved_r;
    bin.semi = p[0].radius+p[1].radius;
    if (eligible()) return 22;
    // Isolate the hard-coded strict upper boundary from the ordinary policy.
    interaction.frozen_radius_factor = 12.0;
    bin.semi = bin.r = 9.99*(p[0].radius+p[1].radius);
    if (!eligible()) return 23;
    bin.semi = bin.r = 10.0*(p[0].radius+p[1].radius);
    if (eligible()) return 24;
    interaction.frozen_radius_factor = 3.0;
    bin.semi = semi_save;
    bin.r = saved_r;

    // The SDAR bypass must preserve a bound Kepler orbit while advancing its
    // phase analytically.  One full period should recover the initial state.
    bin.calcOrbit(1.0);
    const double period = bin.period;
    const PS::F64vec pos0[2] = {p[0].pos, p[1].pos};
    const PS::F64vec vel0[2] = {p[0].vel, p[1].vel};
    HardIntegrator::advanceFrozenKeplerOrbit(bin, period, 1.0);
    for (int k=0; k<3; k++) {
        if (std::fabs(p[0].pos[k]-pos0[0][k])>1.0e-11
            || std::fabs(p[1].pos[k]-pos0[1][k])>1.0e-11
            || std::fabs(p[0].vel[k]-vel0[0][k])>1.0e-9
            || std::fabs(p[1].vel[k]-vel0[1][k])>1.0e-9) return 15;
    }

    HardIntegrator::advanceFrozenKeplerOrbit(bin, 0.25*period, 1.0);
    const PS::F64vec dr = p[1].pos-p[0].pos;
    const PS::F64vec dv = p[1].vel-p[0].vel;
    const double specific_energy = 0.5*(dv*dv)-2.0/std::sqrt(dr*dr);
    if (std::fabs(specific_energy+1.0/bin.semi)>1.0e-9) return 16;

    // Verify the missing tidal thaw path logs exactly once before clearing flags.
    const char* log_path = "build/frozen_reason_test.log";
    interaction.fout_bse.open(log_path);
    if (!interaction.fout_bse) return 31;
    interaction.stellar_evolution_write_flag = true;
    interaction.setFrozenPair(p[0], p[1], true);
    interaction.thawFrozenPair(p[0], p[1], 100.0, "tidal_orbit_change");
    interaction.thawFrozenPair(p[0], p[1], 100.0, "bse_update");
    interaction.fout_bse.close();
    interaction.stellar_evolution_write_flag = false;
    std::ifstream log(log_path);
    std::string line, extra;
    if (!std::getline(log, line) || line.find("reason=tidal_orbit_change")==std::string::npos
        || std::getline(log, extra) || interaction.isFrozenPair(p[0], p[1])) return 32;
    std::cout << "frozen-binary policy checks passed\n";
    return 0;
}
