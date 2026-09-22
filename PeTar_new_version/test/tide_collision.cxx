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
int checkEncounter(double ecc, double factor, double phase, double debug_radius_rsun=0.0) {
    ARInteraction interaction;
    IOParamsBSE io;
    interaction.bse_manager.initial(io, false);
    interaction.gravitational_constant=1;
    interaction.eps_sq=0;
    interaction.stellar_evolution_option=2;
    interaction.stellar_evolution_write_flag=false;
    interaction.tide.gravitational_constant=1;
    interaction.tide.speed_of_light=interaction.bse_manager.getSpeedOfLight();
    interaction.debug_lessmerger_rsun=debug_radius_rsun;
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
    const bool expected=(bin.r<interaction.mergerCheckRadius(p[0],p[1]));
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

// The extra prescription stops at e=0.9, including after capture.
int checkTideScope() {
    TwoBodyTide tide;
    tide.gravitational_constant=1;
    PtclHard p[2];
    AR::BinaryTree<PtclHard> bin;
    bin.setMembers(p,p+1,0,1);
    bin.m1=bin.m2=1;bin.mass=2;
    for (double ecc: {0.5, std::nextafter(0.9,0.0), 0.9}) {
        bin.semi=3/(1-ecc);bin.ecc=ecc;
        const double a=bin.semi;
        if (tide.evolveOrbitDynamicalTide(bin,1.,1.,3.,3.)!=0
            || bin.semi!=a || bin.ecc!=ecc) return 1;
    }
    for (double ecc: {0.900001, 0.95}) {
        bin.semi=3/(1-ecc);bin.ecc=ecc;
        const double a=bin.semi,pold=a*(1-ecc*ecc);
        const double loss=tide.evolveOrbitDynamicalTide(bin,1.,1.,3.,3.);
        if (!(loss>0 && bin.semi<a && bin.ecc>=0.9 && bin.ecc<ecc)) return 2;
        if (std::abs(bin.semi*(1-bin.ecc*bin.ecc)-pold)>1e-12*pold) return 3;
        if (std::abs(loss-(0.5/bin.semi-0.5/a))>1e-12) return 4;
    }
    bin.ecc=1.00001;bin.semi=3/(1-bin.ecc);
    const double loss=tide.evolveOrbitDynamicalTide(bin,1.,1.,3.,3.);
    if (!(loss>0 && bin.semi>0 && bin.ecc>0.9 && bin.ecc<1)) return 5;
    if (!(tide.evolveOrbitDynamicalTide(bin,1.,1.,3.,3.)>0)) return 6;
    // A strong final loss is limited to the boundary, with finite orbital elements.
    bin.ecc=0.900001;bin.semi=3/(1-bin.ecc);
    if (!(tide.evolveOrbitDynamicalTide(bin,1.,1.,3.,3.)>0)
        || bin.ecc!=0.9 || !std::isfinite(bin.semi)) return 7;
    const double a=bin.semi;
    if (tide.evolveOrbitDynamicalTide(bin,1.,1.,3.,3.)!=0 || bin.semi!=a) return 8;
    std::cout<<"tide scope: low-e skipped, high-e bound loss, capture, e=0.9 handoff\n";
    return 0;
}

struct BSETideResult {
    double ecc, period, spin, mass, time;
    int error, tflag;
};

BSETideResult evolveBSETide(double ecc, double cutoff, int tflag) {
    BSEManager manager;
    IOParamsBSE io;
    io.tflag.value=tflag;
    manager.initial(io,false);
    StarParameter stars[2];
    StarParameterOut out[2];
    for (int i=0;i<2;i++) {
        stars[i].initial(1);
        manager.evolveStar(stars[i],out[i],0,true);
    }
    double semi=8/(1-ecc); // detached periapsis, Rsun (default rscale=1)
    double period=std::sqrt(std::pow(semi/215.0954,3)/2)*365.25/manager.year_to_day;
    BinaryEvent event;
    const int error=manager.evolveBinary(stars[0],stars[1],out[0],out[1],
                                        semi,period,ecc,event,0,0.01,cutoff);
    return {ecc,period,stars[0].ospin,stars[0].mt,stars[0].tphys,error,flags_.tflag};
}

bool sameBSE(const BSETideResult& a, const BSETideResult& b) {
    return std::abs(a.ecc-b.ecc)<1e-12
        && std::abs(a.period-b.period)<1e-12
        && std::abs(a.spin-b.spin)<1e-12
        && std::abs(a.mass-b.mass)<1e-12 && a.time==b.time;
}

int checkBSETideSwitch() {
    const auto high=evolveBSETide(0.95,0.9,1);
    const auto off=evolveBSETide(0.95,-1,0);
    const auto original=evolveBSETide(0.95,-1,1);
    std::cout<<std::setprecision(17)<<"BSE high-e: hybrid="<<high.ecc<<" off="<<off.ecc
             <<" unrestricted="<<original.ecc<<"\n";
    if (high.error<0 || high.tflag!=1 || high.time!=0.01
        || !sameBSE(high,off) || sameBSE(high,original)) return 1;
    for (double ecc: {0.89,0.9}) {
        const auto hybrid=evolveBSETide(ecc,0.9,1);
        const auto baseline=evolveBSETide(ecc,-1,1);
        if (hybrid.error<0 || hybrid.tflag!=1 || !sameBSE(hybrid,baseline)) return 2;
    }
    const auto disabled=evolveBSETide(0.89,0.9,0);
    if (disabled.tflag!=0 || !sameBSE(disabled,evolveBSETide(0.89,-1,0))) return 3;
    std::cout<<"BSE switch: high-e tides off, stellar time advances, low-e tides retained, tflag unchanged\n";
    return 0;
}

int checkDynamicalDispatch() {
    ARInteraction interaction;
    IOParamsBSE io;
    interaction.bse_manager.initial(io,false);
    interaction.gravitational_constant=1;
    interaction.eps_sq=0;
    interaction.stellar_evolution_option=2;
    interaction.stellar_evolution_write_flag=false;
    interaction.tide.gravitational_constant=1;
    interaction.tide.speed_of_light=interaction.bse_manager.getSpeedOfLight();
    for (double ecc: {0.89,0.95}) {
        PtclHard p[2];
        for (int i=0;i<2;i++) {
            p[i].id=i+1;p[i].mass=1;p[i].star.initial(1);
            StarParameterOut out;
            interaction.bse_manager.evolveStar(p[i].star,out,0,true);
            p[i].radius=interaction.bse_manager.getMergerRadius(p[i].star);
            p[i].time_record=0;p[i].time_interrupt=1e10;
            p[i].setBinaryPairID(2-i);
            p[i].setBinaryInterruptState(BinaryInterruptState::none);
        }
        AR::BinaryTree<PtclHard> bin;
        bin.setMembers(p,p+1,0,1);
        bin.m1=bin.m2=1;bin.mass=2;
        bin.ecc=ecc;bin.semi=1.6*(p[0].radius+p[1].radius)/(1-ecc);
        bin.incline=bin.rot_horizon=bin.rot_self=0;bin.ecca=2;
        bin.slowdown.setSlowDownFactor(1);
        bin.calcParticles(1);bin.calcCenterOfMass();bin.calcOrbit(1);
        const double a=bin.semi;
        AR::InterruptBinary<PtclHard> intr;
        intr.time_now=0;intr.time_end=1e-6;
        interaction.modifyAndInterruptIter(intr,bin);
        if ((bin.semi<a)!=(ecc>0.9)) return 1;
        if (p[0].mass<=0 || p[1].mass<=0) return 2;
        if (ecc>0.9) {
            if (p[0].getBinaryInterruptState()!=BinaryInterruptState::tide) return 3;
            const double after=bin.semi;
            interaction.modifyAndInterruptIter(intr,bin);
            if (std::abs(bin.semi-after)>1e-12*after) return 4;
        }
    }
    std::cout<<"AR dispatch: high-e bound tides only, no duplicate outgoing-passage loss\n";
    return 0;
}

int main() {
    int failures=0;
    const int tide_scope=checkTideScope(), bse_switch=checkBSETideSwitch();
    std::cout<<"scope_result="<<tide_scope<<" bse_result="<<bse_switch<<"\n";
    failures += (tide_scope!=0) + (bse_switch!=0);
    const int dispatch=checkDynamicalDispatch();
    std::cout<<"dispatch_result="<<dispatch<<"\n";
    failures += (dispatch!=0);
    for (double ecc: {0.5, 0.9, 0.95, 1.5}) {
        for (double factor: {0.5, 1.0-1e-12, 1.0, 1.0+1e-12, 10.0})
            failures += (checkEncounter(ecc, factor, 2.0)!=0);
        failures += (checkEncounter(ecc, 0.5, -2.0)!=0);
        // At peri-centre, factor is the instantaneous separation in units
        // of the summed stellar radius.
        failures += (checkEncounter(ecc, 0.5, 0.0)!=0);
        failures += (checkEncounter(ecc, 1.5, 0.0)!=0);
    }
    // A physical surface overlap is not a merger in the debug mode when the
    // separation stays above its explicit (much smaller) pair threshold.
    failures += (checkEncounter(0.5, 0.5, 0.0, 1.0e-6)!=0);
    failures += (checkEncounter(1.5, 0.5, 0.0, 1.0e-6)!=0);
    return failures ? 1 : 0;
}
