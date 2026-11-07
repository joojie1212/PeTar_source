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


#include <cassert>
#include <cstdio>
#include <sstream>
#include <vector>
#include <unistd.h>

struct TestOrbit {
    double semi=1.0, ecc=0.5, period=1.0;
    void calcOrbit(double) {}
};

int legacyLogTest() {
    const std::string path = "/tmp/petar-change-" + std::to_string(getpid());
    ARInteraction interaction;
    interaction.gravitational_constant = 1.0;
    interaction.bse_manager.vscale = 1.0;
    interaction.change_output.openChange(path, false, true, 17);
    PtclHard p[2];
    for (int i=0; i<2; ++i) {
        p[i].id=i+1;
        p[i].mass=i+1;
        p[i].pos=PS::F64vec(i, 2*i, 3*i);
        p[i].vel=PS::F64vec(4*i, 5*i, 6*i);
    }
    TestOrbit bin;
    // Initial prediction is not a change event, even if requested by caller.
    interaction.refreshBHMergerPrediction(bin, p[0], p[1], 1.0, true);
    const double old_deadline=p[0].time_merger;
    assert(!interaction.isBHMergerOrbitClearlyPerturbed(bin,p[0],p[1]));
    bin.semi=0.5;
    assert(interaction.isBHMergerOrbitClearlyPerturbed(bin,p[0],p[1]));
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],2.0,true);
    assert(p[0].merger_record_time==2.0);
    PtclHard copied;
    copied.DataCopy(p[0]);
    assert(copied.merger_record_time==2.0);
    // A replacement companion must not produce a false before/after pair.
    p[1].id=3;
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],3.0,true);
    interaction.change_output.close();
    interaction.change_output.openChange(path,true,true,17);
    bin.ecc=0.9;
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],4.0,true);
    // An unbound orbit records the invalidated deadline as -1.
    bin.semi=-1.0;
    interaction.refreshBHMergerPrediction(bin,p[0],p[1],5.0,true);
    assert(p[0].merger_record_time==-1.0);
    interaction.change_output.close();
    std::ifstream input(path);
    std::string line;
    std::getline(input,line);
    assert(line.find("time_old time_new id1 id2")==0);
    int rows=0;
    while (std::getline(input,line)) {
        std::istringstream stream(line);
        std::vector<double> values;
        double value;
        while(stream>>value) values.push_back(value);
        assert(values.size()==19);
        if (rows==0) {
            assert(values[0]==1.0 && values[1]==2.0);
            assert(values[2]==1 && values[3]==2);
            assert(values[4]==1 && values[5]==2);
            assert(values[6]==1.0 && values[7]==0.5);
            assert(values[8]==0.5 && values[9]==0.5);
            assert(values[10]==old_deadline);
            assert(values[11]<old_deadline);
        }
        if(rows==1) assert(values[0]==3 && values[1]==4 && values[3]==3);
        if(rows==2) assert(values[8]==-1 && values[11]==-1);
        for(int k=0;k<6;++k) assert(values[12+k]==-(k+1));
        ++rows;
    }
    assert(rows==3);
    std::remove(path.c_str());
    std::cout << "orbit change regression passed\n";
    return 0;
}


int countRows(const std::string& path) {
    std::ifstream f(path); std::string s; int n=-1;
    while (std::getline(f,s)) ++n;
    return n;
}

void clockSeparationTest() {
    const std::string path="/tmp/petar-clock-"+std::to_string(getpid());
    ARInteraction ar; ar.gravitational_constant=1; ar.bse_manager.vscale=1;
    ar.change_output.openChange(path,false,true,17);
    PtclHard p[2];
    for(int i=0;i<2;++i) {p[i].id=i+1;p[i].mass=30;p[i].pos=0;p[i].vel=0;}
    TestOrbit bin; bin.ecc=.9999;
    ar.checkBHMergerDynamicalChange(bin,p[0],p[1],0);
    const double deadline=p[0].time_merger;
    // Repeated tiny internal updates, all below the former absolute-e cutoff,
    // must neither leak into change nor restart an expired clock.
    for(int i=1;i<=2000;++i) {
        bin.ecc-=1e-7;bin.semi*=.99999;
        const double now=deadline*1.1*i/2000.;
        ar.syncBHMergerInternalChange(bin,p[0],p[1],now);
        ar.checkBHMergerDynamicalChange(bin,p[0],p[1],now);
        assert(p[0].time_merger==deadline && p[1].time_merger==deadline);
        assert(p[0].merger_dyn_log_tgw==0.);
    }
    assert(countRows(path)==0);
    assert(p[0].time_merger<deadline*1.1);
    // A tiny *dynamical* eccentricity increase at high e is significant.
    bin.ecc=.9999;bin.semi=1;
    ar.syncBHMergerInternalChange(bin,p[0],p[1],0,true);
    const double before=p[0].time_merger;
    bin.ecc+=1e-5;
    ar.checkBHMergerDynamicalChange(bin,p[0],p[1],0);
    assert(countRows(path)==1 && p[0].time_merger<before);
    // Dynamics may lengthen the clock: never take min(old,new).
    bin.semi*=2;
    const double shortened=p[0].time_merger;
    ar.checkBHMergerDynamicalChange(bin,p[0],p[1],0);
    assert(countRows(path)==2 && p[0].time_merger>shortened);
    // Six small -2% lifetime increments survive interleaved internal updates.
    bin.ecc=.5;bin.semi=1;
    ar.syncBHMergerInternalChange(bin,p[0],p[1],0,true);
    for(int i=0;i<6;++i) {
        bin.semi*=std::exp(-.02/4);
        ar.checkBHMergerDynamicalChange(bin,p[0],p[1],i+1);
        if(i<5) assert(countRows(path)==2);
        const double pending=p[0].merger_dyn_log_tgw;
        bin.semi*=.99;bin.ecc-=.001;
        ar.syncBHMergerInternalChange(bin,p[0],p[1],i+1);
        assert(p[0].merger_dyn_log_tgw==pending);
        PtclHard copy;copy.DataCopy(p[0]);assert(copy.merger_dyn_log_tgw==pending);
    }
    assert(countRows(path)==3);
    // Significant dynamics immediately before BSE remains recorded exactly once.
    bin.semi*=.5;ar.checkBHMergerDynamicalChange(bin,p[0],p[1],7);
    const double dynamic_deadline=p[0].time_merger;
    bin.semi*=.8;ar.syncBHMergerInternalChange(bin,p[0],p[1],7);
    ar.checkBHMergerDynamicalChange(bin,p[0],p[1],7);
    assert(countRows(path)==4 && p[0].time_merger==dynamic_deadline);
    // Internal disruption/capture must not be logged as dynamics.
    bin.semi=-1;bin.ecc=1.1;ar.syncBHMergerInternalChange(bin,p[0],p[1],8);
    assert(p[0].time_merger<0 && p[0].merger_dyn_log_tgw==0);
    bin.semi=1;bin.ecc=.9;ar.syncBHMergerInternalChange(bin,p[0],p[1],9,true);
    assert(p[0].time_merger>=9 && countRows(path)==4);
    // A companion exchange starts a new history without inventing a change.
    p[1].id=3;ar.checkBHMergerDynamicalChange(bin,p[0],p[1],10);
    assert(countRows(path)==4 && p[0].getBinaryPairID()==3);
    ar.change_output.close();std::remove(path.c_str());
    std::cout<<"clock separation: 2000 internal updates, high-e dynamics, delays, cumulative weak perturbations, exchange/capture passed\n";
}

// Exercise the actual Fortran BSE + postProcess + AR interrupt dispatch.
// No orbit is prescribed between callbacks: changes here come only from BSE.
void realBSEClockTest(double eccentricity, int callbacks, bool finish_merger=false) {
    const std::string path="/tmp/petar-real-bse-"+std::to_string(getpid());
    ARInteraction ar; IOParamsBSE io;
    io.rscale.value=4.435e7;io.vscale.value=.9777922217;
    ar.bse_manager.initial(io,false);
    ar.gravitational_constant=.00449830997959438;ar.eps_sq=0;
    ar.stellar_evolution_option=1;ar.stellar_evolution_write_flag=false;
    ar.tide.gravitational_constant=ar.gravitational_constant;
    ar.tide.speed_of_light=ar.bse_manager.getSpeedOfLight();
    ar.change_output.openChange(path,false,true,17);
    PtclHard p[2];
    for(int i=0;i<2;++i) {
        p[i].id=i+1;p[i].mass=30;p[i].star.initial(30,14);
        StarParameterOut out;ar.bse_manager.evolveStar(p[i].star,out,0,true);
        p[i].radius=ar.bse_manager.getMergerRadius(p[i].star);
        p[i].time_record=0;p[i].time_interrupt=0;p[i].dm=0;p[i].spin=0;
        p[i].setBinaryPairID(2-i);p[i].setBinaryInterruptState(BinaryInterruptState::none);
    }
    AR::BinaryTree<PtclHard> bin;bin.setMembers(p,p+1,0,1);
    bin.m1=bin.m2=30;bin.mass=60;bin.ecc=eccentricity;
    // Keep pericentre detached even in the very eccentric case.
    bin.semi=1e-8/(1-eccentricity);
    bin.incline=bin.rot_horizon=bin.rot_self=0;bin.ecca=2;
    bin.slowdown.setSlowDownFactor(1);
    bin.calcParticles(ar.gravitational_constant);bin.calcCenterOfMass();bin.calcOrbit(ar.gravitational_constant);
    ar.checkBHMergerDynamicalChange(bin,p[0],p[1],0);
    const double deadline=p[0].time_merger, initial_a=bin.semi;
    assert(deadline>0);
    int evolved=0;
    for(int i=1;i<=callbacks;++i) {
        p[0].time_interrupt=p[1].time_interrupt=0;
        AR::InterruptBinary<PtclHard> intr;
        intr.time_now=deadline*.2*i/callbacks;intr.time_end=intr.time_now;
        ar.modifyAndInterruptIter(intr,bin);
        if(p[0].star.tphys>0) ++evolved;
        assert(p[0].mass>0 && p[1].mass>0);
        assert(p[0].time_merger==deadline && p[1].time_merger==deadline);
        assert(countRows(path)==0);
    }
    assert(evolved==callbacks && bin.semi<initial_a);
    if (finish_merger) {
        // Keep calling real BSE all the way to the fixed deadline. In particular,
        // repeated deferred contact/coalescence callbacks must not restart it.
        for(int i=1;i<=80;++i) {
            p[0].time_interrupt=p[1].time_interrupt=0;
            AR::InterruptBinary<PtclHard> intr;
            intr.time_now=deadline*(.2+.8*i/80.);intr.time_end=intr.time_now;
            ar.modifyAndInterruptIter(intr,bin);
            if(i<80) {
                assert(p[0].mass>0 && p[1].mass>0);
                assert(p[0].time_merger==deadline);
            }
            else {
                assert(intr.status==AR::InterruptStatus::merge);
                assert(p[0].mass==0 || p[1].mass==0);
            }
            assert(countRows(path)==0);
        }
        ar.change_output.close();std::remove(path.c_str());
        std::cout<<"actual BSE merger at fixed deadline: e="<<eccentricity
                 <<" initial_callbacks="<<callbacks<<" t="<<deadline<<"\n";
        return;
    }
    // A genuine dynamical velocity impulse followed by BSE must be detected
    // BEFORE that internal update, with no duplicated post-BSE change.
    auto cm=(p[0].vel+p[1].vel)*.5;
    p[0].vel=cm+(p[0].vel-cm)*.8;
    p[1].vel=cm+(p[1].vel-cm)*.8;
    bin.calcOrbit(ar.gravitational_constant);
    assert(ar.isBHMergerOrbitClearlyPerturbed(bin,p[0],p[1]));
    p[0].time_interrupt=p[1].time_interrupt=0;
    AR::InterruptBinary<PtclHard> intr;intr.time_now=deadline*.20001;intr.time_end=intr.time_now;
    ar.modifyAndInterruptIter(intr,bin);
    assert(countRows(path)==1);
    ar.change_output.close();std::remove(path.c_str());
    std::cout<<"actual BSE dispatch e="<<eccentricity<<" callbacks="<<callbacks
             <<" deadline="<<deadline<<" a_after/a_initial="<<bin.semi/initial_a
             <<": zero internal change rows; one external impulse row\n";
}

int main() {
    legacyLogTest();clockSeparationTest();
    for(double e: {0.0,0.9,0.9999})
        for(int n: {10,100}) realBSEClockTest(e,n);
    for(double e: {0.0,0.9,0.9999})
        for(int n: {10,100}) realBSEClockTest(e,n,true);
}
