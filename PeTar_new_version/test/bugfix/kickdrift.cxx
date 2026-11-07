#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include "particle_simulator.hpp"
#include "kickdriftstep.hpp"
int main(int argc, char**) {
    if(argc>1) { KickDriftStep unset; unset.nextOneStep(); return 1; }
    for(bool dkd: {false,true}) {
        KickDriftStep step;
        if(dkd) step.setDKDMode();
        step.setStep(2.);
        for(int cycle=0;cycle<2;++cycle) {
            double kick=0.,drift=0.;
            (dkd?drift:kick)+=step.getDtStartContinue();
            // Complete all interior pairs; the last second operation is the end.
            while(true) {
                (dkd?kick:drift)+=dkd?step.getDtKickContinue():step.getDtDriftContinue();
                step.nextContinue();
                if(step.isNextEndPossible()) break;
                (dkd?drift:kick)+=dkd?step.getDtDriftContinue():step.getDtKickContinue();
            }
            (dkd?drift:kick)+=step.getDtEndContinue();
            assert(std::abs(kick-2.)<1e-14 && std::abs(drift-2.)<1e-14);
        }
    }
    std::cout<<"KDK/DKD step regressions passed\n";
}
