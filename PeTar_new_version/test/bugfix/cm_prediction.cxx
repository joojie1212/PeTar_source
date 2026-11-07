#define SMOOTH_CM_USING_VEL_PRED
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <set>
#include <sstream>
#define ASSERT assert
#define DATADUMP(x) ((void)0)
#include "Common/particle_group.h"
#include "Common/binary_tree.h"
#include "Common/list.h"
#include "AR/information.h"
#include "AR/symplectic_step.h"
#include "Hermite/block_time_step.h"
#include "particle_simulator.hpp"
#include "matrix2.hpp"
#include "changeover.hpp"
#include "ptcl.hpp"
#include "energy.hpp"
#include "../../src/profile.hpp"
#include "particle_distribution_generator.hpp"

#include "status.hpp"
struct CMParticle {
    double mass=1.,pot_soft=-1.;
    PS::F64vec pos=0.,vel=0.,acc=0.;
};
int main() {
    Status status;
    status.time=1.; status.pcm.time_record=0.;
    status.pcm.is_center_shift_flag=true;
    status.pcm.vel=PS::F64vec(0.,0.,1.);
    CMParticle p; p.pos=PS::F64vec(0.,0.,5.); p.vel=status.pcm.vel;
    status.calcAndShiftCenterOfMass(&p,1,1);
    assert(status.pcm.pos.z==5.);
    std::cout<<"CM z-axis prediction regression passed\n";
}
