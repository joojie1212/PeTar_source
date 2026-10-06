#include <particle_simulator.hpp>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#if defined(USE_SIMD) && defined(INTRINSIC_X86)
#ifdef CALC_EP_64bit
#define RSQRT_NR_EPJ_X4
#else
#define RSQRT_NR_EPJ_X2
#endif
#include "phantomquad_for_p3t_x86.hpp"
#endif
#include "darkmatter_force.hpp"

int main(int argc, char** argv) {
    PS::Initialize(argc, argv);
    const int rank=PS::Comm::getRank(), nrank=PS::Comm::getNumberOfProc();
    PS::ParticleSystem<DMTreeParticle> particles;
    particles.initialize();
    // All targets start on rank 0: exchange must distribute them correctly.
    particles.setNumberOfParticleLocal(rank==0 ? 80 : 0);
    if (rank==0) for (int i=0; i<80; ++i) {
        particles[i].pos=PS::F64vec(0.13*i, std::sin(i), std::cos(i));
        particles[i].mass=i<8 ? 1.0+0.1*i : 0.0;
        particles[i].active=true;
    }
    PS::DomainInfo domain;
    domain.initialize(); domain.decomposeDomainAll(particles);
    particles.exchangeParticle(domain);
    DMTree tree;
    tree.initialize(80, 0.0, 2, 8);
    const double eps=0.2, G=1.7;
    DMGravityKernel kernel(eps,G);
    tree.calcForceAll(kernel,kernel,particles,domain);
    double error=0;
    for (int i=0; i<particles.getNumberOfParticleLocal(); ++i) {
        PS::F64vec acc(0); double pot=0;
        for (int j=0; j<8; ++j) {
            const auto dr=particles[i].pos-PS::F64vec(0.13*j,std::sin(j),std::cos(j));
            const double ri=1/std::sqrt(dr*dr+eps*eps), mass=1+0.1*j;
            acc-=G*mass*ri*ri*ri*dr; pot-=G*mass*ri;
        }
        const auto f=tree.getForce(i);
        if (!std::isfinite(f.pot) || !std::isfinite(f.acc*f.acc)) PS::Abort(2);
        const auto da=f.acc-acc;
        error=std::max(error,std::sqrt(da*da)/(1+std::sqrt(acc*acc)));
        error=std::max(error,std::abs(f.pot-pot)/(1+std::abs(pot)));
    }
    // Independent cross-energy/self-removal test: one star (m=2) and one DM
    // particle (m=3), separated by one length unit. No 1/2 in cross energy.
    DMEPI target; target.pos=PS::F64vec(1,0,0); target.active=true;
    DMEPJ star; star.pos=PS::F64vec(0); star.mass=2;
    DMForce f; f.clear(); kernel(&target,1,&star,1,&f);
    const double cross=3*f.pot, expected=-G*6/std::sqrt(1+eps*eps);
    error=std::max(error,std::abs(cross/expected-1));
    DMEPJ self; self.pos=target.pos; self.mass=3;
    f.clear(); kernel(&target,1,&self,1,&f);
    error=std::max(error,std::abs(f.pot+G*3/eps)/(G*3/eps));
    // Exercise EP-SP, including a zero-mass moment, independently of theta.
    PS::SPJMonopole sp; sp.mass=star.mass; sp.pos=star.pos;
    f.clear(); kernel(&target,1,&sp,1,&f);
    error=std::max(error,std::abs(3*f.pot/expected-1));
    DMMoment empty; empty.set(); sp.copyFromMoment(empty);
    kernel(&target,1,&sp,1,&f);
    if (!std::isfinite(f.pot)) PS::Abort(3);
    error=PS::Comm::getMaxValue(error);
#if defined(USE_SIMD) && !defined(CALC_EP_64bit)
    const double tolerance=3e-5;
#else
    const double tolerance=2e-11;
#endif
    if (rank==0) std::cout << "ranks=" << nrank << " max_relative_error=" << error << std::endl;
    if (error>tolerance) PS::Abort(4);
    PS::Finalize();
}
