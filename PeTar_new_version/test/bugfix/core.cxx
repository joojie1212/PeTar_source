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

struct Particle {
    Float mass=1.0;
    PS::F64vec pos, vel;
    int id=0;
};
using Tree = COMM::BinaryTree<Particle, COMM::Binary>;

void testTrees() {
    Particle p[6];
    for(int i=0;i<6;++i) { p[i].id=i; p[i].pos=PS::F64vec(i,0,0); p[i].vel=PS::F64vec(0,0.1*i,0); }
    Tree source[5], target[5];
    source[0].setMembers(&p[0],&p[1],0,1);
    source[1].setMembers(&source[0],&p[2],-1,2);
    source[2].setMembers(&p[3],&p[4],3,4);
    source[3].setMembers(&source[2],&p[5],-1,5);
    source[4].setMembers(&source[1],&source[3],-1,-1);
    source[4].getherBinaryTreeIter(target);
    std::set<int> leaves;
    std::set<Tree*> nodes;
    std::function<void(Tree*)> visit = [&](Tree* tree) {
        assert(tree>=target && tree<target+5);
        assert(nodes.insert(tree).second);
        for (int k=0;k<2;++k) {
            if(tree->isMemberTree(k)) visit(tree->getMemberAsTree(k));
            else assert(leaves.insert(tree->getMember(k)->id).second);
        }
    };
    visit(&target[4]);
    assert(leaves.size()==6 && nodes.size()==5);
    for(int n=3;n<=6;++n) {
        int index[6]={0,1,2,3,4,5};
        Tree::generateBinaryTree(target,index,n,p,1.0);
        assert(target[n-2].getMemberN()==n);
    }
}
void testMasslessTreeRebuild() {
    COMM::ParticleGroup<Particle,Particle> group;
    group.setMode(COMM::ListMode::local); group.reserveMem(4);
    for(int i=0;i<4;++i) {
        Particle p; p.id=i; p.pos=PS::F64vec(i,0.,0.);
        p.vel=PS::F64vec(0.,0.1*i,0.); group.addMember(p);
    }
    AR::Information<Particle,Particle> info; info.reserveMem(4);
    info.generateBinaryTree(group,1.);
    for(int live: {2,1,0}) {
        for(int i=0;i<4;++i) group[i].mass=i<live ? 1. : 0.;
        for(int i=0;i<3;++i) {
            info.binarytree[i].semi=13.; info.binarytree[i].r=13.;
            info.binarytree[i].ecc=0.75; info.binarytree[i].stab=0.;
        }
        info.generateBinaryTree(group,1.);
        auto& root=info.getBinaryTreeRoot();
        assert(root.mass==live && root.semi==0. && root.ecc==0.);
        assert(root.period==0. && root.stab==NUMERIC_FLOAT_MAX);
        assert(std::isfinite(root.r) && root.r!=13.);
        for(int i=0;i<3;++i) assert(std::isfinite(root.pos[i]));
        if(live==2) assert(info.binarytree[0].period>0.);
    }
    for(int i=0;i<4;++i) group[i].pos=0.;
    info.generateBinaryTree(group,1.);
    assert(info.getBinaryTreeRoot().r==0.);
}
void testPersistence() {
    AR::Information<Particle, Particle> original, recovered;
    original.ds=0.001234567890123; original.time_offset=13.5;
    FILE* file=tmpfile(); assert(file);
    original.writeBinary(file); rewind(file); recovered.readBinary(file); fclose(file);
    assert(recovered.ds==original.ds && recovered.time_offset==original.time_offset);
    AR::SymplecticStep step, restored;
    for(int order: {2,4,-6,-8,6}) {
        step.initialSymplecticCofficients(order);
        file=tmpfile(); step.writeBinary(file);
        for(int i=0;i<50;++i) { rewind(file); restored.readBinary(file); assert(restored.getOrder()==std::abs(order)); }
        fclose(file);
    }
    step.clear(); restored.clear();
    EnergyAndMomentum e; e.error_cum_pre=3.25; e.ekin=9.0; e.L_ref.x=15.0;
    file=tmpfile(); e.writeBinary(file); assert(ftell(file)==sizeof(e));
    rewind(file); EnergyAndMomentum read; assert(fread(&read,sizeof(read),1,file)==1); fclose(file);
    assert(read.error_cum_pre==3.25 && read.ekin==9.0 && read.L_ref.x==15.0);
}
void testListAndStep() {
    for(auto mode: {COMM::ListMode::local, COMM::ListMode::copy}) {
        COMM::List<int> a,b; a.setMode(mode); a.reserveMem(3);
        int value=42;
        if(mode==COMM::ListMode::copy) a.addMemberAndAddress(value); else a.addMember(value);
        b=a; a=a; assert(a[0]==42 && b[0]==42); a[0]=7; assert(b[0]==42);
    }
    H4::BlockTimeStep4th step; step.setDtRange(1.0,10);
    const double small=step.getDtMin();
    assert(step.calcNextDtLimit(0)==1.0);
    assert(step.calcNextDtLimit(std::nextafter(small,0.0))==small);
    assert(step.calcNextDtLimit(small)==small);
    assert(step.calcNextDtLimit(2*small)==2*small);
    assert(step.calcNextDtLimit(1.0)==1.0);
}
void testNumericalAndArrayFixes() {
    for(double e: {0.,0.3,0.5,0.9,0.99})
        for(int i=0;i<=100;++i) {
            const double mean=2.*std::acos(-1.)*i/100.;
            const double u=COMM::Binary::calcEccAnomaly(mean,e);
            assert(std::abs(u-e*std::sin(u)-mean)<1e-10);
        }
    for(double half: {1e-300,1.,1e300}) {
        PS::MortonKey key; key.initialize(half,0.);
        const auto low=key.getKey(PS::F64vec(-half,-half,-half));
        const auto high=key.getKey(PS::F64vec(half,half,half));
        assert(low<high);
    }
    PS::Matrix2<double> a(2.,3.,5.,7.), b(11.,13.,17.,19.);
    auto c=a*b;
    assert(c.xx==117. && c.xy==99. && c.yx==134. && c.yy==158.);
    SysProfile timing; timing.total.time=10.; timing.total.tbar=3.;
    assert(timing.getMax().total.time==10. && timing.getMax().total.tbar==3.);
    assert(timing.getMin().total.time==10. && timing.getMin().total.tbar==3.);
    SysCounts counts;
    for (auto field: SysCounts::counters) counts.*field=7;
    counts.clear();
    for (auto field: SysCounts::counters) assert((counts.*field).n==0);
    PS::TimeProfile profile;
    assert(profile.calc_force__core__walk_tree==0.);
    assert(profile.calc_force__core__keep_list==0.);
    assert(profile.calc_force__core__copy_ep==0.);
    assert(profile.calc_force__core__dispatch==0.);
    assert(profile.calc_force__core__retrieve==0.);
    PS::ReallocatableArray<int> values;
    values.resizeNoInitialize(8);
    for(int i=0;i<8;++i) values[i]=i*i;
    values.reserve(4096);
    assert(values.size()==8);
    for(int i=0;i<8;++i) assert(values[i]==i*i);
}
void testGeneratorScalingAndSeed() {
    double *m0,*m1; PS::F64vec *p0,*p1,*v0,*v1;
    ParticleDistributionGenerator::makePlummerModel(1.,8,8,m0,p0,v0,-.25,0,1.);
    for (auto factors: {std::array<double,2>{1.,3.},std::array<double,2>{2.,1.}}) {
        const double mass=factors[0],g=factors[1];
        ParticleDistributionGenerator::makePlummerModel(mass,8,8,m1,p1,v1,-.25*g*mass*mass,0,g);
        for(int i=0;i<8;++i) {
            assert(std::abs(m1[i]-mass*m0[i])<1e-14);
            assert((p1[i]-p0[i])*(p1[i]-p0[i])<1e-24);
            const auto delta=v1[i]-std::sqrt(g*mass)*v0[i];
            assert(delta*delta<1e-24);
        }
        delete[] m1;delete[] p1;delete[] v1;
    }
    delete[] m0;delete[] p0;delete[] v0;
    double total;
    ParticleDistributionGenerator::makeKeplerDisk(total,m0,p0,v0,8,8,1.,2.,.01,.01,10.,1.,0.,1.,-1.5,7);
    ParticleDistributionGenerator::makeKeplerDisk(total,m1,p1,v1,8,8,1.,2.,.01,.01,10.,1.,0.,1.,-1.5,7);
    for(int i=0;i<8;++i) {
        assert((p1[i]-p0[i])*(p1[i]-p0[i])==0.);
        assert((v1[i]-v0[i])*(v1[i]-v0[i])==0.);
    }
    delete[] m0;delete[] p0;delete[] v0;
    delete[] m1;delete[] p1;delete[] v1;
}

int main() {
    testGeneratorScalingAndSeed(); testTrees(); testMasslessTreeRebuild(); testPersistence(); testListAndStep(); testNumericalAndArrayFixes();
    ChangeOver c; c.setR(1.,2.);
    assert(std::isfinite(c.calcPotW(1.5)));
    assert(std::abs(c.calcPotW(2.0))<1e-12);
    FILE* file=tmpfile(); c.writeBinary(file); rewind(file); ChangeOver d; d.readBinary(file); fclose(file);
    assert(c.calcPotW(1.5)==d.calcPotW(1.5));
    std::cout << "core regressions passed\n";
}
