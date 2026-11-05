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
#include "changeover.hpp"
#include "ptcl.hpp"
#include "energy.hpp"

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
int main() {
    testTrees(); testPersistence(); testListAndStep();
    ChangeOver c; c.setR(1.,2.);
    assert(std::isfinite(c.calcPotW(1.5)));
    assert(std::abs(c.calcPotW(2.0))<1e-12);
    FILE* file=tmpfile(); c.writeBinary(file); rewind(file); ChangeOver d; d.readBinary(file); fclose(file);
    assert(c.calcPotW(1.5)==d.calcPotW(1.5));
    std::cout << "core regressions passed\n";
}
