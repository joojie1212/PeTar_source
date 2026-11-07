#include "petar.hpp"
#include <cassert>
#include <fstream>
#include <sstream>

using P = EscapeCandidate;
P particle(long long id,double m,PS::F64vec x,PS::F64vec v,long long mate=0) {
    P p;p.id=id;p.mass=m;p.pos=x;p.vel=v;p.partner=mate;return p;
}
double directEnergy(const std::vector<P>& p,double G,double eps2) {
    double e=0;
    for(size_t i=0;i<p.size();++i) {
        e+=.5*p[i].mass*(p[i].vel*p[i].vel);
        for(size_t j=0;j<i;++j) if(p[i].mass>0 && p[j].mass>0)
            e+=p[i].mass*p[j].mass*EscapeRemoval::pairPotential(p[i],p[j],G,eps2);
    }
    return e;
}
void run(const std::string& name,const std::vector<P>& all,std::vector<long long> expected,
         double radius=5,double G=1,double eps2=0,bool check_energy=true,
         std::vector<long long> forced={},PS::F64vec cmvel=0.0) {
    PeTar sim;
    const int rank=PS::Comm::getRank(),nr=PS::Comm::getNumberOfProc();
    // FDPS 7's tree destructor assumes initialize() has allocated its buffers.
    sim.tree_nb.initialize(all.size(),.5,8,32);
    sim.tree_soft.initialize(all.size(),.5,8,32);
    sim.system_soft.initialize();
    std::vector<P> local;
    for(size_t i=0;i<all.size();++i) if(int(i)%nr==rank) local.push_back(all[i]);
    sim.system_soft.setNumberOfParticleLocal(local.size());
    sim.stat.n_real_loc=local.size();sim.stat.n_real_glb=all.size();
    sim.stat.pcm.pos=0.;sim.stat.pcm.vel=cmvel;
    sim.stat.time=1.;sim.input_parameters.gravitational_constant.value=G;
    sim.input_parameters.eps.value=std::sqrt(eps2);sim.input_parameters.write_style.value=0;
    sim.escaper.r_escape_sq=radius*radius;sim.escaper.check_energy_flag=check_energy;
    for(size_t i=0;i<local.size();++i) {
        auto& p=sim.system_soft[i];const auto& a=local[i];
        p.id=a.id;p.mass=a.mass;p.pos=a.pos;p.vel=a.vel;
        p.setBinaryPairID(a.partner);
        p.group_data.artificial.setParticleTypeToSingle();
        if(p.mass==0) p.group_data.artificial.setParticleTypeToUnused();
        // Deliberately wrong old tree potential. It must not affect either
        // the physical decision or the deletion energy.
        p.pot_tot=-100.;p.pot_soft=1e30;
        if(std::find(forced.begin(),forced.end(),p.id)!=forced.end()) sim.remove_list.push_back(i);
    }
    const double before=directEnergy(all,G,eps2);
    sim.stat.energy.etot_ref=before;sim.stat.energy.etot_sd_ref=before;
    sim.removeParticles();
    std::vector<P> remainlocal;
    for(int i=0;i<sim.system_soft.getNumberOfParticleLocal();++i) {
        const auto& p=sim.system_soft[i];remainlocal.push_back(particle(p.id,p.mass,p.pos,p.vel,p.getBinaryPairID()));
    }
    auto remain=EscapeRemoval::gather(remainlocal);
    std::vector<long long> removed;
    for(const auto& p:all)
        if(std::none_of(remain.begin(),remain.end(),[&](const P& r){return r.id==p.id;})) removed.push_back(p.id);
    std::sort(removed.begin(),removed.end());std::sort(expected.begin(),expected.end());
    assert(removed==expected);
    const double after=directEnergy(remain,G,eps2),delta=sim.stat.energy.de_change_cum;
    const double scale=std::max(1.,std::abs(before)+std::abs(after));
    assert(std::abs(after-before-delta)<2e-12*scale);
    assert(std::abs(after-sim.stat.energy.etot_ref)<2e-12*scale);
    assert(PS::Comm::getMaxValue(delta)==PS::Comm::getMinValue(delta));
    assert(sim.stat.n_remove_glb==int(expected.size()));
    if(forced.empty()) {
        int npositive=0;for(const auto& p:all) if(p.mass>0 && std::find(expected.begin(),expected.end(),p.id)!=expected.end()) ++npositive;
        assert(sim.stat.n_escape_glb==npositive);
    }
    if(rank==0) std::cout<<std::setprecision(17)<<"PASS "<<name<<" ranks="<<nr
        <<" removed="<<expected.size()<<" delta="<<delta
        <<" residual="<<(after-before-delta)/scale<<std::endl;
}

int main(int argc,char**argv) {
    PeTar::initialFDPS(argc,argv);
    if(argc==2) {
        std::ifstream f(argv[1]);int n;double G,radius;f>>n>>G>>radius;
        std::vector<P> p(n);std::vector<long long> expected;
        for(auto& a:p) f>>a.id>>a.partner>>a.mass>>a.pos.x>>a.pos.y>>a.pos.z>>a.vel.x>>a.vel.y>>a.vel.z;
        long long id;while(f>>id) expected.push_back(id);
        run("historical_escape_state",p,expected,radius,G);
    }
    else {
        auto core=particle(1,100,PS::F64vec(0,0,0),PS::F64vec(0));
        const double v=std::sqrt(1./.004);
        auto a=particle(2,1,PS::F64vec(10,-.001,0),PS::F64vec(10+v,0,0),3);
        auto b=particle(3,1,PS::F64vec(10,.001,0),PS::F64vec(10-v,0,0),2);
        // Old per-member rule using cached phi=-100 removes only a.
        assert(.5*(a.vel*a.vel)-100>0 && .5*(b.vel*b.vel)-100<0);
        run("binary_atomic_old_half_delete",{core,a,b},{2,3});
        auto aa=a,bb=b;aa.vel.x-=10;bb.vel.x-=10;
        run("bound_cm_stays",{core,aa,bb},{});
        run("radius_only_atomic",{core,aa,bb},{2,3},5,1,0,false);
        auto boosted=std::vector<P>{core,a,b};
        for(auto& p:boosted) p.vel+=PS::F64vec(100,20,-30);
        run("moving_cluster_frame",boosted,{2,3},5,1,0,true,{},PS::F64vec(100,20,-30));
        run("single_escape",{core,particle(2,1,PS::F64vec(10,0,0),PS::F64vec(10,0,0))},{2});
        aa=a;bb=b;aa.pos.y=bb.pos.y=0;aa.pos.x=4.99;bb.pos.x=5.03;
        aa.vel=PS::F64vec(10,3,0);bb.vel=PS::F64vec(10,-3,0);
        run("pair_cm_out_member_inside",{core,aa,bb},{2,3});
        aa.pos.x=4.97;bb.pos.x=5.01;
        run("pair_cm_inside_member_out",{core,aa,bb},{});
        aa=a;bb=b;aa.vel=PS::F64vec(100,0,0);bb.vel=PS::F64vec(0);
        run("unbound_stale_pair_ids",{core,aa,bb},{2});
        aa.partner=999;
        run("missing_partner",{core,aa},{2});
        aa=a;bb=b;aa.pos.x=bb.pos.x=1;
        run("explicit_remove_promotes_pair",{core,aa,bb},{2,3},5,1,0,true,{2});
        run("softened_pair_energy",{core,a,b},{2,3},5,1,1e-6);
        aa=a;bb=b;aa.id=4;bb.id=5;aa.partner=5;bb.partner=4;
        aa.pos.x=bb.pos.x=20;
        run("two_pairs_cross_energy",{core,a,b,aa,bb},{2,3,4,5});
        auto z=particle(9,0,PS::F64vec(2,0,0),PS::F64vec(0));
        run("unused_zero_mass",{core,z},{9});
    }
    PeTar::finalizeFDPS();
}
