#include "petar.hpp"
#include <vector>
#include <iostream>

// Independent invariants: component cutoff masses and the analytical
// orbit-averaged second moment, including unequal masses and rotated orbits.
void checkMoments(double m1,double m2,double ecc,double angle) {
    COMM::BinaryTree<FPSoft,COMM::Binary> bin;
    bin.m1=m1; bin.m2=m2; bin.mass=m1+m2; bin.semi=0.2; bin.ecc=ecc;
    bin.incline=angle; bin.rot_horizon=0.3; bin.rot_self=-0.5;
    bin.pos=PS::F64vec(1.,-2.,3.); bin.vel=PS::F64vec(2.,3.,4.);
    PseudoParticleMultipoleManager manager;
    std::vector<FPSoft> sources(manager.getParticleN());
    manager.createSampleParticles(sources.data(),bin);
    double mass[2]={0,0}, moment[3][3]={}; PS::F64vec dipole(0.);
    for(int i=0;i<int(sources.size());++i) {
        auto& p=sources[i]; mass[i%2]+=p.mass;
        const auto dr=p.pos-bin.pos; dipole+=p.mass*dr;
        for(int j=0;j<3;++j)for(int k=0;k<3;++k)moment[j][k]+=p.mass*dr[j]*dr[k];
        assert((p.vel-bin.vel)*(p.vel-bin.vel)==0.);
    }
    for(int i=0;i<2;++i)assert(std::abs(mass[i]/(i==0?m1:m2)-1.)<1e-12);
    assert(std::sqrt(dipole*dipole)<1e-12*bin.mass);
    PS::F64vec x(1.,0.,0.),y(0.,1.,0.);
    bin.rotateToOriginalFrame(&x.x); bin.rotateToOriginalFrame(&y.x);
    double mu=m1*m2/(m1+m2), scale=mu*bin.semi*bin.semi;
    for(int j=0;j<3;++j)for(int k=0;k<3;++k) {
        double ref=scale*((0.5+2*ecc*ecc)*x[j]*x[k]+0.5*(1-ecc*ecc)*y[j]*y[k]);
        assert(std::abs(moment[j][k]-ref)<1e-12*scale);
    }
}

struct Particles: std::vector<FPSoft> {
    using std::vector<FPSoft>::vector;
    int getNumberOfParticleLocal() const {return size();}
};
struct Neighbors {
    std::vector<EPJSoft> sources;
    int getNeighborListOneParticle(const FPSoft&,EPJSoft*& out) {out=sources.data(); return sources.size();}
};

void checkCoincidentSources() {
    ForceSoft::grav_const=1.; EPISoft::eps=0.; EPISoft::r_out=2.;
    Particles p(2); Neighbors tree; tree.sources.resize(2);
    for(int i=0;i<2;++i) {
        p[i].id=i+1; p[i].mass=i?1028.7321917245959/3:10.;
        p[i].pos=0.; p[i].acc=PS::F64vec(i,0.,0.);
        p[i].changeover.setR(i?2.:1.,i?4.:2.);
        p[i].changeover.r_scale_next=i?1.:2.;
        p[i].group_data.artificial.setParticleTypeToArtificial(i+1);
        p[i].r_search=4.;
    }
    for(int i=0;i<2;++i) {
        p[i].pot_soft=p[i].pot_tot=-p[1-i].mass/EPISoft::r_out;
        tree.sources[i].copyFromFP(p[i]);
#ifdef KDKDK_4TH
        p[i].acorr=-2.*p[1-i].mass/8.*(p[i].acc-p[1-i].acc);
#endif
    }
    ArtificialParticleManager ap;
    // Exercise the actual neighbor correction, including the EPJSoft overload.
    SystemHard::correctForceWithCutoffTreeNeighborOMP<Particles,FPSoft,Neighbors,EPJSoft>(p,tree,2,ap);
    for(int i=0;i<2;++i) {
        assert(std::isfinite(p[i].acc.x)); assert(p[i].acc.x==i);
        assert(std::abs(p[i].pot_tot)<1e-12 && std::abs(p[i].pot_soft)<1e-12);
    }
#ifdef KDKDK_4TH
    SystemHard::correctForceWithCutoffTreeNeighborOMP<Particles,FPSoft,Neighbors,EPJSoft>(p,tree,2,ap,true);
    for(const auto& q:p)assert(std::sqrt(q.acorr*q.acorr)<1e-12);
#endif
}

int main() {
    for(double ratio:{0.001,1.,1028.7321917245959/30.,1000.})
        for(double ecc:{0.,0.7,0.99})for(double angle:{0.,0.8})checkMoments(30.,30.*ratio,ecc,angle);
    checkCoincidentSources();
    std::cout<<"quadrupole cutoff masses, moments and coincident-source corrections passed\n";
}
