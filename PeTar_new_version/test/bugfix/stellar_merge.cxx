#include "bse_interface.h"
#include <algorithm>

extern "C" void hrdiag_(double*,double*,double*,double*,double*,double*,
    double*,double*,double*,double*,double*,int*,double*,double*,double*,
    double*,double*,double*,double*,double*,int*);

static double inertia(StarParameter s, BSEManager& manager) {
    double tm,tn,ts[20],ls[10],gb[10],menv,renv,k2,frac,total,mco;
    int ecs,kw=static_cast<int>(s.kw);
    double age=s.tphys-s.epoch;
    star_(&kw,&s.m0,&s.mt,&tm,&tn,ts,ls,gb,manager.zpars);
    hrdiag_(&s.m0,&age,&s.mt,&tm,&tn,ts,ls,gb,manager.zpars,
            &s.r,&s.lum,&kw,&s.mc,&s.rc,&menv,&renv,&k2,&frac,&total,&mco,&ecs);
    return k2*s.r*s.r*(s.mt-s.mc)+0.21*s.rc*s.rc*s.mc;
}

static void close(double x, double y) {
    assert(std::isfinite(x) && std::isfinite(y));
    assert(std::abs(x-y)<1e-8*std::max({1.,std::abs(x),std::abs(y)}));
}
int main() {
    IOParamsBSE params;
    BSEManager manager; manager.initial(params,false);
    // Real MIX path: two main-sequence stars, followed by an independent
    // zero-duration evolution to check the returned stellar structure.
    StarParameter a,b; a.initial(1.); b.initial(0.8);
    StarParameterOut oa{},ob{};
    manager.evolveStar(a,oa,0.); manager.evolveStar(b,ob,0.);
    a.ospin=17.; b.ospin=3.;
    const double initial_spin_momentum=a.ospin*inertia(a,manager);
    const double old_radius=a.r;
    double semi=-10., ecc=1.5;
    manager.merge(a,b,oa,ob,semi,ecc);
    assert(a.mt>1. && b.mt==0. && b.kw==15);
    close(a.ospin*inertia(a,manager),initial_spin_momentum);
    assert(a.r!=old_radius && a.lum>0. && oa.tm>0.);
    assert(b.r==0. && b.mc==0. && b.rc==0. && b.lum==0.);
    assert(b.ospin==0. && ob.menv==0. && ob.renv==0.);
    StarParameter reference=a; StarParameterOut out{};
    manager.evolveStar(reference,out,0.);
    close(a.r,reference.r); close(a.lum,reference.lum);
    close(a.mc,reference.mc); close(a.rc,reference.rc);
    close(oa.menv,out.menv); close(oa.renv,out.renv);
    std::cout<<"MIX remnant mass/radius/luminosity: "
             <<a.mt<<" "<<a.r<<" "<<a.lum<<"\n";
    // Real CE path: a core-helium-burning donor overflows at pericentre.
    a.initial(5.); b.initial(1.);
    oa=StarParameterOut{}; ob=StarParameterOut{};
    manager.evolveStar(a,oa,100.); manager.evolveStar(b,ob,100.);
    assert(a.kw==4);
    semi=-10.; ecc=1.5;
    manager.merge(a,b,oa,ob,semi,ecc);
    assert(a.mt>0. && b.mt==0.);
    reference=a; out=StarParameterOut{};
    manager.evolveStar(reference,out,0.);
    close(a.r,reference.r); close(a.lum,reference.lum);
    close(a.mc,reference.mc); close(a.rc,reference.rc);
    close(oa.menv,out.menv); close(oa.renv,out.renv);
    assert(std::isfinite(a.ospin) && b.ospin==0.);
    std::cout<<"CE remnant mass/radius/luminosity: "
             <<a.mt<<" "<<a.r<<" "<<a.lum<<"\n";
    // A compact remnant must have a radius consistent with its new mass.
    a.initial(10.,14); b.initial(5.,14);
    oa=StarParameterOut{}; ob=StarParameterOut{};
    manager.evolveStar(a,oa,0.); manager.evolveStar(b,ob,0.);
    semi=-10.; ecc=1.5;
    manager.merge(a,b,oa,ob,semi,ecc);
    assert(a.kw==14 && a.mt==15. && b.mt==0.);
    close(a.r,4.24e-6*a.mt); close(a.mc,a.mt); close(a.rc,a.r);
    std::cout<<"stellar merger structure regressions passed\n";
}
