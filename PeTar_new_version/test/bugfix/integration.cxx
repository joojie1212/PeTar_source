#include "petar.hpp"
#include <memory>
#include <sstream>

FPSoft sampleParticle(int id) {
    FPSoft p;
    p.id=id; p.mass=1.; p.pos=PS::F64vec(2.,0.,0.); p.vel=PS::F64vec(1.,0.,0.);
    p.acc=0.; p.pot_tot=p.pot_soft=0.; p.n_ngb=1; p.binary_state=0;
    p.r_search=2.; p.changeover.setR(1.,2.);
    p.group_data.artificial.setParticleTypeToSingle();
#ifdef STELLAR_EVOLUTION
    p.radius=0.; p.dm=0.; p.time_record=0.; p.time_interrupt=100.;
#endif
#ifdef BSE_BASE
    p.star.initial(1.); p.star.lum=1.;
#endif
#ifdef BHMERGER
    p.spin=0.; p.time_merger=976.5625;
#endif
    return p;
}

void testRestartAndPrint() {
    FILE* file=tmpfile(); assert(file);
    FPSoft a=sampleParticle(1), b=sampleParticle(2), read;
    a.writeAscii(file); b.time_merger=1777.25; b.writeAscii(file); rewind(file);
    read.readAscii(file); assert(read.id==1 && read.time_merger==976.5625);
    read.readAscii(file); assert(read.id==2 && read.time_merger==1777.25);
    rewind(file);
    char line[8192]; assert(fgets(line,sizeof(line),file)); fclose(file);
    char* suffix=strstr(line,"#BHMERGER_V1"); assert(suffix); *suffix='\0';
    file=tmpfile(); fputs(line,file); fputc('\n',file); rewind(file);
    read.time_merger=99.; read.readAscii(file); assert(read.time_merger==-1.0); fclose(file);
    PtclHard h; h.DataCopy(a); h.id_cluster=17; h.adr_org=23;
    std::ostringstream out; h.print(out);
    assert(out.str().find("id_cluster=17")!=std::string::npos);
    H4::ParticleH4<PtclHard> h4;
    h4.acc0[0]=1; h4.acc0[1]=2; h4.acc0[2]=3;
    h4.acc1[0]=4; h4.acc1[1]=5; h4.acc1[2]=6;
    std::ostringstream h4out; h4.print(h4out);
    assert(h4out.str().find("acc1=4,5,6")!=std::string::npos);
}

void testNeighborAndIdMap() {
    for (double shift: {0.,1e8}) {
        EPISoft i; EPJSoft j;
        i.pos=PS::F64vec(shift+3.5,0,0); j.pos=PS::F64vec(shift+4.5,0,0);
        i.r_search=j.r_search=2.; j.mass=1.; j.id=12;
        ForceSoft force; force.clear();
        SearchNeighborEpEpSimd()( &i,1,&j,1,&force );
        assert(force.n_ngb==1);
    }
    auto tree=std::make_unique<PeTar::TreeForce>();
    tree->initialize(2, 0.5, 8, 64);
    tree->epj_sorted_.resizeNoInitialize(2);
    tree->epj_sorted_[0].id=101; tree->epj_sorted_[1].id=202;
    tree->tp_glb_.resizeNoInitialize(3);
    tree->tp_glb_[0].adr_ptcl_=PS::U32(1)<<31;
    tree->tp_glb_[1].adr_ptcl_=0; tree->tp_glb_[2].adr_ptcl_=1;
    assert(tree->getEpjFromId(101)==&tree->epj_sorted_[0]);
    assert(tree->getEpjFromId(202)==&tree->epj_sorted_[1]);
    assert(tree->getEpjFromId(303)==nullptr);
}

void testSimdHeapBuffer() {
    static_assert(alignof(PhantomGrapeQuad) >= 64,
                  "heap-backed SIMD workspace must retain cache-line alignment");
    EPISoft i;
    EPJSoft j;
    ForceSoft force_simd, force_scalar;
    i.pos=PS::F64vec(0.,0.,0.); i.r_search=2.; i.type=1;
    j.pos=PS::F64vec(1.,0.,0.); j.r_search=2.; j.mass=1.; j.id=2;
    force_simd.clear(); force_scalar.clear();
    EPISoft::eps=0.; EPISoft::r_out=2.;
    CalcForceEpEpWithLinearCutoffSimd()(&i,1,&j,1,&force_simd);
    CalcForceEpEpWithLinearCutoffNoSimd()(&i,1,&j,1,&force_scalar);
    assert(std::isfinite(force_simd.acc.x));
    assert(std::isfinite(force_simd.pot));
    assert(std::abs(force_simd.acc.x-force_scalar.acc.x)<1.e-6);
    assert(std::abs(force_simd.acc.y-force_scalar.acc.y)<1.e-6);
    assert(std::abs(force_simd.acc.z-force_scalar.acc.z)<1.e-6);
    assert(std::abs(force_simd.pot-force_scalar.pot)<1.e-6);
    assert(force_simd.n_ngb==force_scalar.n_ngb);
}

void testEscaperReduction() {
    for (int threads: {1,4}) {
        omp_set_num_threads(threads);
        auto p=std::make_unique<PeTar>();
        p->system_soft.initialize(); p->system_soft.setNumberOfParticleLocal(4096);
        p->stat.n_real_loc=4096; p->n_interrupt_glb=0;
        p->stat.pcm.pos=0.; p->stat.pcm.vel=PS::F64vec(3.,0.,0.);
        p->stat.energy.clear(); p->escaper.r_escape_sq=1.; p->escaper.check_energy_flag=false;
        p->input_parameters.write_style.value=0;
        for(int i=0;i<4096;++i) p->system_soft[i]=sampleParticle(i+1);
        p->removeParticles();
        assert(p->stat.n_real_loc==0);
#ifdef RECORD_CM_IN_HEADER
        const double expected=-4096.*8.;
#else
        const double expected=-4096.*0.5;
#endif
        assert(p->stat.energy.etot_ref==expected);
        assert(p->stat.energy.etot_sd_ref==expected);
        assert(p->stat.energy.de_change_cum==expected);
    }
}

void testBSEUnitsAndMass() {
    IOParamsBSE params; BSEManager manager; manager.initial(params,false);
    StarParameter a,b; a.initial(1.0,14); b.initial(0.5,14);
    a.tphys=b.tphys=100.; a.epoch=b.epoch=0.;
    double expected=-1.;
    for(double scale: {1.,2.,0.5}) {
        manager.tscale=scale;
        assert(manager.getTimeStepStar(a)*scale==500.);
        double semi=1.,ecc=0.2; int type=1;
        const double dt=manager.getTimeStepBinary(a,b,semi,ecc,type)*scale;
        if(expected<0.) expected=dt;
        assert(std::abs(dt-expected)<1e-12*std::max(1.,expected));
    }
    manager.tscale=1.;
    a.initial(1.); b.initial(1.); StarParameterOut oa,ob;
    manager.evolveStar(a,oa,0.); manager.evolveStar(b,ob,1.);
    const double old1=a.mt,old2=b.mt;
    double semi=100., period=1., ecc=0.1; BinaryEvent event{};
    int status=manager.evolveBinary(a,b,oa,ob,semi,period,ecc,event,1,0.01);
    assert(status>=0);
    assert(std::abs(oa.dm-(a.mt-old1))<1e-14);
    assert(std::abs(ob.dm-(b.mt-old2))<1e-14);
}

int main(int argc,char** argv) {
    PeTar::initialFDPS(argc,argv); omp_set_num_threads(4);
    testRestartAndPrint(); testNeighborAndIdMap(); testSimdHeapBuffer();
    testEscaperReduction(); testBSEUnitsAndMass();
    PeTar::finalizeFDPS(); std::cout<<"integration regressions passed\n";
}
