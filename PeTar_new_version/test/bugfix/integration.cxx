#include "petar.hpp"
#include <memory>
#include <sstream>
#include <unistd.h>

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
    a.merger_semi=2.; a.merger_ecc=0.25;
    a.merger_record_time=12.; a.merger_dyn_log_tgw=0.01;
    a.writeAscii(file); b.time_merger=1777.25; b.writeAscii(file); rewind(file);
    read.readAscii(file); assert(read.id==1 && read.time_merger==976.5625);
    assert(read.merger_semi==2. && read.merger_ecc==0.25);
    assert(read.merger_record_time==12. && read.merger_dyn_log_tgw==0.01);
    read.readAscii(file); assert(read.id==2 && read.time_merger==1777.25);
    rewind(file);
    char line[8192]; assert(fgets(line,sizeof(line),file)); fclose(file);
    char* suffix=strstr(line,"#BHMERGER_V2"); assert(suffix); *suffix='\0';
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

void testRestartClock() {
    struct Orbit { double semi=2.,ecc=0.25,period=1.; void calcOrbit(double) {} } bin;
    FPSoft a=sampleParticle(1),b=sampleParticle(2),r1,r2;
    a.setBinaryPairID(2); b.setBinaryPairID(1);
    for(FPSoft* p: {&a,&b}) {
        p->merger_semi=bin.semi; p->merger_ecc=bin.ecc;
        p->merger_record_time=12.; p->merger_dyn_log_tgw=0.01;
    }
    FILE* f=tmpfile(); assert(f);
    a.writeAscii(f); b.writeAscii(f); rewind(f);
    r1.readAscii(f); r2.readAscii(f); fclose(f);
    PtclHard p1,p2; p1.DataCopy(r1); p2.DataCopy(r2);
    ARInteraction interaction; interaction.gravitational_constant=1.;
    interaction.checkBHMergerDynamicalChange(bin,p1,p2,13.);
    assert(p1.time_merger==a.time_merger && p2.time_merger==b.time_merger);
    assert(p1.merger_record_time==13. && p1.merger_dyn_log_tgw==0.01);
}

void testNeighborAndIdMap() {
    // Legacy SIMD counting is only translation-robust at representable scales.
    // Large-offset F32 robustness (B09) was deliberately reverted for speed.
    for (double shift: {0.,16.}) {
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
    tree->tc_glb_.resizeNoInitialize(1);
    tree->tc_glb_[0].n_ptcl_=3; tree->tc_glb_[0].adr_ptcl_=0;
    tree->spj_sorted_.resizeNoInitialize(1);
    tree->spj_sorted_[0].mass=7.; tree->spj_sorted_[0].pos=0.;
    for(int i=0;i<2;++i) {
        tree->epj_sorted_[i].mass=i+1.; tree->epj_sorted_[i].pos=0.;
    }
    tree->makeInteractionListLongForZeroTheta(PS::TagWithoutCutoff(),0);
    assert(tree->epj_for_force_[0].size()==2 && tree->spj_for_force_[0].size()==1);
    assert(tree->epj_for_force_[0][0].id==101 && tree->epj_for_force_[0][1].id==202);
    assert(tree->spj_for_force_[0][0].mass==7.);
    tree->tc_glb_[0].mom_.mass=10.; tree->tc_glb_[0].mom_.pos=0.;
    std::ostringstream check;
    PS::CheckCalcMomentLongGlobalTree(tree->epj_sorted_.getPointer(),
        tree->spj_sorted_.getPointer(),3,tree->tc_glb_.getPointer(),
        tree->tp_glb_.getPointer(),1e-12,check);
    assert(check.str().find("PASS")!=std::string::npos);

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
        p->system_soft.initialize(); p->system_soft.setNumberOfParticleLocal(256);
        p->stat.n_real_loc=256; p->n_interrupt_glb=0;
        p->stat.pcm.pos=0.; p->stat.pcm.vel=PS::F64vec(3.,0.,0.);
        p->stat.energy.clear(); p->escaper.r_escape_sq=1.; p->escaper.check_energy_flag=false;
        p->input_parameters.write_style.value=0;
        for(int i=0;i<256;++i) {
            p->system_soft[i]=sampleParticle(i+1);
            p->system_soft[i].pos.x=2.+i;
        }
        p->removeParticles();
        assert(p->stat.n_real_loc==0);
#ifdef RECORD_CM_IN_HEADER
        double expected=-256.*8.;
#else
        double expected=-256.*0.5;
#endif
        // Removing all particles also removes each pair potential once.
        for(int d=1;d<256;++d)
            expected+=p->input_parameters.gravitational_constant.value*(256-d)/double(d);
        assert(std::abs(p->stat.energy.etot_ref-expected)<1e-10);
        assert(std::abs(p->stat.energy.etot_sd_ref-expected)<1e-10);
        assert(std::abs(p->stat.energy.de_change_cum-expected)<1e-10);
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

void testOutputProtectionAndTitles() {
    const auto dir=std::filesystem::temp_directory_path()/
        ("petar-output-"+std::to_string(getpid()));
    std::filesystem::create_directories(dir);
    auto p=std::make_unique<PeTar>();
    p->system_soft.initialize(); p->system_soft.setNumberOfParticleLocal(1);
    p->system_soft[0]=sampleParticle(1);
    p->stat.n_real_loc=p->stat.n_all_loc=p->stat.n_real_glb=1;
    p->input_parameters.print_flag=false;
    p->input_parameters.fname_snp.value=(dir/"data").string();
    p->input_parameters.fname_inp.value=(dir/"alias").string();
    { std::ofstream f(dir/"data.0"); f<<"original snapshot\n"; }
    std::filesystem::create_hard_link(dir/"data.0",dir/"alias");
    p->fstatus.open(dir/"status");
    p->input_parameters.write_style.value=1; p->file_header.nfile=-1;
    p->output(true);
    std::ifstream original(dir/"data.0"); std::string line;
    std::getline(original,line); assert(line=="original snapshot");
    p->output(); assert(std::filesystem::exists(dir/"data.1"));
    p->fstatus.close(); p->fstatus.open(dir/"columns");
    p->input_parameters.write_style.value=2; p->output(); p->fstatus.close();
    std::ofstream expected(dir/"expected");
    p->stat.printColumnTitle(expected, WRITE_WIDTH);
    p->system_soft[0].printColumnTitle(expected, WRITE_WIDTH);
    std::ifstream columns(dir/"columns"); std::getline(columns,line);
    expected.close();
    std::ifstream expected_file(dir/"expected"); std::string expected_line;
    std::getline(expected_file,expected_line); assert(line==expected_line);
    p.reset(); std::filesystem::remove_all(dir);
}

int main(int argc,char** argv) {
    PeTar::initialFDPS(argc,argv); omp_set_num_threads(4);
    testRestartAndPrint(); testRestartClock(); testNeighborAndIdMap(); testSimdHeapBuffer();
    testEscaperReduction(); testBSEUnitsAndMass(); testOutputProtectionAndTitles();
    PeTar::finalizeFDPS(); std::cout<<"integration regressions passed\n";
}
