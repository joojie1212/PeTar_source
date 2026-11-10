#include "petar.hpp"
#include <memory>
#include <sstream>
#include <unistd.h>
#include <random>
#include "ar_force_reference.hpp"

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
    struct Orbit { double semi=2.,ecc=0.25,period=1.; int calls=0; void calcOrbit(double) { ++calls; } } bin;
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
    assert(bin.calls==1);
    bin.semi=-1.; bin.calls=0;
    interaction.checkBHMergerDynamicalChange(bin,p1,p2,14.);
    assert(bin.calls==1); // refresh consumes the already computed orbit
    assert(p1.time_merger==-1. && p2.time_merger==-1.);
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


void testARPairKernelAndCache() {
#pragma omp parallel for
    for (int worker=0; worker<4; ++worker) {
        std::mt19937_64 local(7214+worker);
        std::uniform_real_distribution<double> draw(-1.,1.);
        ARInteraction fast; DirectedForceReference ref;
        fast.gravitational_constant=ref.gravitational_constant=0.73;
        for (int n: {1,2,3,8,32}) for(int trial=0; trial<40; ++trial) {
            std::vector<PtclHard> p(n);
            std::vector<AR::Force> f(n), r(n);
            for (int i=0;i<n;++i) {
                p[i].mass=(i==0 && trial%3==0)?0.:std::exp(draw(local)*8.);
                p[i].pos=PS::F64vec(draw(local),draw(local),draw(local));
                // Includes unequal inner radii with equal outer radii.
                p[i].changeover.setR(0.1+0.02*i, trial%2?3.:2.+0.1*i);
            }
            Float ef,er;
            Float kf=fast.calcInnerAccPotAndGTKickInv(f.data(),ef,p.data(),n);
            Float kr=ref.calcInnerAccPotAndGTKickInv(r.data(),er,p.data(),n);
            auto equal=[](Float a,Float b) {
                assert(std::abs(a-b)<=2e-14*std::max(Float(1.),std::max(std::abs(a),std::abs(b))));
            };
            equal(ef,er); equal(kf,kr);
            for(int i=0;i<n;++i) for(int k=0;k<3;++k) {
                equal(f[i].acc_in[k],r[i].acc_in[k]);
#ifdef AR_TTL
                equal(f[i].gtgrad[k],r[i].gtgrad[k]);
#endif
            }
        }
        for(double e: {0.,0.1,0.9,0.9999,0.99999}) for(double a: {1e-9,0.1,1e5}) {
            const double logf=e>PetersMergerTime::eccentricity_switch
                ? std::log(768.0/425.0)+3.5*(std::log1p(-e)+std::log1p(e))
                : std::log(PetersMergerTime::eccentricityFactor(e));
            const double expected=4.*std::log(a)+logf;
            for(int repeat=0;repeat<3;++repeat)
                assert(ARInteraction::bhMergerLogLifetimeShape(a,e)==expected);
        }
    }
}



// A resolved neighbor must reproduce direct member forces at the requested
// epoch, including backward prediction and slowdown-scaled orbital motion.
void testResolvedGroupPerturbation() {
#pragma omp parallel for
    for(int worker=0;worker<4;++worker) {
        using Group=AR::TimeTransformedSymplecticIntegrator<PtclHard,PtclH4,ARPerturber,ARInteraction,H4::ARInformation<PtclHard>>;
        AR::TimeTransformedSymplecticManager<ARInteraction> manager;
        manager.time_error_max=1e-10; manager.time_step_min=1e-12;
        manager.energy_error_relative_max=1e-9; manager.slowdown_pert_ratio_ref=1e-4;
        manager.slowdown_timescale_max=100.; manager.step_count_max=10000;
        manager.step.initialSymplecticCofficients(6);
        manager.interaction.gravitational_constant=1.; manager.interaction.eps_sq=0.;
        manager.interaction.stellar_evolution_option=0;
        manager.interaction.stellar_evolution_write_flag=false;
        Group source; source.manager=&manager;
        source.particles.setMode(COMM::ListMode::copy); source.particles.reserveMem(2);
        source.info.reserveMem(2); source.info.r_break_crit=10.;
        source.perturber.r_neighbor_crit_sq=100.;
        PtclHard origin[2];
        for(int k=0;k<2;++k) {
            origin[k].DataCopy(sampleParticle(100+k)); origin[k].mass=k+1.;
            origin[k].pos=PS::F64vec(k?0.2:-0.4,0,0);
            origin[k].vel=PS::F64vec(0,k?0.2:-0.4,0);
            origin[k].changeover.setR(0.2+0.1*k,2.+k);
            source.particles.addMemberAndAddress(origin[k]);
        }
        source.reserveIntegratorMem();
        source.info.generateBinaryTree(source.particles,1.);
        source.initialIntegration(7.);
        source.info.getBinaryTreeRoot().slowdown.setSlowDownFactor(4.);
        // Deliberately desynchronize the Hermite CM epoch from the AR epoch.
        source.particles.cm.time=5.; source.particles.cm.pos=PS::F64vec(1.,0.3,0.);
        source.particles.cm.vel=PS::F64vec(0.1,-0.2,0.05);
        for(int d=0;d<3;++d) source.particles.cm.acc0[d]=source.particles.cm.acc1[d]=0.;
        auto adr=H4::NBAdr<PtclHard>::fromGroup(source,2);
        H4::NBAdr<PtclHard> copied; copied=adr;
        ARPerturber pert; pert.reserveMem(2); pert.n_neighbor_group=1;
        pert.neighbor_address.addMember(copied);
        PtclH4 target_cm; target_cm.mass=3.; target_cm.time=7.;
        target_cm.pos=target_cm.vel=PS::F64vec(0.);
        for(int d=0;d<3;++d) target_cm.acc0[d]=target_cm.acc1[d]=0.;
        PtclHard target[2]; AR::Force force[2];
        for(int k=0;k<2;++k) {
            target[k].DataCopy(sampleParticle(200+k)); target[k].mass=k+1.;
            target[k].pos=PS::F64vec(k?0.1:-0.2,0.15,0.);
            target[k].changeover.setR(0.1,2.5);
        }
        for(double t: {6.9,7.,7.1}) {
            Float positions[6]; copied.predict_group_members(copied.group_integrator,positions,t);
            for(int k=0;k<2;++k) for(int d=0;d<3;++d)
                assert(std::abs(positions[3*k+d]-source.particles[k].pos[d]
                       -(t-7.)*source.particles[k].vel[d]/4.)<1e-14);
            manager.interaction.calcAccPert(force,target,2,target_cm,pert,t);
            PS::F64vec direct[2]={PS::F64vec(0.),PS::F64vec(0.)};
            const auto cm=source.particles.cm.pos+(t-5.)*source.particles.cm.vel;
            for(int i=0;i<2;++i) for(int j=0;j<2;++j) {
                PS::F64vec pos(positions[3*j],positions[3*j+1],positions[3*j+2]);
                auto dr=cm+pos-target[i].pos; double r=std::sqrt(dr*dr);
                direct[i]+=source.particles[j].mass*ChangeOver::calcAcc0WTwo(
                    target[i].changeover,source.particles[j].changeover,r)*dr/(r*r*r);
            }
            auto mean=(direct[0]+2.*direct[1])/3.;
            for(int i=0;i<2;++i) for(int d=0;d<3;++d) {
                assert(std::abs(force[i].acc_pert[d]-(direct[i][d]-mean[d]))<2e-13);
                assert(std::abs(force[0].acc_pert[d]+2.*force[1].acc_pert[d])<2e-13);
            }
        }
        // Removing all neighbors must not reuse the previous call's sources.
        pert.resetNeighbor(); manager.interaction.calcAccPert(force,target,2,target_cm,pert,7.);
        for(auto& f:force) for(double a:f.acc_pert) assert(a==0.);
    }
}

void testTreeSourceMassGuard() {
    Ptcl::group_data_mode=GroupDataMode::artificial;
    std::vector<FPSoft> p(3);
    for(int i=0;i<3;++i)p[i]=sampleParticle(i+1);
    p[0].group_data.artificial.setParticleTypeToMember(7.);
    p[0].setParticleCMAddress(10); p[0].mass=7.;
    p[1].group_data.artificial.setParticleTypeToMember(3.);
    p[1].mass=3.;
    const auto metadata=p[0].group_data;
    {
        TreeSourceMassGuard guard; guard.mask(p,p.size());
        assert(p[0].mass==0 && p[1].mass==3 && p[2].mass==1);
        {TreeSourceMassGuard nested; nested.mask(p,p.size());}
        assert(p[0].mass==0);
    }
    assert(p[0].mass==7 && p[0].group_data.artificial.getMassBackup()==metadata.artificial.getMassBackup());
}

void testSynchronizedPartitionCache() {
    ForceSoft::grav_const=1.; EPISoft::eps=0.; EPISoft::r_out=1.;
    struct ParticleArray: std::vector<FPSoft> {
        ParticleArray(): std::vector<FPSoft>(3) {}
        int getNumberOfParticleLocal() const { return size(); }
    } p;
    struct NeighborTree {
        EPJSoft neighbors[3];
        int getNeighborListOneParticle(const FPSoft&, EPJSoft*& out) {
            out=neighbors; return 3;
        }
    } tree;
    ArtificialParticleManager artificial;
    for(int i=0;i<3;++i) {
        p[i]=sampleParticle(i+1); p[i].mass=i+1;
        p[i].pos=PS::F64vec(i==2?32.:0.5*i,0.,0.);
    }
    auto rebuild=[&]() {
        for(int i=0;i<3;++i) {
            p[i].pot_tot=p[i].pot_soft=-p[i].mass; p[i].acc=0.;
            for(int j=0;j<3;++j)if(i!=j) {
                auto dr=p[i].pos-p[j].pos;
                double r=std::max(1.,std::sqrt(dr*dr));
                p[i].pot_tot-=p[j].mass/r; p[i].pot_soft-=p[j].mass/r;
                p[i].acc-=p[j].mass/(r*r*r)*dr;
            }
            tree.neighbors[i].copyFromFP(p[i]);
        }
        SystemHard::correctForceWithCutoffTreeNeighborOMP<ParticleArray,FPSoft,NeighborTree,EPJSoft>(p,tree,3,artificial);
    };
    rebuild(); auto before=p;
    p[0].changeover.r_scale_next=3.; p[0].changeover.updateWithRScale();
    rebuild();
    for(int i=0;i<3;++i) {
        assert(std::abs(p[i].pot_tot-before[i].pot_tot)<1e-12);
        double dh=0.; PS::F64vec ah_old(0.),ah_new(0.);
        for(int j=0;j<3;++j)if(i!=j) {
            auto dr=p[i].pos-p[j].pos; double r=std::sqrt(dr*dr);
            dh-=p[j].mass/r*(ChangeOver::calcPotWTwo(p[i].changeover,p[j].changeover,r)
                           -ChangeOver::calcPotWTwo(before[i].changeover,before[j].changeover,r));
            ah_old-=p[j].mass/(r*r*r)*ChangeOver::calcAcc0WTwo(before[i].changeover,before[j].changeover,r)*dr;
            ah_new-=p[j].mass/(r*r*r)*ChangeOver::calcAcc0WTwo(p[i].changeover,p[j].changeover,r)*dr;
        }
        assert(std::abs(p[i].pot_soft-before[i].pot_soft+dh)<1e-12);
        auto da=p[i].acc+ah_new-before[i].acc-ah_old;
        assert(std::sqrt(da*da)<1e-12);
    }
}

// Compare callback accounting to independently summed whole-system energies.
void testMergerPartitionAccounting() {
    HardManager manager;
    manager.r_in_base=0.1; manager.r_out_base=1.;
    manager.ar_manager.interaction.gravitational_constant=1.;
    manager.h4_manager.interaction.eps_sq=0.;
    HardIntegrator integrator;
    integrator.manager=&manager;
    Ptcl::mean_mass_inv=1.;
    std::vector<PtclHard> p(4);
    for(int i=0;i<4;++i) {
        p[i].DataCopy(sampleParticle(i+1));
        p[i].mass=i+1;
        p[i].pos=PS::F64vec(i==3?100.:0.4*i,0.,0.);
        p[i].changeover.setR(0.1*(i+1),1.*(i+1));
    }
    integrator.merger_initial_state=p;
    auto energy=[](const std::vector<PtclHard>& q,int kind) {
        double e=0;
        for(int i=0;i<int(q.size());++i)for(int j=0;j<i;++j) {
            if(q[i].mass==0 || q[j].mass==0)continue;
            auto dr=q[i].pos-q[j].pos; double r=std::sqrt(dr*dr);
            double w=ChangeOver::calcPotWTwo(q[i].changeover,q[j].changeover,r);
            if(kind==2 && (i==3 || j==3))continue;
#ifdef AR_CHANGEOVER
            if(kind==2) w=ChangeOver::calcPotWTwo(q[i].changeover,q[j].changeover,r);
#else
            if(kind==2) w=1.;
#endif
            e-=q[i].mass*q[j].mass/r*(kind==1?1.-w:w);
        }
        return e;
    };
    // Non-merging reconstruction: conserve pair mass and CM, but move its
    // members relative to a nearby group-external particle.
    {
        const auto before=p;
        p[3].pos=PS::F64vec(0.6,0.3,0.);
        const auto event_before=p;
        const double pending=p[0].changeover.r_scale_next;
        {
            MergerPartitionEvent event(p[0],p[1],1.,false);
            event.state=event_before; event.in_group={1,1,1,0};
            event.first_index=0; event.second_index=1;
            event.bound.owner=&integrator;
            event.bound.finish=&HardIntegrator::finishMergerPartition;
            p[0].pos.x+=0.1; p[1].pos.x-=0.05;
        } // Destructor must dispatch even though neither particle disappeared.
        const double hard=energy(p,0)-energy(event_before,0);
        const double internal=energy(p,2)-energy(event_before,2);
        const double soft=energy(p,1)-energy(event_before,1);
        assert(std::abs(hard-internal)>1e-4);
        assert(std::abs(integrator.merger_partition_correction-(hard-internal))<1e-11);
        assert(std::abs(integrator.merger_soft_correction-soft)<1e-11);
        assert(p[0].changeover.r_scale_next==pending);
        integrator.merger_partition_correction=integrator.merger_soft_correction=0;
        p=before;
    }
    double expected_hard=0, expected_soft=0;
    for(int victim: {1,2}) {
        const auto before=p;
        MergerPartitionEvent event(p[0],p[victim],1.,false);
        event.state=before; event.in_group={1,1,1,0};
        event.first_index=0; event.second_index=victim;
        event.group_potential=-0.7;
        p[0].mass+=p[victim].mass-0.1;
        p[0].pos.x+=0.02; p[victim].mass=0.;
        const double old_r=p[0].changeover.getRout();
        HardIntegrator::finishMergerPartition(&integrator,event);
        expected_hard+=energy(p,0)-energy(before,0)-(energy(p,2)-energy(before,2))-0.07;
        expected_soft+=energy(p,1)-energy(before,1);
        for(int i: {0,victim}) {
            double phi=0;
            const auto& initial=integrator.merger_initial_state;
            for(int j=0;j<4;++j)if(j!=i) {
                auto dr=initial[i].pos-initial[j].pos; double r=std::sqrt(dr*dr);
                phi-=initial[j].mass/r*(1.-ChangeOver::calcPotWTwo(initial[i].changeover,initial[j].changeover,r));
            }
            expected_soft-=(p[i].mass-before[i].mass)*phi;
        }
        assert(std::abs(integrator.merger_partition_correction-expected_hard)<1e-11);
        assert(std::abs(integrator.merger_soft_correction-expected_soft)<1e-11);
        assert(p[0].changeover.getRout()==old_r);
        assert(std::abs(old_r*p[0].changeover.r_scale_next-std::cbrt(p[0].mass))<1e-12);
    }
    integrator.manager=nullptr;
}

// Minimal active-group view for the real Hermite energy kernel (no SOFT_PERT
// in this test build). Sparse slots and non-contiguous member indices matter.
struct EnergyTestParticles: std::vector<PtclHard> {
    int getSize() const { return int(size()); }
};
struct EnergyTestGroup {
    EnergyTestParticles particles;
    struct { std::vector<int> particle_index; } info;
};
void testLocalGroupEnergy() {
    HermiteInteraction interaction;
    interaction.gravitational_constant=1.; interaction.eps_sq=0.;
    PtclHard p[6];
    for(int i=0;i<6;++i) {
        p[i].DataCopy(sampleParticle(i+1));
        p[i].vel=PS::F64vec(0.); p[i].mass=1.;
        p[i].changeover.setR(0.1,1.);
    }
    EnergyTestGroup groups[4]; int indices[2]={3,1};
    groups[3].info.particle_index={0,1};
    groups[3].particles.assign({p[0],p[1]});
    groups[3].particles[0].pos=PS::F64vec(-0.5e-10,0.,0.);
    groups[3].particles[1].pos=PS::F64vec(0.5e-10,0.,0.);
    H4::HermiteEnergy e;
    for(double center: {0.,1e-4,1.}) {
        for(int i=0;i<2;++i)p[i].pos=groups[3].particles[i].pos+PS::F64vec(center,0.,0.);
        interaction.calcEnergy(e,p,2,groups,indices,1,0);
        assert(std::abs(e.epot/(-1e10)-1.)<1e-14);
    }
    // Two active groups, one with three members, one external single.
    groups[3].info.particle_index={0,2,4};
    groups[1].info.particle_index={3,1};
    const int owner[6]={0,1,0,1,0,-1};
    PS::F64vec local[6];
    for(int i=0;i<6;++i) {
        local[i]=PS::F64vec(.03*i,.01*i*i,0.);
        p[i].pos=local[i]+PS::F64vec(owner[i]==0?1e4:owner[i]==1?1e4+0.5:1e4+2.,0.,0.);
        p[i].changeover.setR(.1+.01*i,1.+.1*i);
    }
    for(int k=0;k<2;++k) {
        auto& g=groups[indices[k]];g.particles.clear();
        for(int i:g.info.particle_index) {g.particles.push_back(p[i]);g.particles.back().pos=local[i];}
    }
    for(double eps2: {0.,.01}) {
        interaction.eps_sq=eps2;
        double expected=0.;
        for(int i=0;i<6;++i)for(int j=0;j<i;++j) {
            const auto dr=owner[i]>=0 && owner[i]==owner[j]?local[i]-local[j]:p[i].pos-p[j].pos;
            const double r=std::sqrt(dr*dr+eps2);
            expected-=ChangeOver::calcPotWTwo(p[i].changeover,p[j].changeover,r)/r;
        }
        interaction.calcEnergy(e,p,6,groups,indices,2,0);
        assert(std::abs(e.epot-expected)<1e-12);
    }
    // Dropping all groups must ignore stale scratch membership.
    double expected=0.;
    for(int i=0;i<6;++i)for(int j=0;j<i;++j) {
        const auto dr=p[i].pos-p[j].pos;const double r=std::sqrt(dr*dr+interaction.eps_sq);
        expected-=ChangeOver::calcPotWTwo(p[i].changeover,p[j].changeover,r)/r;
    }
    interaction.calcEnergy(e,p,6,groups,indices,0,0);
    assert(std::abs(e.epot-expected)<1e-12);
}

void testHardEnergyWarning() {
    const double eps=1e-3;
    // Nearly cancelling energy: retain the residual, suppress only its warning.
    HardEnergyWarning small(-0.155972,3759470.,-3759466.66394,3.33606);
    assert(small.floor_used && !small.significant(eps));
    // The real unbooked orbital event must remain detectable.
    HardEnergyWarning missed(-16576.7,38762900.,-52849400.,-14086500.);
    assert(!missed.floor_used && missed.significant(eps));
    assert(HardEnergyWarning(100.,1e6,-1e6,0.).significant(eps));
    assert(!HardEnergyWarning(0.,0.,0.,0.).significant(eps));
    assert(HardEnergyWarning(std::numeric_limits<double>::quiet_NaN(),1.,-1.,0.).significant(eps));
    assert(HardEnergyWarning(1.,std::numeric_limits<double>::infinity(),-1.,0.).significant(eps));
    // The scale remains finite even when the unscaled component sum overflows.
    HardEnergyWarning large(1.,1e308,-1e308,0.);
    assert(std::isfinite(large.denominator));
}

int main(int argc,char** argv) {
    PeTar::initialFDPS(argc,argv); omp_set_num_threads(4);
    testTreeSourceMassGuard();
    testSynchronizedPartitionCache();
    testMergerPartitionAccounting();
    testHardEnergyWarning();
    testLocalGroupEnergy();
    testARPairKernelAndCache();
    testResolvedGroupPerturbation();
    testRestartAndPrint(); testRestartClock(); testNeighborAndIdMap(); testSimdHeapBuffer();
    testEscaperReduction(); testBSEUnitsAndMass(); testOutputProtectionAndTitles();
    PeTar::finalizeFDPS(); std::cout<<"integration regressions passed\n";
}
