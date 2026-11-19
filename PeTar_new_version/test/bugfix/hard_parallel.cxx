// Reuse the production particle fixtures and AR/BSE reference tests.
#define PETAR_BUGFIX_NO_MAIN
#include "integration.cxx"

void initializeH4Fields(PtclH4& p) {
    p.time=p.dt=0.;
    for(int d=0;d<3;++d) p.acc0[d]=p.acc1[d]=0.;
}

void testMergeLevels() {
    std::mt19937 rng(5291);
    for (int n: {0,1,3,31,256,1025}) {
        PS::ReallocatableArray<int> values;
        values.resizeNoInitialize(n+4);
        std::vector<int> reference(n);
        values[0]=values[1]=-9876; values[n+2]=values[n+3]=9876;
        for(int i=0;i<n;++i) values[i+2]=reference[i]=rng()%37;
        std::sort(reference.begin(),reference.end());
        PS::MergeSortOmp(values,2,n+2,std::less<int>());
        for(int i=0;i<n;++i) assert(values[i+2]==reference[i]);
        assert(values[0]==-9876 && values[1]==-9876);
        assert(values[n+2]==9876 && values[n+3]==9876);
    }
}

void testNeighborSnapshot() {
#ifdef HARD_PARALLEL
    PtclH4 p; p.DataCopy(sampleParticle(1)); initializeH4Fields(p);
    H4::Neighbor<PtclHard> source, target;
    target.reserveMem(4); source.r_neighbor_crit_sq=target.r_neighbor_crit_sq=4.;
    // Both directions matter: ignore the live flag, including when it has
    // already been changed by another target in this force evaluation.
    for(bool saved: {false,true}) {
        target.clearNoFreeMem(); source.initial_step_flag=!saved;
        target.checkAndAddNeighborSingle(1.,p,source,0,saved);
        assert(target.initial_step_flag==saved);
    }
#endif
}

// Match the original serial drift/evolution order, including real SSE mass
// loss and stellar-type transitions, rather than just testing a no-event step.
void testStellarDrift(std::ostream& digest) {
    Ptcl::search_factor=1.; Ptcl::r_search_min=0.01;
    for(int evolution: {0,1}) {
        HardManager manager;
        IOParamsBSE params;
        auto& interaction=manager.ar_manager.interaction;
        interaction.bse_manager.initial(params,false);
        interaction.stellar_evolution_option=evolution;
        interaction.stellar_evolution_write_flag=false;
        SystemHard system; system.manager=&manager;
        system.setTimeOrigin(0.); system.initializeForOneCluster(12);
        std::vector<PtclH4> reference(12);
        for(int i=0;i<12;++i) {
            auto p=sampleParticle(i+1);
            p.star.initial(1.+(i%3));
            StarParameterOut out;
            interaction.bse_manager.evolveStar(p.star,out,0.);
            p.mass=interaction.bse_manager.getMass(p.star);
            p.time_record=p.time_interrupt=0.;
            p.pos=PS::F64vec(i,0.2*i,-0.1*i);
            p.vel=PS::F64vec(0.01*i,0.02,-0.01);
            reference[i].DataCopy(p); system.getPtcl()[i].DataCopy(p);
            initializeH4Fields(reference[i]); initializeH4Fields(system.getPtcl()[i]);
        }
        const double dt=1000.;
        double expected_energy=0.; int changed=0, changed_type=0;
        for(auto& p:reference) {
            p.pos+=p.vel*dt;
            const auto before=p;
            if(interaction.modifyOneParticle(p,0.,dt))
                expected_energy+=0.5*(p.mass*(p.vel*p.vel)-before.mass*(before.vel*before.vel));
            p.Ptcl::calcRSearch(dt);
            changed+=(p.mass!=before.mass);
            changed_type+=(p.star.kw!=before.star.kw);
        }
        if(evolution) assert(changed>0 && changed_type>0);
        system.driveForOneClusterOMP(dt);
        assert(system.getTimeOrigin()==dt);
        assert(system.energy.de_change_cum==expected_energy);
        assert(system.energy.de_sd_change_cum==expected_energy);
        for(int i=0;i<12;++i) {
            const auto& p=system.getPtcl()[i]; const auto& r=reference[i];
            assert(p.mass==r.mass && p.star.kw==r.star.kw && p.dm==r.dm);
            assert(p.time_record==r.time_record && p.time_interrupt==r.time_interrupt);
            assert(p.r_search==r.r_search);
            for(int d=0;d<3;++d) assert(p.pos[d]==r.pos[d] && p.vel[d]==r.vel[d]);
            digest<<p.mass<<' '<<p.star.kw<<' '<<p.r_search<<' '<<p.time_record<<'\n';
        }
        digest<<expected_energy<<'\n';
    }
}

void testHermiteTasks(std::ostream& digest) {
    HardIntegrator hard;
    HardManager manager;
    manager.setEpsSq(0.); manager.setGravitationalConstant(1.);
    manager.setDtRange(0.125,30);
    manager.h4_manager.step.eta_2nd=0.01; manager.h4_manager.step.eta_4th=0.01;
    auto& ar=manager.ar_manager;
    ar.energy_error_relative_max=1e-9; ar.slowdown_pert_ratio_ref=1e-4;
    ar.slowdown_timescale_max=1.; ar.step_count_max=10000;
    ar.step.initialSymplecticCofficients(6);
    ar.interaction.stellar_evolution_option=0; ar.interaction.stellar_evolution_write_flag=false;
    auto& h=hard.h4_int;
    h.manager=&manager.h4_manager; h.ar_manager=&ar; h.step=manager.h4_manager.step;
    h.particles.setMode(COMM::ListMode::local); h.particles.reserveMem(32);
    for(int i=0;i<32;++i) {
        PtclH4 p; p.DataCopy(sampleParticle(i+1)); initializeH4Fields(p);
        p.pos=PS::F64vec(0.01*i,0.03*(i%3),0.02*(i%5));
        p.vel=PS::F64vec(0.001*i,0.01,-0.02);
        h.particles.addMember(p);
    }
    h.groups.setMode(COMM::ListMode::local); h.groups.reserveMem(32);
    h.reserveIntegratorMem(); h.initialSystemSingle(0.);
#pragma omp parallel
    {
#pragma omp single
        h.initialIntegration();
    }
    for(int i=0;i<32;++i) {
        const auto& p=h.particles[i];
        H4::ForceH4 reference; reference.clear();
        for(int j=0;j<32;++j) if(j!=i)
            manager.h4_manager.interaction.calcAccJerkPairSingleSingle(reference,p,h.particles[j]);
        for(int d=0;d<3;++d) {
            assert(p.acc0[d]==reference.acc0[d] && p.acc1[d]==reference.acc1[d]);
            digest<<p.acc0[d]<<' '<<p.acc1[d]<<' ';
        }
        assert(h.neighbors[i].n_neighbor_single==31);
        digest<<p.dt<<'\n';
    }
}


void testClusterCorrection(std::ostream& digest) {
    constexpr int n=64;
    PS::ParticleSystem<FPSoft> particles; particles.initialize(); particles.setNumberOfParticleLocal(n);
    PS::ReallocatableArray<PS::S32> addresses, counts;
    addresses.resizeNoInitialize(n); counts.resizeNoInitialize(1); counts[0]=n;
    Ptcl::group_data_mode=GroupDataMode::artificial; Ptcl::mean_mass_inv=1.;
    ForceSoft::grav_const=1.; EPISoft::eps=0.; EPISoft::r_out=4.;
    for(int i=0;i<n;++i) {
        particles[i]=sampleParticle(i+1); particles[i].pos=PS::F64vec(0.2*i,0.3*(i%3),0.4*(i%5));
        particles[i].changeover.setR(0.1,4.); addresses[i]=i;
    }
    HardManager manager; manager.r_in_base=0.1; manager.r_out_base=4.;
    manager.setGravitationalConstant(1.);
    SystemHard hard; hard.manager=&manager;
    hard.setPtclForIsolatedMultiClusterOMP(particles,addresses,counts);
    hard.findGroupsAndCreateArtificialParticlesOMP<PS::ParticleSystem<FPSoft>,FPSoft>(particles,0.01);
    assert(hard.getGroupNumberOfMemberList()[0]==0);
    hard.correctForceWithCutoffClusterOMP(particles);
    for(int i=0;i<n;++i) {
        assert(std::isfinite(particles[i].pot_tot));
        for(int d=0;d<3;++d) digest<<particles[i].acc[d]<<' ';
        digest<<particles[i].pot_tot<<' '<<particles[i].pot_soft<<'\n';
    }
}

void testSinglePerturbers() {
    ARInteraction interaction; interaction.gravitational_constant=0.73; interaction.eps_sq=0.001;
    PtclH4 source[7], cm; cm.DataCopy(sampleParticle(1000));
    cm.pos=cm.vel=PS::F64vec(0.); initializeH4Fields(cm);
    ARPerturber pert; pert.reserveMem(7); pert.n_neighbor_single=7;
    for(int j=0;j<7;++j) {
        source[j].DataCopy(sampleParticle(100+j)); initializeH4Fields(source[j]);
        source[j].pos=PS::F64vec(0.1+0.3*j,0.05*j,-0.2*j);
        source[j].vel=0.; source[j].mass=0.5+j;
        source[j].changeover.setR(0.2+0.1*j,2.+0.1*j);
        for(int d=0;d<3;++d) source[j].acc0[d]=source[j].acc1[d]=0.;
        pert.neighbor_address.addMember(H4::NBAdr<PtclHard>(&source[j],j));
    }
    for(int n: {1,2,3,4,5,8,9}) {
        std::vector<PtclHard> target(n); std::vector<AR::Force> force(n);
        std::vector<PS::F64vec> reference(n,PS::F64vec(0.));
        PS::F64vec mean(0.); double mass=0.;
        for(int i=0;i<n;++i) {
            target[i].DataCopy(sampleParticle(i+1)); target[i].mass=1.+i;
            target[i].pos=PS::F64vec(-0.1*i,0.13*i,-0.04*i);
            target[i].changeover.setR(0.1,2.);
            for(int j=0;j<7;++j) {
                auto dr=source[j].pos-target[i].pos;
                double r2=dr*dr+interaction.eps_sq, r=std::sqrt(r2);
                const double w=ChangeOver::calcAcc0WTwo(target[i].changeover,source[j].changeover,r);
                reference[i]+=(interaction.gravitational_constant*source[j].mass/(r*r2)*w)*dr;
            }
            mean+=target[i].mass*reference[i]; mass+=target[i].mass;
        }
        mean/=mass; cm.mass=mass;
        interaction.calcAccPert(force.data(),target.data(),n,cm,pert,0.);
        for(int i=0;i<n;++i) for(int d=0;d<3;++d) {
            const double expected=reference[i][d]-mean[d];
            assert(std::abs(force[i].acc_pert[d]-expected)<2e-13*std::max(1.,std::abs(expected)));
        }
    }
}

int main(int argc,char** argv) {
    PeTar::initialFDPS(argc,argv);
    const int threads=argc>1 ? std::atoi(argv[1]) : 1;
    omp_set_dynamic(0); omp_set_num_threads(threads);
    std::ofstream digest(argc>2 ? argv[2] : "/tmp/petar-hard-parallel.digest");
    digest<<std::hexfloat;
    testMergeLevels(); testNeighborSnapshot();
    testStellarDrift(digest); testHermiteTasks(digest); testClusterCorrection(digest);
    testSinglePerturbers();
    testARPairKernelAndCache(); testResolvedGroupPerturbation(); testBSEUnitsAndMass();
    PeTar::finalizeFDPS();
    std::cout<<"hard parallel regressions passed (threads="<<threads<<")\n";
}
