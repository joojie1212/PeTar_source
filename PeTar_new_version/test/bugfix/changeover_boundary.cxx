#include "petar.hpp"
#include <vector>
#include <cmath>

// Run with an ordinary restart input; no physical integration is performed.
// Tests the real FDPS tree and isolated/connected cluster correction paths.
int main(int argc,char** argv) {
    PeTar::initialFDPS(argc,argv);
    {
        PeTar p;
        if(p.readParameters(argc,argv)<0)return 2;
        p.readDataFromFile(); p.initialParameters();
        p.domainDecompose(); p.exchangeParticle();
        const double dt=p.input_parameters.dt_soft.value;
        p.dt_manager.setStep(dt);
        p.prepareRemnantSearch(dt); p.treeNeighborSearch();
        Ptcl::group_data_mode=GroupDataMode::cm;
        p.searchCluster(); p.createGroup(dt);
#ifndef ORBIT_SAMPLING
        // Check the actual group builder (including remote hard members), not
        // just the quadrupole generator: every old cutoff must carry the
        // physical member mass before the first kick.
        auto& ap=p.hard_manager.ap_manager;
        int checked=0;
        for(int start=p.stat.n_real_loc;start<p.stat.n_all_loc;start+=ap.getArtificialParticleN()) {
            auto* arti=&p.system_soft[start];
            if(ap.getMemberN(arti)!=2)continue;
            const int cm=start+ap.getIndexOffsetCM();
            std::vector<const Ptcl*> members;
            auto collect=[&](auto& particles) {
                for(int j=0;j<particles.size();++j) {
                    const auto& q=particles[j];
                    if(q.group_data.artificial.isMember() && q.getParticleCMAddress()==cm)members.push_back(&q);
                }
            };
            collect(p.system_hard_isolated.getPtcl());
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            collect(p.system_hard_connected.getPtcl());
#endif
            assert(members.size()==2);
            auto* orbital=ap.getOrbitalParticles(arti);
            for(const auto* member:members) {
                double real_mass=0.,source_mass=0.;
                auto same_cutoff=[&](const auto& q) {
                    return q.changeover.getRin()==member->changeover.getRin() &&
                           q.changeover.getRout()==member->changeover.getRout();
                };
                for(const auto* q:members)if(same_cutoff(*q))real_mass+=q->group_data.artificial.getMassBackup();
                for(int j=0;j<ap.getOrbitalParticleN();++j)if(same_cutoff(orbital[j]))source_mass+=orbital[j].mass;
                assert(std::abs(source_mass/real_mass-1.)<1e-12);
            }
            ++checked;
        }
        assert(PS::Comm::getSum(checked)>0);
#endif
        p.treeSoftForce(); p.externalForce(); p.treeForceCorrectChangeover();
        struct CutoffTarget {double rin,rout;};
        std::vector<CutoffTarget> artificial_targets;
        std::vector<ChangeOver> artificial_old;
        int pending_orbital=0;
        for(int i=p.stat.n_real_loc;i<p.stat.n_all_loc;++i) {
            const auto& q=p.system_soft[i];
            artificial_old.push_back(q.changeover);
            artificial_targets.push_back({q.changeover.getRin()*q.changeover.r_scale_next,
                                          q.changeover.getRout()*q.changeover.r_scale_next});
            pending_orbital += q.mass>0 && q.changeover.r_scale_next!=1.0;
        }
        // This input must exercise an actual source-cutoff transition, rather
        // than only rebuilding an already synchronized representation.
        assert(PS::Comm::getSum(pending_orbital)>0);
        p.correctForceChangeOverUpdate();
        if(std::getenv("PETAR_TEST_STALE_ORBITAL")) {
            for(int i=p.stat.n_real_loc;i<p.stat.n_all_loc;++i)
                p.system_soft[i].changeover=artificial_old[i-p.stat.n_real_loc];
        }
        auto check_cutoffs=[&]() {
            for(int i=p.stat.n_real_loc;i<p.stat.n_all_loc;++i) {
                const auto& c=p.system_soft[i].changeover;
                const auto& target=artificial_targets[i-p.stat.n_real_loc];
                assert(c.r_scale_next==1.0);
                assert(std::abs(c.getRin()/target.rin-1.0)<1e-12);
                assert(std::abs(c.getRout()/target.rout-1.0)<1e-12);
            }
        };
        check_cutoffs();
        struct Cache {double soft,total; PS::F64vec acc;};
        std::vector<Cache> before;
        int members=0;
        for(int i=0;i<p.stat.n_real_loc;++i) {
            const auto& q=p.system_soft[i];
            before.push_back({q.pot_soft,q.pot_tot,q.acc});
            members+=q.group_data.artificial.isMember() && q.getParticleCMAddress()>0;
        }
        assert(PS::Comm::getSum(members)>0);
        auto check=[&]() {
            check_cutoffs();
            double error=0;
            for(int i=0;i<p.stat.n_real_loc;++i) {
                const auto& q=p.system_soft[i]; const auto& b=before[i];
                error=std::max(error,std::abs(q.pot_soft-b.soft)/std::max(1.,std::abs(b.soft)));
                error=std::max(error,std::abs(q.pot_tot-b.total)/std::max(1.,std::abs(b.total)));
                for(int k=0;k<3;++k)error=std::max(error,std::abs(q.acc[k]-b.acc[k])/std::max(1.,std::abs(b.acc[k])));
            }
            error=PS::Comm::getMaxValue(error);
            if(error>1e-11) {std::cerr<<"boundary cache error="<<error<<std::endl; PS::Abort();}
            if(PS::Comm::getRank()==0)std::cout<<"boundary cache relative error="<<error<<std::endl;
        };
        auto rebuild_and_check=[&]() {
            std::vector<std::pair<const Ptcl*,double>> masses;
            auto save=[&](auto& a,int n) {for(int i=0;i<n;++i)masses.emplace_back(&a[i],a[i].mass);};
            save(p.system_soft,p.stat.n_real_loc);
            save(p.system_hard_isolated.getPtcl(),p.system_hard_isolated.getPtcl().size());
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            save(p.system_hard_connected.getPtcl(),p.system_hard_connected.getPtcl().size());
#endif
            if(std::getenv("PETAR_TEST_UNMASKED")) {
                p.treeSoftForce(); p.externalForce(); p.treeForceCorrectChangeover();
            } else p.rebuildSoftForceWithChangeover();
            for(const auto& m:masses)assert(m.first->mass==m.second);
            check();
        };
        rebuild_and_check(); // tree-source state, before kick
        p.kick(0.0);
        rebuild_and_check(); // hard masses restored, no output
        p.writeBackHardParticles();
        rebuild_and_check(); // output has restored soft masses as well
        p.correctForceChangeOverUpdate(); // termination must precede CM conversion
        check();
#ifdef CLUSTER_VELOCITY
        p.setParticleGroupDataToCMData();
#endif
        p.system_soft.setNumberOfParticleLocal(p.stat.n_real_loc);
    }
    PeTar::finalizeFDPS();
}
