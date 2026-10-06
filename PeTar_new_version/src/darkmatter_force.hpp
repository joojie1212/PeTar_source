#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

// FDPS carries massless targets along with sources so LET construction covers
// both components. Source and target order is retained by getForce().
struct DMForce {
    PS::F64vec acc;
    PS::F64 pot;
    void clear() { acc = 0.0; pot = 0.0; }
};

struct DMTreeParticle {
    PS::F64vec pos;
    PS::F64 mass;
    bool active;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec& p) { pos = p; }
};

struct DMEPI {
    PS::F64vec pos;
    bool active;
    PS::F64vec getPos() const { return pos; }
    void copyFromFP(const DMTreeParticle& p) { pos = p.pos; active = p.active; }
};

struct DMEPJ {
    PS::F64vec pos;
    PS::F64 mass;
    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec& p) { pos = p; }
    PS::F64 getCharge() const { return mass; }
    void copyFromFP(const DMTreeParticle& p) { pos = p.pos; mass = p.mass; }
};

// FDPS's stock moment divides by mass unconditionally. Probe-only cells must
// remain finite, including on ranks containing no massive source particles.
struct DMMoment : PS::MomentMonopole {
    DMMoment() = default;
    DMMoment(const PS::MomentMonopole& m): PS::MomentMonopole(m) {}
    void set() { if (mass > 0.0) pos /= mass; else pos = 0.0; }
};

typedef PS::TreeForForce<PS::SEARCH_MODE_LONG, DMForce, DMEPI, DMEPJ,
                        DMMoment, DMMoment, PS::SPJMonopole> DMTree;

// Same softened monopole law for EP-EP and EP-SP. Include self in the kernel
// (zero acceleration) and remove -G*m/eps from DM self-potential on writeback.
struct DMGravityKernel {
    PS::F64 eps2, grav;
    DMGravityKernel(PS::F64 eps, PS::F64 g): eps2(eps*eps), grav(g) {}

    template<class Source>
    void operator()(const DMEPI* pi, const PS::S32 ni,
                    const Source* pj, const PS::S32 nj, DMForce* force) const {
#if defined(USE_SIMD) && defined(INTRINSIC_X86)
#if defined(CALC_EP_64bit)
        static thread_local PhantomGrapeQuad64Bit pg;
#else
        static thread_local PhantomGrapeQuad pg;
#endif
        // Tile both lists; FDPS group/leaf limits may exceed SIMD buffers.
        std::vector<PS::S32> indices;
        indices.reserve(ni);
        for (PS::S32 i=0; i<ni; ++i) if (pi[i].active) indices.push_back(i);
        pg.set_eps2(eps2);
        pg.set_r_crit2(0.0);
        for (std::size_t first=0; first<indices.size(); first+=pg.NIMAX) {
            const PS::S32 count=std::min(std::size_t(pg.NIMAX), indices.size()-first);
            for (PS::S32 k=0; k<count; ++k) {
                const auto x=pi[indices[first+k]].pos;
                pg.set_xi_one(k, x.x, x.y, x.z, 0.0);
            }
            // Leave one valid slot for Phantom-GRAPE's source prefetch.
            for (PS::S32 first_j=0; first_j<nj; first_j+=pg.NJMAX-1) {
                const PS::S32 count_j=std::min(PS::S32(pg.NJMAX-1), nj-first_j);
                for (PS::S32 j=0; j<count_j; ++j) {
                    const auto x=pj[first_j+j].getPos();
                    pg.set_epj_one(j, x.x, x.y, x.z, pj[first_j+j].getCharge(), 0.0);
                }
                // run_epj() masks ALL coincident pairs, including different
                // particles. The zero-cutoff kernel preserves softened pairs.
                pg.set_epj_one(count_j, 0.0, 0.0, 0.0, 0.0, 0.0);
                pg.run_epj_for_p3t_with_linear_cutoff(count, count_j);
                for (PS::S32 k=0; k<count; ++k) {
                    PS::F64 ax=0, ay=0, az=0, pot=0;
                    pg.accum_accp_one(k, ax, ay, az, pot);
                    auto& f=force[indices[first+k]];
                    f.acc += grav*PS::F64vec(ax, ay, az);
                    f.pot += grav*pot;
                }
            }
        }
#else
        for (PS::S32 i=0; i<ni; ++i) if (pi[i].active) {
            PS::F64vec acc(0.0);
            PS::F64 pot=0.0;
            for (PS::S32 j=0; j<nj; ++j) {
                const PS::F64 mass=pj[j].getCharge();
                if (mass == 0.0) continue;
                const PS::F64vec dr=pi[i].pos-pj[j].getPos();
                const PS::F64 ri=1.0/std::sqrt(dr*dr+eps2);
                acc -= (mass*ri*ri*ri)*dr;
                pot -= mass*ri;
            }
            force[i].acc += grav*acc;
            force[i].pot += grav*pot;
        }
#endif
    }
};
