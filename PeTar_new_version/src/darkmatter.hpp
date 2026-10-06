#pragma once

#ifdef DARKMATTER

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>
#include "darkmatter_force.hpp"

//! DM snapshots have their own header because softening is a component-wide
//! restart parameter rather than a per-particle property.
class DarkMatterFileHeader {
public:
    long long int nfile;
    long long int n_body;
    double time;
    double eps;

    DarkMatterFileHeader(): nfile(0), n_body(0), time(0.0), eps(0.0) {}

    int readAscii(FILE* fp) {
        const int n = fscanf(fp, "%lld %lld %lf %lf\n",
                             &nfile, &n_body, &time, &eps);
        if (n != 4 || !(eps > 0.0) || !std::isfinite(eps)) {
            std::cerr << "Error: DARKMATTER header requires: fid n time softening"
                      << std::endl;
            std::abort();
        }
        return static_cast<int>(n_body);
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp, "%lld %lld %26.17e %26.17e\n",
                nfile, n_body, time, eps);
    }

    int readBinary(FILE* fp) {
        if (fread(this, sizeof(DarkMatterFileHeader), 1, fp) != 1) {
            std::cerr << "Error: failed to read DARKMATTER header" << std::endl;
            std::abort();
        }
        if (!(eps > 0.0) || !std::isfinite(eps)) {
            std::cerr << "Error: invalid DARKMATTER softening " << eps << std::endl;
            std::abort();
        }
        return static_cast<int>(n_body);
    }

    void writeBinary(FILE* fp) const {
        if (fwrite(this, sizeof(DarkMatterFileHeader), 1, fp) != 1) {
            std::cerr << "Error: failed to write DARKMATTER header" << std::endl;
            std::abort();
        }
    }
};

class DarkMatterParticle {
public:
    PS::F64 mass;
    PS::F64vec pos;
    PS::F64vec vel;
    PS::S64 id;
    PS::F64vec acc;
    PS::F64 pot_dm;
    PS::F64 pot_star;
#ifdef GALPY
    // Recomputed before every kick, including restart. Preserve the legacy
    // snapshot prefix; DM uses the stellar header's common frame offsets.
    PS::F64 pot_ext = 0.0;
    PS::F64vec acc_ext = PS::F64vec(0.0);
#endif

    DarkMatterParticle(): mass(0.0), pos(0.0), vel(0.0), id(0),
                          acc(0.0), pot_dm(0.0), pot_star(0.0) {}

    PS::F64vec getPos() const { return pos; }
    void setPos(const PS::F64vec& p) { pos = p; }
    PS::F64 getCharge() const { return mass; }

    void clearForce() {
        acc = PS::F64vec(0.0);
        pot_dm = 0.0;
        pot_star = 0.0;
#ifdef GALPY
        pot_ext = 0.0;
        acc_ext = PS::F64vec(0.0);
#endif
    }

    void writeAscii(FILE* fp) const {
        fprintf(fp,
                "%26.17e %26.17e %26.17e %26.17e "
                "%26.17e %26.17e %26.17e %lld "
                "%26.17e %26.17e %26.17e %26.17e %26.17e\n",
                mass, pos.x, pos.y, pos.z, vel.x, vel.y, vel.z,
                static_cast<long long>(id),
                acc.x, acc.y, acc.z, pot_dm, pot_star);
    }

    void readAscii(FILE* fp) {
        long long id_tmp = 0;
        clearForce();
        const int n = fscanf(fp,
                             "%lf %lf %lf %lf %lf %lf %lf %lld "
                             "%lf %lf %lf %lf %lf",
                             &mass, &pos.x, &pos.y, &pos.z,
                             &vel.x, &vel.y, &vel.z, &id_tmp,
                             &acc.x, &acc.y, &acc.z, &pot_dm, &pot_star);
        id = id_tmp;
        if (n != 13) {
            std::cerr << "Error: a DARKMATTER snapshot row requires 13 columns; got "
                      << n << std::endl;
            std::abort();
        }
        validate();
    }

    void writeBinary(FILE* fp) const {
        if (fwrite(this, snapshotSize(), 1, fp) != 1) {
            std::cerr << "Error: failed to write DARKMATTER particle" << std::endl;
            std::abort();
        }
    }

    void readBinary(FILE* fp) {
        clearForce();
        if (fread(this, snapshotSize(), 1, fp) != 1) {
            std::cerr << "Error: failed to read DARKMATTER particle" << std::endl;
            std::abort();
        }
        validate();
    }

    void validate() const {
        if (!(mass > 0.0) || id <= 0 || !std::isfinite(mass)
            || !std::isfinite(pos.x) || !std::isfinite(pos.y) || !std::isfinite(pos.z)
            || !std::isfinite(vel.x) || !std::isfinite(vel.y) || !std::isfinite(vel.z)) {
            std::cerr << "Error: invalid DARKMATTER particle id=" << id
                      << " mass=" << mass << std::endl;
            std::abort();
        }
    }

    static std::size_t snapshotSize() {
#ifdef GALPY
        return offsetof(DarkMatterParticle, pot_ext);
#else
        return sizeof(DarkMatterParticle);
#endif
    }
};


template <class T>
inline void darkMatterAllGather(const T* local, const PS::S32 n_local,
                                std::vector<T>& global) {
    const PS::S32 n_proc = PS::Comm::getNumberOfProc();
    std::vector<PS::S32> counts(n_proc);
    std::vector<PS::S32> displs(n_proc);
    PS::Comm::allGather(&n_local, 1, counts.data());
    PS::S32 n_global = 0;
    for (PS::S32 i=0; i<n_proc; ++i) {
        displs[i] = n_global;
        n_global += counts[i];
    }
    global.resize(n_global);
    T dummy;
    T* send = n_local > 0 ? const_cast<T*>(local) : &dummy;
    T* recv = n_global > 0 ? global.data() : &dummy;
    PS::Comm::allGatherV(send, n_local, recv, counts.data(), displs.data());
}

#endif // DARKMATTER
