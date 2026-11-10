#pragma once
#include "Common/binary_tree.h"
#include <iostream>

//! PseudoParticle Multipole method for counter force from binaries
class PseudoParticleMultipoleManager{
public:
    //! check paramters
    bool checkParams() {
        return true;
    }

    //! create pseudoparticle multipole particles
    /*! Create three quadrupole nodes, each split into two member-weighted sources.
        Keeping both sources at the same position preserves the original mass
        distribution and quadrupole. The alternating member index used to assign
        changeovers must also describe the source mass: equal M/3 sources with
        alternating old member radii give incorrect forces before the radii are
        committed. After commitment each pair is equivalent to one M/3 source.
        @param[out] _ptcl_artificial: array of getParticleN() sample particles
        @param[in] _bin: binary orbit 
     */
    template <class Tptcl>
    void createSampleParticles(Tptcl* _ptcl_artificial,
                               COMM::BinaryTree<Tptcl,COMM::Binary> &_bin) {
        PS::F64 m12 = _bin.mass;
        PS::F64 mu = _bin.m1*_bin.m2/m12;
        PS::F64 prefactor = _bin.semi*std::sqrt(mu/m12);
        PS::F64 ecc2 = _bin.ecc*_bin.ecc;

        PS::F64 beta = prefactor*std::sqrt((1-ecc2)/4.0);
        PS::F64 alpha= prefactor*std::sqrt(3*ecc2+0.75);

        PS::F64vec nodes[3] = {PS::F64vec(0, 2*beta, 0),
                              PS::F64vec(alpha, -beta, 0),
                              PS::F64vec(-alpha, -beta, 0)};
        for (int k=0; k<3; ++k) {
            _bin.rotateToOriginalFrame(&(nodes[k].x));
            nodes[k] += _bin.pos;
            for (int member=0; member<2; ++member) {
                auto& p = _ptcl_artificial[2*k+member];
                // Leaf masses have already been masked for the tree. Use the
                // physical component masses stored in the binary orbit.
                p.mass = (member==0 ? _bin.m1 : _bin.m2)/3.0;
                p.pos = nodes[k];
                p.vel = _bin.vel;
            }
        }

#ifdef ARTIFICIAL_PARTICLE_DEBUG
        PS::F64vec dv = (nodes[0] + nodes[1] + nodes[2])/3 - _bin.pos;
        assert(dv*dv<1e-10);
#endif
    }

    //! get particle number 
    static PS::S32 getParticleN() {
        return 6;
    }

    //! write class data to file with binary format
    /*! @param[in] _fp: FILE type file for output
     */
    void writeBinary(FILE *_fp) {
    }

    //! read class data to file with binary format
    /*! @param[in] _fp: FILE type file for reading
     */
    void readBinary(FILE *_fin) {
    }    

    //! print parameters
    void print(std::ostream & _fout) const{
    }
};
