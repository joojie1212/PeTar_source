#pragma once
/* Four-target SIMD kernel for the AR perturbation force of one cluster.
 *
 * Two source kinds are handled:
 *   singles4      - point-source perturbers (NBAdr::Single)
 *   groupMembers4 - a group perturber resolved into its members, matching the
 *                   source distribution Hermite uses (upstream 2026-11 energy
 *                   batch: the AR force used to see only the group CM)
 *
 * Bitwise requirement.  The scalar code being replaced is
 *
 *   single:
 *     xi  = pi.pos[d] + xcm[d]
 *     dr  = xp[j][d] - xi
 *     r2  = dr0*dr0 + dr1*dr1 + dr2*dr2 + eps_sq      (j ascending)
 *     r   = sqrt(r2)
 *     k   = changeover weight of the larger r_out
 *     acc += ((G*m_j)/(r*r2) * k) * dr
 *
 *   group member:
 *     dr  = (xp[j][d] - xcm[d]) + (pos_k[d] - pi.pos[d])
 *     r2  = dr0*dr0 + dr1*dr1 + dr2*dr2 + eps_sq      (j ascending, k ascending)
 *     r   = sqrt(r2)
 *     w   = changeover weight of the larger r_out
 *     acc += (G*m_k*w/(r*r2)) * dr
 *
 * Every step is pinned with explicit intrinsics so the compiler cannot
 * re-associate or re-contract it, and each lane keeps one accumulator that is
 * only ever updated in source order.
 *
 * Lane packing: block of 4 targets; lane l is target i0+min(l,nb-1), i.e. the
 * last valid target is duplicated into the padding lanes and never written
 * back.  The accumulation over sources is never split or reordered.
 */
#include "Common/Float.h"
#include <algorithm>

#if defined(__AVX2__) && defined(__FMA__) && defined(__x86_64__)
#  include <immintrin.h>
/* The kernel implements the default changeover polynomial
 * (ChangeOver::calcAcc0W without INTEGRATED_CUTOFF_FUNCTION).  Builds that
 * select the integrated cutoff keep the scalar loops in ARInteraction. */
#  if !defined(INTEGRATED_CUTOFF_FUNCTION)
#    define PETAR_AR_SIMD4 1
#  endif
#endif

namespace ARSimd {

//! number of targets per vector block
static const int BLOCK = 4;

#ifdef PETAR_AR_SIMD4
static_assert(sizeof(Float) == 8, "AR 4-lane kernel requires double precision");

//! select _a where mask is set, _b elsewhere (blendv_pd takes _b as destination)
static inline __m256d choose(const __m256d _mask, const __m256d _a, const __m256d _b) {
    return _mm256_blendv_pd(_b, _a, _mask);
}

//! \f$ W_0(x) = (x-1)^4 (1 + 4x + 10x^2 + 20x^3 + 35 c x^4) \f$, x clamped to [0,1]
static inline __m256d changeoverWeight(const __m256d _r,
                                       const __m256d _rin, const __m256d _rout,
                                       const __m256d _norm, const __m256d _coff) {
    const __m256d vone  = _mm256_set1_pd(1.0);
    const __m256d vzero = _mm256_setzero_pd();
    (void)_rout; // only the caller's selection uses it; kept for signature symmetry
    __m256d u = _mm256_mul_pd(_mm256_sub_pd(_r, _rin), _norm);
    u = choose(_mm256_cmp_pd(u, vone,  _CMP_LT_OQ), u, vone);   // upper clamp
    u = choose(_mm256_cmp_pd(u, vzero, _CMP_GT_OQ), u, vzero);  // lower clamp

    const __m256d um  = _mm256_sub_pd(u, vone);
    const __m256d um2 = _mm256_mul_pd(um, um);
    const __m256d um4 = _mm256_mul_pd(um2, um2);
    const __m256d u2  = _mm256_mul_pd(u, u);
    const __m256d u3  = _mm256_mul_pd(u2, u);
    const __m256d u4  = _mm256_mul_pd(u2, u2);

    __m256d s = _mm256_fmadd_pd(_mm256_set1_pd(4.0),  u,  vone);
    s = _mm256_fmadd_pd(_mm256_set1_pd(10.0), u2, s);
    s = _mm256_fmadd_pd(_mm256_set1_pd(20.0), u3, s);
    s = _mm256_fmadd_pd(_mm256_mul_pd(_mm256_set1_pd(35.0), _coff), u4, s);
    return _mm256_mul_pd(um4, s);
}

//! r2 = dr0*dr0 then fma(dr1,dr1) then fma(dr2,dr2) then +eps (fixed order)
static inline __m256d distanceSquared(const __m256d _d0, const __m256d _d1,
                                      const __m256d _d2, const __m256d _eps) {
    __m256d t = _mm256_mul_pd(_d0, _d0);
    t = _mm256_fmadd_pd(_d1, _d1, t);
    t = _mm256_fmadd_pd(_d2, _d2, t);
    return _mm256_add_pd(t, _eps);
}
#endif

//! accumulate the single-source perturbation force for one block of 4 targets.
/*!
  @param[in] _xp:        predicted perturber positions, _xp[j][0..2]
  @param[in] _m:         perturber masses
  @param[in] _rin,_rout,_norm,_coff: cached per-single-source changeover parameters
  @param[in] _n_single:  number of single perturbers
  @param[in] _xib:       per-lane target positions relative to the c.m., [3][4]
  @param[in] _rin_i,...: per-lane target changeover parameters (padding lanes copy the last target)
  @param[in] _eps_sq,_G: softening squared and gravitational constant
  @param[out] _out:      per-lane accumulated acceleration, [3][4]
 */
inline void singles4(const Float (*_xp)[3], const Float* _m,
                     const Float* _rin, const Float* _rout,
                     const Float* _norm, const Float* _coff,
                     const int _n_single,
                     const Float _xib[3][4],
                     const Float* _rin_i, const Float* _rout_i,
                     const Float* _norm_i, const Float* _coff_i,
                     const Float _eps_sq, const Float _G,
                     Float _out[3][4]) {
#ifdef PETAR_AR_SIMD4
    const __m256d veps = _mm256_set1_pd((double)_eps_sq);

    const __m256d xi0 = _mm256_loadu_pd(&_xib[0][0]);
    const __m256d xi1 = _mm256_loadu_pd(&_xib[1][0]);
    const __m256d xi2 = _mm256_loadu_pd(&_xib[2][0]);
    const __m256d t_rout = _mm256_loadu_pd(_rout_i);
    const __m256d t_rin  = _mm256_loadu_pd(_rin_i);
    const __m256d t_norm = _mm256_loadu_pd(_norm_i);
    const __m256d t_coff = _mm256_loadu_pd(_coff_i);

    __m256d ax = _mm256_setzero_pd();
    __m256d ay = ax, az = ax;

    for (int j = 0; j < _n_single; j++) {
        const __m256d dx = _mm256_sub_pd(_mm256_set1_pd((double)_xp[j][0]), xi0);
        const __m256d dy = _mm256_sub_pd(_mm256_set1_pd((double)_xp[j][1]), xi1);
        const __m256d dz = _mm256_sub_pd(_mm256_set1_pd((double)_xp[j][2]), xi2);

        const __m256d r2 = distanceSquared(dx, dy, dz, veps);
        const __m256d r  = _mm256_sqrt_pd(r2);

        // ChangeOver::calcAcc0WTwo: strictly rout_i > rout_j keeps i, ties keep j
        const __m256d own = _mm256_cmp_pd(t_rout, _mm256_set1_pd((double)_rout[j]), _CMP_GT_OQ);
        const __m256d ci  = choose(own, t_rin,  _mm256_set1_pd((double)_rin[j]));
        const __m256d cn  = choose(own, t_norm, _mm256_set1_pd((double)_norm[j]));
        const __m256d cc  = choose(own, t_coff, _mm256_set1_pd((double)_coff[j]));
        const __m256d k   = changeoverWeight(r, ci, t_rout, cn, cc);

        const __m256d r3 = _mm256_mul_pd(r, r2);
        const __m256d gm = _mm256_set1_pd((double)(_G * _m[j]));
        const __m256d a  = _mm256_mul_pd(_mm256_div_pd(gm, r3), k);

        ax = _mm256_fmadd_pd(a, dx, ax);
        ay = _mm256_fmadd_pd(a, dy, ay);
        az = _mm256_fmadd_pd(a, dz, az);
    }
    _mm256_storeu_pd(&_out[0][0], ax);
    _mm256_storeu_pd(&_out[1][0], ay);
    _mm256_storeu_pd(&_out[2][0], az);
#else
    for (int l = 0; l < BLOCK; l++) {
        const Float xi0 = _xib[0][l], xi1 = _xib[1][l], xi2 = _xib[2][l];
        Float ax = 0.0, ay = 0.0, az = 0.0;
        for (int j = 0; j < _n_single; j++) {
            const Float dx = _xp[j][0] - xi0;
            const Float dy = _xp[j][1] - xi1;
            const Float dz = _xp[j][2] - xi2;
            const Float r2 = dx*dx + dy*dy + dz*dz + _eps_sq;
            const Float r  = std::sqrt(r2);
            const bool own = _rout_i[l] > _rout[j];
            const Float ci   = own ? _rin_i[l]  : _rin[j];
            const Float cn   = own ? _norm_i[l] : _norm[j];
            const Float cc   = own ? _coff_i[l] : _coff[j];
            Float x = (r - ci)*cn;
            x = (x < 1.0) ? x : 1.0;
            x = (x > 0.0) ? x : 0.0;
            const Float x1 = x - 1.0;
            const Float x2 = x1*x1;
            const Float x4 = x2*x2;
            const Float s2 = x*x;
            const Float s3 = s2*x;
            const Float s4 = s2*s2;
            const Float k  = x4*(1.0 + 4.0*x + 10.0*s2 + 20.0*s3 + 35.0*cc*s4);
            const Float a  = (_G*_m[j])/(r*r2)*k;
            ax += a*dx; ay += a*dy; az += a*dz;
        }
        _out[0][l] = ax; _out[1][l] = ay; _out[2][l] = az;
    }
#endif
}

//! accumulate the member-resolved group perturbation force for one block of 4 targets.
/*!
  The scalar form being reproduced is, for every group j (ascending) and every
  member k of that group (ascending):

      dr[d] = (xp[j][d] - xcm[d]) + (member_pos[d] - pi.pos[d])

  so the two subtractions stay separate; folding them into one would change the
  rounding at large common offsets.

  @param[in] _gdx:       per group, the predicted CM minus the system c.m., [g][3]
  @param[in] _gcount:    member count per group
  @param[in] _members:   per group, pointer to its member array (mass, changeover)
  @param[in] _gpos:      flat member positions gathered by the predictor, 3*n
  @param[in] _goff:      per group, offset of its members in _gpos
  @param[in] _n_group:   number of group perturbers
  @param[in] _pipos:     per-lane target positions in the cluster frame, [3][4]
  @param[in] _rin_i,...: per-lane target changeover parameters
  @param[in,out] _out:  adds to the block result already accumulated by singles4()
 */
template <class Tmember>
inline void groupMembers4(const Float (*_gdx)[3], const int* _gcount,
                          const Tmember* const* _members,
                          const Float* _gpos, const int* _goff,
                          const int _n_group,
                          const Float _pipos[3][4],
                          const Float* _rin_i, const Float* _rout_i,
                          const Float* _norm_i, const Float* _coff_i,
                          const Float _eps_sq, const Float _G,
                          Float _out[3][4]) {
#ifdef PETAR_AR_SIMD4
    const __m256d veps   = _mm256_set1_pd((double)_eps_sq);
    const __m256d pipos0 = _mm256_loadu_pd(&_pipos[0][0]);
    const __m256d pipos1 = _mm256_loadu_pd(&_pipos[1][0]);
    const __m256d pipos2 = _mm256_loadu_pd(&_pipos[2][0]);
    const __m256d t_rout = _mm256_loadu_pd(_rout_i);
    const __m256d t_rin  = _mm256_loadu_pd(_rin_i);
    const __m256d t_norm = _mm256_loadu_pd(_norm_i);
    const __m256d t_coff = _mm256_loadu_pd(_coff_i);

    // accumulate into the block result already produced by singles4()
    __m256d ax = _mm256_loadu_pd(&_out[0][0]);
    __m256d ay = _mm256_loadu_pd(&_out[1][0]);
    __m256d az = _mm256_loadu_pd(&_out[2][0]);

    for (int j = 0; j < _n_group; j++) {
        const __m256d gdx0 = _mm256_set1_pd((double)_gdx[j][0]);
        const __m256d gdx1 = _mm256_set1_pd((double)_gdx[j][1]);
        const __m256d gdx2 = _mm256_set1_pd((double)_gdx[j][2]);
        const Tmember* mem = _members[j];
        const Float* pos = _gpos + 3*_goff[j];
        for (int k = 0; k < _gcount[j]; k++) {
            if (mem[k].mass <= 0.0) continue;
            const __m256d dx = _mm256_add_pd(gdx0, _mm256_sub_pd(_mm256_set1_pd((double)pos[3*k+0]), pipos0));
            const __m256d dy = _mm256_add_pd(gdx1, _mm256_sub_pd(_mm256_set1_pd((double)pos[3*k+1]), pipos1));
            const __m256d dz = _mm256_add_pd(gdx2, _mm256_sub_pd(_mm256_set1_pd((double)pos[3*k+2]), pipos2));

            const __m256d r2 = distanceSquared(dx, dy, dz, veps);
            const __m256d r  = _mm256_sqrt_pd(r2);

            const __m256d m_rout = _mm256_set1_pd((double)mem[k].changeover.getRout());
            const __m256d own    = _mm256_cmp_pd(t_rout, m_rout, _CMP_GT_OQ);
            const __m256d ci  = choose(own, t_rin,  _mm256_set1_pd((double)mem[k].changeover.getRin()));
            const __m256d cn  = choose(own, t_norm, _mm256_set1_pd((double)mem[k].changeover.getNorm()));
            const __m256d cc  = choose(own, t_coff, _mm256_set1_pd((double)mem[k].changeover.getCoff()));
            const __m256d w   = changeoverWeight(r, ci, m_rout, cn, cc);

            const __m256d gmw = _mm256_mul_pd(_mm256_set1_pd((double)(_G*mem[k].mass)), w);
            const __m256d f   = _mm256_div_pd(gmw, _mm256_mul_pd(r, r2));

            ax = _mm256_fmadd_pd(f, dx, ax);
            ay = _mm256_fmadd_pd(f, dy, ay);
            az = _mm256_fmadd_pd(f, dz, az);
        }
    }
    _mm256_storeu_pd(&_out[0][0], ax);
    _mm256_storeu_pd(&_out[1][0], ay);
    _mm256_storeu_pd(&_out[2][0], az);
#else
    for (int l = 0; l < BLOCK; l++) {
        const Float px = _pipos[0][l], py = _pipos[1][l], pz = _pipos[2][l];
        Float ax = _out[0][l], ay = _out[1][l], az = _out[2][l];
        for (int j = 0; j < _n_group; j++) {
            const Tmember* mem = _members[j];
            const Float* pos = _gpos + 3*_goff[j];
            for (int k = 0; k < _gcount[j]; k++) {
                if (mem[k].mass <= 0.0) continue;
                const Float dx = _gdx[j][0] + (pos[3*k+0] - px);
                const Float dy = _gdx[j][1] + (pos[3*k+1] - py);
                const Float dz = _gdx[j][2] + (pos[3*k+2] - pz);
                const Float r2 = dx*dx + dy*dy + dz*dz + _eps_sq;
                const Float r  = std::sqrt(r2);
                const bool own = _rout_i[l] > mem[k].changeover.getRout();
                const Float ci   = own ? _rin_i[l]  : mem[k].changeover.getRin();
                const Float cn   = own ? _norm_i[l] : mem[k].changeover.getNorm();
                const Float cc   = own ? _coff_i[l] : mem[k].changeover.getCoff();
                Float x = (r - ci)*cn;
                x = (x < 1.0) ? x : 1.0;
                x = (x > 0.0) ? x : 0.0;
                const Float x1 = x - 1.0;
                const Float x2 = x1*x1;
                const Float x4 = x2*x2;
                const Float s2 = x*x;
                const Float s3 = s2*x;
                const Float s4 = s2*s2;
                const Float w  = x4*(1.0 + 4.0*x + 10.0*s2 + 20.0*s3 + 35.0*cc*s4);
                const Float f  = (_G*mem[k].mass)*w/(r*r2);
                ax += f*dx; ay += f*dy; az += f*dz;
            }
        }
        _out[0][l] = ax; _out[1][l] = ay; _out[2][l] = az;
    }
#endif
}

} // namespace ARSimd
