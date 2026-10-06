#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

// Orbit-averaged quadrupole evolution (Peters 1964, Phys. Rev. 136, B1224).
// This module has no particle, BSE, slowdown, or mass-loss side effects.
// All inputs must use the same units as G and c. The stopping separation
// is a caller-supplied modelling boundary, not a prediction of full GR.
namespace BHPeters {
struct Orbit {
    double a, e;
};
struct Result {
    Orbit orbit;
    double elapsed;
    bool reached_cutoff;
};

inline double beta(double G, double c, double m1, double m2) {
    if (!(G > 0 && c > 0 && m1 > 0 && m2 > 0)
        || !std::isfinite(G+c+m1+m2))
        throw std::invalid_argument("Invalid Peters masses or units");
    const double b = (64.0/5.0)*G*G*G*m1*m2*(m1+m2)/std::pow(c,5);
    if (!(b > 0) || !std::isfinite(b))
        throw std::invalid_argument("Unrepresentable Peters coefficient");
    return b;
}

// x=log(a): de/dx is regular even when inspiral accelerates dramatically.
inline double dedx(double e) {
    const double e2=e*e;
    return (19.0/12.0)*e*(1-e2)*(1+(121.0/304.0)*e2)
        /(1+(73.0/24.0)*e2+(37.0/96.0)*e2*e2);
}
inline double minus_dtdx(double a, double e, double b) {
    const double e2=e*e;
    return std::pow(a,4)*std::pow((1-e)*(1+e),3.5)
        /(b*(1+(73.0/24.0)*e2+(37.0/96.0)*e2*e2));
}
inline Result step(Orbit o, double h, double b) {
    const double k1=dedx(o.e);
    const double e2=o.e+h*k1/2;
    const double k2=dedx(e2);
    const double e3=o.e+h*k2/2;
    const double k3=dedx(e3);
    const double e4=o.e+h*k3;
    const double k4=dedx(e4);
    const double amid=o.a*std::exp(h/2), afinal=o.a*std::exp(h);
    const double dt=-h*(minus_dtdx(o.a,o.e,b)
        +2*minus_dtdx(amid,e2,b)+2*minus_dtdx(amid,e3,b)
        +minus_dtdx(afinal,e4,b))/6;
    return {{afinal,std::max(0.0,o.e+h*(k1+2*k2+2*k3+k4)/6)},dt,false};
}

// Stop only when the WHOLE osculating ellipse is within cutoff. In particular
// a tiny pericentre does not authorize deletion of a distant particle.
// dt=+infinity requests a lifetime prediction. Recompute it from the current
// perturbed orbit; never retain min(old_deadline,new_deadline).
inline Result advance(Orbit o, double dt, double b, double cutoff,
                      double log_step=0.01) {
    if (!(o.a>0 && o.e>=0 && o.e<1 && dt>=0 && b>0 && cutoff>0
          && log_step>0 && log_step<=0.05)
        || !std::isfinite(o.a+o.e+b+cutoff+log_step))
        throw std::invalid_argument("Invalid bound Peters orbit or interval");
    double elapsed=0;
    for (int n=0;n<100000;++n) {
        if (o.a*(1+o.e)<=cutoff) return {o,elapsed,true};
        if (elapsed>=dt) return {o,elapsed,false};
        double h=-log_step;
        Result next=step(o,h,b);
        const bool crosses_cutoff=next.orbit.a*(1+next.orbit.e)<=cutoff;
        const bool crosses_time=next.elapsed>dt-elapsed;
        if (crosses_cutoff || crosses_time) {
            double lo=0, hi=log_step;
            // Locate the first event in this monotonic segment.
            for (int i=0;i<60;++i) {
                const double mid=(lo+hi)/2;
                Result trial=step(o,-mid,b);
                if (trial.orbit.a*(1+trial.orbit.e)<=cutoff
                    || trial.elapsed>dt-elapsed) hi=mid;
                else lo=mid;
            }
            next=step(o,-hi,b);
            const bool reached=next.orbit.a*(1+next.orbit.e)
                <=cutoff*(1+16*std::numeric_limits<double>::epsilon());
            return {next.orbit,elapsed+next.elapsed,reached};
        }
        elapsed+=next.elapsed;
        o=next.orbit;
    }
    throw std::runtime_error("Peters integration exceeded step limit");
}
} // namespace BHPeters
