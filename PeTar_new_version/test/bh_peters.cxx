#include "../src/bh_peters.hpp"
#include <cassert>
#include <iostream>

double invariant(BHPeters::Orbit o) {
    return o.a*(1-o.e*o.e)/std::pow(o.e,12.0/19.0)
        /std::pow(1+(121.0/304.0)*o.e*o.e,870.0/2299.0);
}
int main() {
    using namespace BHPeters;
    const double inf=std::numeric_limits<double>::infinity();
    // Independent circular solution: a(t)^4=a0^4-4 beta t.
    auto circular=advance({10,0},inf,1,1);
    assert(circular.reached_cutoff);
    assert(std::abs(circular.elapsed/2499.75-1)<2e-8);
    auto partial=advance({10,0},100,1,1);
    assert(!partial.reached_cutoff);
    assert(std::abs(partial.orbit.a/std::pow(9600.0,0.25)-1)<2e-8);
    // Peters' independently derived a(e) invariant and step convergence.
    for (double e: {0.01,0.5,0.9,0.999}) {
        Orbit initial{10,e};
        auto full=advance(initial,inf,1,0.1);
        auto fine=advance(initial,inf,1,0.1,0.005);
        assert(full.reached_cutoff);
        assert(std::abs(invariant(full.orbit)/invariant(initial)-1)<2e-6);
        assert(std::abs(full.elapsed/fine.elapsed-1)<2e-7);
        auto half=advance(initial,full.elapsed/2,1,0.1);
        auto rest=advance(half.orbit,inf,1,0.1);
        assert(std::abs((half.elapsed+rest.elapsed)/full.elapsed-1)<2e-7);
        assert(half.orbit.a<initial.a && half.orbit.e<initial.e);
    }
    // Wider perturbed orbits MUST produce a later deadline.
    assert(advance({20,0.5},inf,1,0.1).elapsed
           >advance({10,0.5},inf,1,0.1).elapsed);
    // A small pericentre alone must not be interpreted as coalescence.
    assert(!advance({10,0.999},0,1,0.1).reached_cutoff);
    bool rejected=false;
    try { advance({-10,1.1},1,1,0.1); }
    catch (const std::invalid_argument&) { rejected=true; }
    assert(rejected);
    std::cout << "Peters circular, eccentric, convergence, split-time, perturbation and domain tests passed\n";
}
