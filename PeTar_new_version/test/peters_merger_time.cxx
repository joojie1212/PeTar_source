#include <algorithm>
#include <cmath>
#include <iostream>

#include "peters_merger_time.hpp"

int main() {
    struct Reference { double eccentricity, factor; };
    const Reference references[] = {
        {0.0,    1.0},
        {0.1,    0.9640015943645863},
        {0.5,    0.3577982552068474},
        {0.9,    0.003407024596059112},
        {0.99,   1.6685681838299192e-6},
        {0.9999, 2.0028669646730718e-13},
    };
    double maximum_relative_error = 0.0;
    for (const Reference& reference: references) {
        const double factor = PetersMergerTime::eccentricityFactor(
            reference.eccentricity);
        const double relative_error = std::abs(factor/reference.factor-1.0);
        maximum_relative_error = std::max(maximum_relative_error, relative_error);
        if (!std::isfinite(factor) || relative_error > 0.002) return 1;
    }

    const double at_switch = PetersMergerTime::eccentricityFactor(0.9999);
    const double above_switch = PetersMergerTime::eccentricityFactor(
        std::nextafter(0.9999, 1.0));
    const double switch_jump = above_switch/at_switch-1.0;
    if (!(switch_jump > 0.02 && switch_jump < 0.021)) return 2;
    if (!std::isnan(PetersMergerTime::eccentricityFactor(-0.1))
        || !std::isnan(PetersMergerTime::eccentricityFactor(1.0))) return 3;

    std::cout << "Peters table maximum sampled relative error: "
              << maximum_relative_error << '\n'
              << "asymptotic switch jump: " << switch_jump << '\n';
    return 0;
}
