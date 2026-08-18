#pragma once

namespace BHMerger {

struct RemnantParameters {
    double vkick[3]{};
    double vkick_nor[3]{};
    double spin{};
    double spindirection[3]{};
    double spindirection_nor[3]{};
    double mfin{};
};

int initializePython();
int finalizePython();

bool calculateRemnant(double m1, double m2,
                      double dx, double dy, double dz,
                      double vrx, double vry, double vrz,
                      double spin1x, double spin1y, double spin1z,
                      double spin2x, double spin2y, double spin2z,
                      double gravitational_constant,
                      double speed_of_light);

RemnantParameters getRemnant();

} // namespace BHMerger
