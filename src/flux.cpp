#include "flux.hpp"
#include <algorithm>
#include <cmath>

U3 flux_phys(const U3& U) {
    P3 W = cons2prim(U);
    double H = (U. E + W.p) / W.r;
    return { W.r*W.u, W.r*W.u*W.u + W.p, W.r*H*W.u };
}

U3 flux_rusanov(const U3& UL, const U3& UR) {
    P3 WL = cons2prim(UL);
    P3 WR = cons2prim(UR);

    double aL = std::sqrt(gamma_gas * WL.p / WL.r);
    double aR = std::sqrt(gamma_gas * WR.p / WR.r);
    double smax = std::max(std::abs(WL.u) + aL, std::abs(WR.u) + aR);

    U3 FL = flux_phys(UL);
    U3 FR = flux_phys(UR);

    return {
        0.5 * (FL.r + FR.r) - 0.5 * smax * (UR.r - UL.r),
        0.5 * (FL.ru + FR.ru) - 0.5 * smax * (UR.ru - UL.ru),
        0.5 * (FL.E + FR.E) - 0.5 * smax * (UR.E - UL.E)
    };
}