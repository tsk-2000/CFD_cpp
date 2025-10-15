#pragma once
#include <cmath>
struct U3 { double r, ru, E; };
struct P3 { double r, u, p; };
inline constexpr double gamma_gas = 1.4;

inline P3 cons2prim(const U3& U) {
    double u = U.ru / U.r;
    double p = (gamma_gas - 1.0) * (U.E - 0.5 * U.r * u * u);
    return {U.r, u, p};
}
inline U3 prim2cons(const P3& W){
    double E = W.p/(gamma_gas - 1.0) + 0.5 * W.r * W.u * W.u;
    return { W.r, W.r * W.u, E };
}