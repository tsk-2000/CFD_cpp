#include "timeint.hpp"
#include "flux.hpp"
#include <algorithm>
#include <cmath>

double compute_dt_cfl(const Grid1D& g, const std::vector<U3>& U, double CFL)
{
    double smax = 1e-300;
    for (int i = g.ilo(); i < g.ihi(); ++i) {
        auto W = cons2prim(U[i]);
        double a = std::sqrt(gamma_gas * W.p / W.r);
        smax = std::max(smax, std::abs(W.u) + a);
    }
    return CFL * g.dx / smax;
}

void step_forward_euler(const Grid1D& g, std::vector<U3>& U, std::vector<U3>& F, double dt)
{
    std::vector<U3> Un = U;
    const double c = dt / g.dx;
    for (int i = g.ilo(); i < g.ihi(); ++i) {
        U[i].r = Un[i].r - c * (F[i].r - F[i-1].r);
        U[i].ru = Un[i].ru - c * (F[i].ru - F[i-1].ru);
        U[i].E = Un[i].E - c * (F[i].E - F[i-1].E);
    }
}