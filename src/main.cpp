#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include "state.hpp"
#include "grid.hpp"
#include "flux.hpp"
#include "timeint.hpp"
// #include <vector>
// #include <cblas.h>   // BLAS
// #include <lapacke.h> // LAPACK (Cインターフェース)

int main() {

    auto g = Grid1D::make(2001, 1, -0.5, 0.5);

    std::vector<U3> U(g.Ntot);

    for (int i = g.begin(); i < g.end(); ++i) {
        if (g.xc[i] < 0.0) U[i] = prim2cons({1.0, 0.0, 1.0});
        else               U[i] = prim2cons({0.125, 0.0, 0.1});
    }


    const double CFL = 0.5;
    const double t_end = 1.0;
    double t = 0.0;
    int iteration = 0;
    const int dump_every = 1;
    const std::string h5file = "./debug/output/out.h5";

    for (int i = g.ilo(); i < g.ihi(); ++i) {
        std::string output = "./debug/output/output_" + std::to_string(i) + ".csv";
        std::ofstream file(output);
        auto W = cons2prim(U[i]);
        file << "t,r,u,p\n";
        file.close();
    }

    while (t < t_end) {
        g.apply_neumann(U);

        std::vector<U3> F(g.Ntot - 1);
        for (int f = g.begin(); f < g.end() - 1; ++f) {
            F[f] = flux_rusanov(U[f], U[f+1]);
        }

        double dt = compute_dt_cfl(g, U, CFL);
        if (t + dt > t_end) dt = t_end - t;

        step_forward_euler(g, U, F, dt);
        t += dt;
        iteration += 1;

        if (iteration % dump_every == 0) {
            for (int i = g.ilo(); i < g.ihi(); ++i) {
                std::string output = "./debug/output/output_" + std::to_string(i) + ".csv";
                std::ofstream file(output, std::ios::app);
                auto W = cons2prim(U[i]);
                // std::printf("%.8f, %.8f, %.8f, %.8f\n", t, W.r, W.u, W.p);
                file << t << "," << W.r << "," << W.u << "," << W.p << "\n";
                file.close();
            }
        }
    }

    std::printf("# x, rho, u, p\n");
    for (int i = g.ilo(); i < g.ihi(); ++i) {
        auto W = cons2prim(U[i]);
        std::printf("%.8f, %.8f, %.8f, %.8f\n", g.xc[i], W.r, W.u, W.p);
    }
    std::printf("# iteration = %d\n", iteration);
    
    // std::string output = "./debug/output.csv";
    // std::ofstream ofs(output);
    // ofs << "# x, rho, u, p\n";
    // ofs << std::fixed << std::setprecision(8);

    // for (int i = g.ilo(); i < g.ihi(); ++i) {
    //     auto W = cons2prim(U[i]);
    //     ofs << g.xc[i] << ", " << W.r << ", " << W.u << ", " << W.p << "\n";
    // }

    // ofs.close();

    // for (int i : {g.ilo(), g.ilo() + 50, g.ilo()+100, g.ihi()-1}) {
    //     auto 
    //     W = cons2prim(U[i]);
    //     std::printf("i = %d x = %.4f rho = %.4f\n", i, g.xc[i], W.r);
    // }
    return 0;
}
