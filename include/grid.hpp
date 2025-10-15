#pragma once
#include <vector>
#include <cassert>
#include <algorithm>


struct Grid1D {
    int Nc;
    int ng;
    int Ntot;
    double x0, x1;
    double dx;
    std::vector<double> xc;

    static Grid1D make(int Nc_, int ng_, double x0_, double x1_) {
        assert(Nc_ > 0 && ng_ >= 0 && x1_ > x0_);
        Grid1D g;
        g.Nc = Nc_;
        g.ng = ng_;
        g.Ntot = Nc_ + 2 * ng_;
        g.x0 = x0_;
        g.x1 = x1_;
        g.dx = (x1_ - x0_) / Nc_;
        g.xc.resize(g.Ntot);

        for (int i = 0; i < g.Ntot; ++i) {
            g.xc[i] = g.x0 + ( (i - g.ng) + 0.5 ) * g.dx;
        }
        return g;
    }

    inline int idx(int i_phys) const { return i_phys + ng; }

    inline int ilo() const { return ng; }
    inline int ihi() const { return ng + Nc;}

    inline int begin() const { return 0; }
    inline int end() const {  return Ntot; }

    inline double xf(int i_face) const {
        return x0 + i_face * dx;
    }

    template< class ArrayLike >
    void apply_neumann(ArrayLike& A) const {
        A[0] = A[1];
        A[Ntot-1] = A[Ntot-2];
    }
};