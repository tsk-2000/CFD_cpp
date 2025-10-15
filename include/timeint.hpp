#pragma once
#include <vector>
#include "grid.hpp"
#include "state.hpp"

double compute_dt_cfl(const Grid1D& g, const std::vector<U3>& U, double CFL);

void step_forward_euler(const Grid1D& g, std::vector<U3>& U, std::vector<U3>& F, double dt);