#pragma once
#include "state.hpp"

U3 flux_phys(const U3& U);

U3 flux_rusanov(const U3& UL, const U3& UR);