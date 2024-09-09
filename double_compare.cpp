#include "double_compare.hpp"

#include <limits>
#include <cmath>

bool is_equal(double x, double y) {
    return std::fabs(x - y) < std::numeric_limits<double>::epsilon();
}