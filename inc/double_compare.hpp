#ifndef DOUBLE_COMPARE_HPP
#define DOUBLE_COMPARE_HPP

namespace Compare {
    
const double epsilon = 1.0e-10;

bool is_equal(const double x, const double y);
bool is_greater_or_equal(const double x, const double y);
bool is_less_or_equal(const double x, const double y);

} // namespace Compare

#endif // DOUBLE_COMPARE_HPP