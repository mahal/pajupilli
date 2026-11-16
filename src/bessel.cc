#include "bessel.h"

#include <cmath>

#if !defined(__STDCPP_MATH_SPEC_FUNCS__) && !defined(__cpp_lib_math_special_functions)
#define PAJUPILLI_OWN_BESSEL 1
#endif

namespace pajupilli {

#ifdef PAJUPILLI_OWN_BESSEL
namespace {
constexpr double kPi = 3.14159265358979323846;

// The power series: sum_k s^k (x/2)^(2k+m) / (k! (k+m)!), s = -1 for J and +1
// for I. For I every term is positive; for J it is used only below x = 12,
// where the alternating terms stay small.
double series(int m, double x, double s) {
  const double q = 0.25 * x * x;
  double t = 1;
  for (int k = 1; k <= m; ++k) t *= 0.5 * x / k;
  double sum = t;
  for (int k = 1; k < 500; ++k) {
    t *= s * q / (double(k) * (k + m));
    sum += t;
    if (std::fabs(t) <= 1e-17 * std::fabs(sum)) break;
  }
  return sum;
}
}  // namespace

double besselJ(int m, double x) {
  if (m < 0) return (m & 1) ? -besselJ(-m, x) : besselJ(-m, x);
  if (x < 12) return series(m, x, -1);
  if (x < m) {
    // Miller: J_(k-1) = 2k/x J_k - J_(k+1) downwards from far above m (stable
    // that way), scaled at the end by 1 = J_0 + 2 (J_2 + J_4 + ...).
    const int top = 2 * ((m + 30 + int(std::sqrt(160.0 * m))) / 2);
    double above = 0, j = 1e-30, ans = 0, even = 0;
    for (int k = top; k > 0; --k) {
      const double below = 2 * k / x * j - above;
      above = j;
      j = below;  // J_(k-1)
      if (k - 1 == m) ans = j;
      if (k - 1 > 0 && (k - 1) % 2 == 0) even += j;
      if (std::fabs(j) > 1e250) j *= 1e-250, above *= 1e-250, ans *= 1e-250, even *= 1e-250;
    }
    return ans / (j + 2 * even);
  }
  // Bessel's integral, J_m(x) = 1/pi int_0^pi cos(m t - x sin t) dt, by the
  // trapezoid rule: the integrand is periodic, so it converges exponentially
  // once there are more points than x + m.
  const int n = int(x + m) + 40;
  double sum = 0.5 * (1 + ((m & 1) ? -1 : 1));  // the ends, t = 0 and pi
  for (int i = 1; i < n; ++i) {
    const double t = kPi * i / n;
    sum += std::cos(m * t - x * std::sin(t));
  }
  return sum / n;
}

double besselI(int m, double x) { return series(m < 0 ? -m : m, x, 1); }

#else

double besselJ(int m, double x) { return std::cyl_bessel_j(double(m), x); }
double besselI(int m, double x) { return std::cyl_bessel_i(double(m), x); }

#endif

}  // namespace pajupilli
