// Bessel functions of whole-number order: J_m(x) and I_m(x), x >= 0, for the
// drum heads, the cymbals and the sci-fi sounds. The C++17 std::cyl_bessel_j
// and _i where the library has them (libstdc++), our own where it does not
// (libc++ on macOS).
#pragma once

namespace pajupilli {

double besselJ(int m, double x);
double besselI(int m, double x);

}  // namespace pajupilli
