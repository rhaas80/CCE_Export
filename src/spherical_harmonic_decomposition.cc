#include "spherical_harmonic_decomposition.hh"

namespace CCE_export {

using std::sin, std::cos, std::pow;

// Copied from Multipole
#define idx(xx, yy)                                                            \
  (assert((xx) <= nx), assert((xx) >= 0), assert((yy) <= ny),                  \
   assert((yy) >= 0), ((xx) + (yy) * (nx + 1)))

// Copied from Multipole
CCTK_REAL Simpson2DIntegral(CCTK_REAL const *f, int nx, int ny, CCTK_REAL hx,
                            CCTK_REAL hy) {
  CCTK_REAL integrand_sum = 0;
  int ix = 0, iy = 0;

  assert(nx > 0);
  assert(ny > 0);
  assert(f);
  assert(nx % 2 == 0);
  assert(ny % 2 == 0);

  int px = nx / 2;
  int py = ny / 2;

  // Corners
  integrand_sum +=
      f[idx(0, 0)] + f[idx(nx, 0)] + f[idx(0, ny)] + f[idx(nx, ny)];

  // Edges
  for (iy = 1; iy <= py; iy++)
    integrand_sum += 4 * f[idx(0, 2 * iy - 1)] + 4 * f[idx(nx, 2 * iy - 1)];

  for (iy = 1; iy <= py - 1; iy++)
    integrand_sum += 2 * f[idx(0, 2 * iy)] + 2 * f[idx(nx, 2 * iy)];

  for (ix = 1; ix <= px; ix++)
    integrand_sum += 4 * f[idx(2 * ix - 1, 0)] + 4 * f[idx(2 * ix - 1, ny)];

  for (ix = 1; ix <= px - 1; ix++)
    integrand_sum += 2 * f[idx(2 * ix, 0)] + 2 * f[idx(2 * ix, ny)];

  // Interior
  for (iy = 1; iy <= py; iy++)
    for (ix = 1; ix <= px; ix++)
      integrand_sum += 16 * f[idx(2 * ix - 1, 2 * iy - 1)];

  for (iy = 1; iy <= py - 1; iy++)
    for (ix = 1; ix <= px; ix++)
      integrand_sum += 8 * f[idx(2 * ix - 1, 2 * iy)];

  for (iy = 1; iy <= py; iy++)
    for (ix = 1; ix <= px - 1; ix++)
      integrand_sum += 8 * f[idx(2 * ix, 2 * iy - 1)];

  for (iy = 1; iy <= py - 1; iy++)
    for (ix = 1; ix <= px - 1; ix++)
      integrand_sum += 4 * f[idx(2 * ix, 2 * iy)];

  return (1.0 / 9.0) * hx * hy * integrand_sum;
}

static inline int Sphere_Index(int it, int ip, int ntheta) {
  return it + (ntheta + 1) * ip;
}

// Copied from Multipole
void Integrate(int array_size, int ntheta, int nphi, vector<CCTK_REAL> &array1r,
               vector<CCTK_REAL> &array1i, vector<CCTK_REAL> &array2r,
               vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph, CCTK_REAL *outre,
               CCTK_REAL *outim) {

  int il = Sphere_Index(0, 0, ntheta);
  int iu = Sphere_Index(1, 0, ntheta);
  CCTK_REAL dth = th[iu] - th[il];
  iu = Sphere_Index(0, 1, ntheta);
  CCTK_REAL dph = ph[iu] - ph[il];

  // Construct an array for the real integrand
  static std::vector<CCTK_REAL> fr;
  static std::vector<CCTK_REAL> fi;
  fr.resize(array_size);
  fi.resize(array_size);

  // the below calculations take the integral of conj(array1)*array2*sin(th)
  for (int i = 0; i < array_size; i++) {
    fr[i] = (array1r[i] * array2r[i]) * sin(th[i]);
    fi[i] = (-1 * array1i[i] * array2r[i]) * sin(th[i]);
  }

  if (nphi % 2 != 0 || ntheta % 2 != 0) {
    CCTK_WARN(
        CCTK_WARN_ABORT,
        "The Simpson integration method requires even ntheta and even nphi");
  }
  *outre = Simpson2DIntegral(fr.data(), ntheta, nphi, dth, dph);
  *outim = Simpson2DIntegral(fi.data(), ntheta, nphi, dth, dph);
}

void Decompose_Spherical_Harmonics(
    vector<CCTK_REAL> &th, vector<CCTK_REAL> &phi,
    vector<CCTK_REAL> &sphere_values, vector<CCTK_REAL> &re_data,
    vector<CCTK_REAL> &im_data, vector<vector<CCTK_REAL> > &re_ylms,
    vector<vector<CCTK_REAL> > &im_ylms, int array_size, int lmax, int ntheta,
    int nphi) {
  for (int l = 0; l < lmax + 1; l++) {
    for (int m = -l; m < l + 1; m++) {
      int mode_index = l_m_to_index(l, m);
      Integrate(array_size, ntheta, nphi, re_ylms.at(mode_index),
                im_ylms.at(mode_index), sphere_values, th, phi,
                &re_data.at(mode_index), &im_data.at(mode_index));
    }
  }
}

static inline CCTK_REAL factorial(CCTK_REAL x) {
  CCTK_REAL answer = 1;

  while (x > 0) {
    answer *= x;
    x -= 1;
  }
  return answer;
}

static inline int imin(int a, int b) { return a < b ? a : b; }

static inline int imax(int a, int b) { return a > b ? a : b; }

static inline double combination(int n, int m) {
  // Binomial coefficient is undefined if these conditions do not hold
  assert(n >= 0);
  assert(m >= 0);
  assert(m <= n);

  return factorial(n) / (factorial(m) * factorial(n - m));
}

void Compute_Ylms(vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
                  vector<vector<CCTK_REAL> > &re_ylms,
                  vector<vector<CCTK_REAL> > &im_ylms, int lmax,
                  int array_size) {
  const CCTK_REAL PI = acos(-1.0);
  for (int l = 0; l < lmax + 1; l++) {
    for (int m = -l; m < l + 1; m++) {
      int ylm_index = l_m_to_index(l, m);
      for (int array_index = 0; array_index < array_size; array_index++) {
        double all_coeff = 0, sum = 0;
        all_coeff = pow(-1.0, m);
        all_coeff *= sqrt(factorial(l + m) * factorial(l - m) * (2 * l + 1) /
                          (4. * PI * factorial(l) * factorial(l)));
        sum = 0.;
        for (int i = imax(m, 0); i <= imin(l + m, l); i++) {
          double sum_coeff = combination(l, i) * combination(l, i - m);
          sum += sum_coeff * pow(-1.0, l - i) *
                 pow(cos(th[array_index] / 2.), 2 * i - m) *
                 pow(sin(th[array_index] / 2.), 2 * (l - i) + m);
        }
        re_ylms.at(ylm_index).at(array_index) =
            all_coeff * sum * cos(m * ph[array_index]);
        im_ylms.at(ylm_index).at(array_index) =
            all_coeff * sum * sin(m * ph[array_index]);
      }
    }
  }
}

} // namespace CCE_export
