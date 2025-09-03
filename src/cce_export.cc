#include "cce_export.hh"
#include "interpolate.hh"
#include "h5_export.hh"
#include "spherical_harmonic_decomposition.hh"
#include <vector>
#include <sys/stat.h>
#include <iomanip>
#include <map>

namespace CCE_export {

using std::vector;
using std::string;
using std::ostringstream;
using std::ios;

void CCE_Export(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (out_every == 0 || cctk_iteration % out_every != 0)
    return;

  const int ntheta = 120;
  const int nphi = 240;
  const int array_size = (ntheta + 1) * (nphi + 1);
  const CCTK_REAL PI = acos(-1.0);
  const int lmax = 8;
  const int mode_count = l_m_to_index(lmax, lmax) + 1;

  // Variables on the spheres
  // extrinsic curvature, 3d vector (3, 3, array_size)
  vector<vector<vector<CCTK_REAL> > > k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dx_k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dy_k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dz_k(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  // metric, 3d vector (3, 3, array_size)
  vector<vector<vector<CCTK_REAL> > > g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dx_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dy_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dz_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dr_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  vector<vector<vector<CCTK_REAL> > > dt_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(array_size)));
  // shift (beta), 2d vector (3, array_size)
  vector<vector<CCTK_REAL> > beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dx_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dy_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dz_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dr_beta(3, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > dt_beta(3, vector<CCTK_REAL>(array_size));
  // lapse (alpha), vector (array_size)
  vector<CCTK_REAL> alpha(array_size);
  vector<CCTK_REAL> dx_alpha(array_size);
  vector<CCTK_REAL> dy_alpha(array_size);
  vector<CCTK_REAL> dz_alpha(array_size);
  vector<CCTK_REAL> dr_alpha(array_size);
  vector<CCTK_REAL> dt_alpha(array_size);

  // theta and phi points on the sphere
  vector<CCTK_REAL> th(array_size);
  vector<CCTK_REAL> ph(array_size);
  // x, y, z unit vectors toward theta and phi points on the sphere
  vector<CCTK_REAL> xhat(array_size);
  vector<CCTK_REAL> yhat(array_size);
  vector<CCTK_REAL> zhat(array_size);
  // x, y, z values of points on the sphere
  vector<CCTK_REAL> xs(array_size);
  vector<CCTK_REAL> ys(array_size);
  vector<CCTK_REAL> zs(array_size);

  // Spherical harmonic values (Y_lm), vector (mode_count)
  vector<vector<CCTK_REAL> > re_ylms(mode_count, vector<CCTK_REAL>(array_size));
  vector<vector<CCTK_REAL> > im_ylms(mode_count, vector<CCTK_REAL>(array_size));

  // Variables decomposed into spherical harmonics
  // metric, 3d vector (3, 3, mode_count)
  vector<vector<vector<CCTK_REAL> > > re_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
  vector<vector<vector<CCTK_REAL> > > im_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
  vector<vector<vector<CCTK_REAL> > > re_dr_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
  vector<vector<vector<CCTK_REAL> > > im_dr_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
  vector<vector<vector<CCTK_REAL> > > re_dt_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
  vector<vector<vector<CCTK_REAL> > > im_dt_g(
      3, vector<vector<CCTK_REAL> >(3, vector<CCTK_REAL>(mode_count)));
  // shift (beta), 2d vector (3, mode_count)
  vector<vector<CCTK_REAL> > re_beta(3, vector<CCTK_REAL>(mode_count));
  vector<vector<CCTK_REAL> > im_beta(3, vector<CCTK_REAL>(mode_count));
  vector<vector<CCTK_REAL> > re_dr_beta(3, vector<CCTK_REAL>(mode_count));
  vector<vector<CCTK_REAL> > im_dr_beta(3, vector<CCTK_REAL>(mode_count));
  vector<vector<CCTK_REAL> > re_dt_beta(3, vector<CCTK_REAL>(mode_count));
  vector<vector<CCTK_REAL> > im_dt_beta(3, vector<CCTK_REAL>(mode_count));
  // lapse (alpha), vector (mode_count)
  vector<CCTK_REAL> re_alpha(mode_count);
  vector<CCTK_REAL> im_alpha(mode_count);
  vector<CCTK_REAL> re_dr_alpha(mode_count);
  vector<CCTK_REAL> im_dr_alpha(mode_count);
  vector<CCTK_REAL> re_dt_alpha(mode_count);
  vector<CCTK_REAL> im_dt_alpha(mode_count);

  // Compute the theta and phi points as well as the corresponding x, y, z unit
  // vectors Based on the number of theta and phi points desired (ntheta, nphi)
  for (int theta_index = 0; theta_index <= ntheta; theta_index++) {
    for (int phi_index = 0; phi_index <= nphi; phi_index++) {
      const int array_index = theta_index + (ntheta + 1) * phi_index;

      th.at(array_index) = theta_index * PI / (ntheta);
      ph.at(array_index) = phi_index * 2 * PI / nphi;
      xhat.at(array_index) = sin(th.at(array_index)) * cos(ph.at(array_index));
      yhat.at(array_index) = sin(th.at(array_index)) * sin(ph.at(array_index));
      zhat.at(array_index) = cos(th.at(array_index));
    }
  }

  // loop through the desired radii
  for (int r = 0; r < nradii; r++) {

    // Extract the metric, shift, and lapse data on sphere of desired radius
    Extract_Metric_Shift_Lapse_On_Sphere(
        CCTK_PASS_CTOC, k, dx_k, dy_k, dz_k, g, dx_g, dy_g, dz_g, dr_g, dt_g,
        beta, dx_beta, dy_beta, dz_beta, dr_beta, dt_beta, alpha, dx_alpha,
        dy_alpha, dz_alpha, dr_alpha, dt_alpha, th, ph, xhat, yhat, zhat, xs,
        ys, zs, ntheta, nphi, array_size, r);

    // Decompose into spherical harmonics
    Compute_Ylms(th, ph, re_ylms, im_ylms, lmax, array_size);

    // Decompose g, dr_g, dt_g
    for (int i = 0; i < 3; i++) {
      for (int j = i; j < 3; j++) {
        Decompose_Spherical_Harmonics(th, ph, g.at(i).at(j), re_g.at(i).at(j),
                                      im_g.at(i).at(j), re_ylms, im_ylms,
                                      array_size, lmax, ntheta, nphi);
        Decompose_Spherical_Harmonics(
            th, ph, dr_g.at(i).at(j), re_dr_g.at(i).at(j), im_dr_g.at(i).at(j),
            re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
        Decompose_Spherical_Harmonics(
            th, ph, dt_g.at(i).at(j), re_dt_g.at(i).at(j), im_dt_g.at(i).at(j),
            re_ylms, im_ylms, array_size, lmax, ntheta, nphi);
      }
    }

    // Decompose beta, dr_beta, dt_beta
    for (int i = 0; i < 3; i++) {
      Decompose_Spherical_Harmonics(th, ph, beta.at(i), re_beta.at(i),
                                    im_beta.at(i), re_ylms, im_ylms, array_size,
                                    lmax, ntheta, nphi);
      Decompose_Spherical_Harmonics(th, ph, dr_beta.at(i), re_dr_beta.at(i),
                                    im_dr_beta.at(i), re_ylms, im_ylms,
                                    array_size, lmax, ntheta, nphi);
      Decompose_Spherical_Harmonics(th, ph, dt_beta.at(i), re_dt_beta.at(i),
                                    im_dt_beta.at(i), re_ylms, im_ylms,
                                    array_size, lmax, ntheta, nphi);
    }

    // Decompose alpha, dr_alpha, dt_alpha
    Decompose_Spherical_Harmonics(th, ph, alpha, re_alpha, im_alpha, re_ylms,
                                  im_ylms, array_size, lmax, ntheta, nphi);
    Decompose_Spherical_Harmonics(th, ph, dr_alpha, re_dr_alpha, im_dr_alpha,
                                  re_ylms, im_ylms, array_size, lmax, ntheta,
                                  nphi);
    Decompose_Spherical_Harmonics(th, ph, dt_alpha, re_dt_alpha, im_dt_alpha,
                                  re_ylms, im_ylms, array_size, lmax, ntheta,
                                  nphi);

    // Store output in h5 file
    if (CCTK_MyProc(cctkGH) == 0) {
      Output_Decomposed_Metric_Data(
          CCTK_PASS_CTOC, re_g, im_g, re_dr_g, im_dr_g, re_dt_g, im_dt_g,
          re_beta, im_beta, re_dr_beta, im_dr_beta, re_dt_beta, im_dt_beta,
          re_alpha, im_alpha, re_dr_alpha, im_dr_alpha, re_dt_alpha,
          im_dt_alpha, radius[r], lmax);
    }
  }
}
} // namespace CCE_export
