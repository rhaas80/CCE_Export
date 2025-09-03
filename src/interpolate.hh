#ifndef INTERPOLATE_HH
#define INTERPOLATE_HH

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include "vector"
#include "string"

namespace CCE_export {

using std::vector;
using std::string;

void Interpolate_On_Sphere_With_Derivatives(
    CCTK_ARGUMENTS, vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys,
    vector<CCTK_REAL> &zs, std::string name, vector<CCTK_REAL> &sphere_values,
    vector<CCTK_REAL> &sphere_dx, vector<CCTK_REAL> &sphere_dy,
    vector<CCTK_REAL> &sphere_dz, CCTK_INT array_size);

void Interpolate_On_Sphere(CCTK_ARGUMENTS, vector<CCTK_REAL> &xs,
                           vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs,
                           std::string name, vector<CCTK_REAL> &sphere_values,
                           CCTK_INT array_size);

void Extract_Metric_Shift_Lapse_On_Sphere(
    CCTK_ARGUMENTS, vector<vector<vector<CCTK_REAL> > > &k,
    vector<vector<vector<CCTK_REAL> > > &dx_k,
    vector<vector<vector<CCTK_REAL> > > &dy_k,
    vector<vector<vector<CCTK_REAL> > > &dz_k,
    vector<vector<vector<CCTK_REAL> > > &g,
    vector<vector<vector<CCTK_REAL> > > &dx_g,
    vector<vector<vector<CCTK_REAL> > > &dy_g,
    vector<vector<vector<CCTK_REAL> > > &dz_g,
    vector<vector<vector<CCTK_REAL> > > &dr_g,
    vector<vector<vector<CCTK_REAL> > > &dt_g, vector<vector<CCTK_REAL> > &beta,
    vector<vector<CCTK_REAL> > &dx_beta, vector<vector<CCTK_REAL> > &dy_beta,
    vector<vector<CCTK_REAL> > &dz_beta, vector<vector<CCTK_REAL> > &dr_beta,
    vector<vector<CCTK_REAL> > &dt_beta, vector<CCTK_REAL> &alpha,
    vector<CCTK_REAL> &dx_alpha, vector<CCTK_REAL> &dy_alpha,
    vector<CCTK_REAL> &dz_alpha, vector<CCTK_REAL> &dr_alpha,
    vector<CCTK_REAL> &dt_alpha, vector<CCTK_REAL> &th, vector<CCTK_REAL> &ph,
    vector<CCTK_REAL> &xhat, vector<CCTK_REAL> &yhat, vector<CCTK_REAL> &zhat,
    vector<CCTK_REAL> &xs, vector<CCTK_REAL> &ys, vector<CCTK_REAL> &zs,
    int ntheta, int nphi, int array_size, int rad_index);

} // namespace CCE_export

#endif
