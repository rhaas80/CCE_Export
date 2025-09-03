#ifndef H5_EXPORT_HH
#define H5_EXPORT_HH

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "utils.hh"

#include <string.h>
#include <vector>
#include <map>
#include <iomanip>

#define H5_USE_16_API
#include <hdf5.h>

namespace CCE_export {

using std::vector;
using std::string;

void Create_Dataset(string datasetname, CCTK_REAL *data, int mode_count);

void Output_Decomposed_Metric_Data(
    CCTK_ARGUMENTS, vector<vector<vector<CCTK_REAL> > > &re_g,
    vector<vector<vector<CCTK_REAL> > > &im_g,
    vector<vector<vector<CCTK_REAL> > > &re_dr_g,
    vector<vector<vector<CCTK_REAL> > > &im_dr_g,
    vector<vector<vector<CCTK_REAL> > > &re_dt_g,
    vector<vector<vector<CCTK_REAL> > > &im_dt_g,
    vector<vector<CCTK_REAL> > &re_beta, vector<vector<CCTK_REAL> > &im_beta,
    vector<vector<CCTK_REAL> > &re_dr_beta,
    vector<vector<CCTK_REAL> > &im_dr_beta,
    vector<vector<CCTK_REAL> > &re_dt_beta,
    vector<vector<CCTK_REAL> > &im_dt_beta, vector<CCTK_REAL> &re_alpha,
    vector<CCTK_REAL> &im_alpha, vector<CCTK_REAL> &re_dr_alpha,
    vector<CCTK_REAL> &im_dr_alpha, vector<CCTK_REAL> &re_dt_alpha,
    vector<CCTK_REAL> &im_dt_alpha, float radius, int lmax);

} // namespace CCE_export

#endif
