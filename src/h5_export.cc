#include "h5_export.hh"

#if defined __cpp_lib_filesystem && __cpp_lib_filesystem < 201703L
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#else
#include <filesystem>
namespace fs = std::filesystem;
#endif

namespace CCE_export {

using std::string;
using std::ostringstream;
using std::map;
using std::ios;
using std::setprecision;

#define HDF5_ERROR(fn_call)                                                    \
  do {                                                                         \
    hid_t _error_code = fn_call;                                               \
                                                                               \
    if (_error_code < 0) {                                                     \
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,                      \
                 "HDF5 call '%s' returned error code %d", #fn_call,            \
                 static_cast<int>(_error_code));                               \
    }                                                                          \
  } while (0)

static bool dataset_exists(hid_t file, const string &dataset_name) {
  // To test whether a dataset exists, the recommended way in API 1.6
  // is to use H5Gget_objinfo, but this prints an error to stderr if
  // the dataset does not exist.  We explicitly avoid this by wrapping
  // the call in H5E_BEGIN_TRY/H5E_END_TRY statements.  In 1.8,
  // H5Gget_objinfo is deprecated, and H5Lexists does the job.  See
  // http://www.mail-archive.com/hdf-forum@hdfgroup.org/msg00125.html

  bool exists;
  H5E_BEGIN_TRY {
    exists = H5Gget_objinfo(file, dataset_name.c_str(), 1, NULL) >= 0;
  }
  H5E_END_TRY;
  return exists;
}

void Create_Dataset(CCTK_ARGUMENTS, hid_t file, string datasetname,
                    CCTK_REAL *data, int lmax) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int mode_count = l_m_to_index(lmax, lmax) + 1;

  hid_t dataset = -1;

  if (dataset_exists(file, datasetname)) {
    HDF5_ERROR(dataset = H5Dopen(file, datasetname.c_str()));
  } else {
    hsize_t dims[2] = {0, static_cast<hsize_t>(2 * mode_count + 1)};
    hsize_t maxdims[2] = {H5S_UNLIMITED,
                          static_cast<hsize_t>(2 * mode_count + 1)};
    hid_t dataspace = H5Screate_simple(2, dims, maxdims);

    hid_t cparms = -1;
    hsize_t chunk_dims[2] = {static_cast<hsize_t>(hdf5_chunk_size),
                             static_cast<hsize_t>(2 * mode_count + 1)};
    cparms = H5Pcreate(H5P_DATASET_CREATE);
    HDF5_ERROR(H5Pset_chunk(cparms, 2, chunk_dims));

    dataset = H5Dcreate(file, datasetname.c_str(), H5T_NATIVE_DOUBLE, dataspace,
                        cparms);
    H5Pclose(cparms);

    // Create the legend

    // Create variable-length string datatype
    hid_t str_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_type, H5T_VARIABLE);

    // Create dataspace for array of strings
    hsize_t legend_dims[] = {static_cast<hsize_t>(2 * mode_count + 1)};
    hid_t space = H5Screate_simple(1, legend_dims, NULL);

    // Create attribute
    hid_t attr = H5Acreate(dataset, "Legend", str_type, space, H5P_DEFAULT);

    // Write array of strings
    vector<const char *> legend(2 * mode_count + 1);
    legend.at(0) = "time";
    for (int l = 0; l < lmax + 1; l++) {
      for (int m = -l; m < l + 1; m++) {
        int mode_index = l_m_to_index(l, m);
        ostringstream re_label;
        re_label << "Re(" << l << "," << m << ")";
        legend.at(2 * mode_index + 1) = strdup(re_label.str().c_str());
        ostringstream im_label;
        im_label << "Im(" << l << "," << m << ")";
        legend.at(2 * mode_index + 2) = strdup(im_label.str().c_str());
      }
    }
    H5Awrite(attr, str_type, legend.data());

    // Clean up by closing the attribute, dataspace, and datatype
    H5Aclose(attr);
    H5Sclose(space);
    H5Tclose(str_type);

    for (size_t i = 0; i < legend.size(); i++) {
      if (legend[i] != nullptr && strcmp(legend[i], "time") != 0) {
        free((void*)legend[i]);
      }
    }
  }

  hid_t filespace = H5Dget_space(dataset);

  hsize_t filedims[2];
  hsize_t maxdims[2];
  HDF5_ERROR(H5Sget_simple_extent_dims(filespace, filedims, maxdims));

  filedims[0] += 1;
  hsize_t size[2] = {filedims[0], filedims[1]};
  HDF5_ERROR(H5Dextend(dataset, size));
  HDF5_ERROR(H5Sclose(filespace));

  /* Select a hyperslab  */
  hsize_t offset[2] = {filedims[0] - 1, 0};
  // hsize_t dims2[2] = {1, 3};
  hsize_t dims2[2] = {1, (hsize_t)(2 * mode_count + 1)};
  filespace = H5Dget_space(dataset);
  HDF5_ERROR(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims2,
                                 NULL));

  hid_t memdataspace = H5Screate_simple(2, dims2, NULL);

  /* Write the data to the hyperslab  */
  HDF5_ERROR(H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memdataspace, filespace,
                      H5P_DEFAULT, data));

  HDF5_ERROR(H5Dclose(dataset));
  HDF5_ERROR(H5Sclose(filespace));
  HDF5_ERROR(H5Sclose(memdataspace));
}

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
    vector<CCTK_REAL> &im_dt_alpha, float rad, int lmax) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int mode_count = l_m_to_index(lmax, lmax) + 1;

  const char *my_out_dir = strcmp(out_dir, "") ? out_dir : io_out_dir;
  if (CCTK_CreateDirectory(0755, my_out_dir) < 0)
    CCTK_VWarn(CCTK_WARN_ABORT, __LINE__, __FILE__, CCTK_THORNSTRING,
               "CCE_Export output directory %s could not be created",
               my_out_dir);

  static map<string, bool> checked;

  ostringstream basename;
  basename << base_file_name << "R" << setiosflags(ios::fixed) << setprecision(2)
           << rad << "." << extension;
  string output_name = (fs::path(my_out_dir) / basename.str()).string();

  hid_t file;

  if (!fs::exists(output_name) ||
      (!checked[output_name] && IO_TruncateOutputFiles(cctkGH))) {
    HDF5_ERROR(file =
	       H5Fcreate(output_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT));
  } else {
    HDF5_ERROR(file = H5Fopen(output_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
  }

  checked[output_name] = true;

  // store metric data
  for (int i = 0; i < 3; i++) {
    for (int j = i; j < 3; j++) {
      string components;
      if (i == 0 && j == 0) {
        components = "xx";
      }
      if (i == 0 && j == 1) {
        components = "xy";
      }
      if (i == 0 && j == 2) {
        components = "xz";
      }
      if (i == 1 && j == 1) {
        components = "yy";
      }
      if (i == 1 && j == 2) {
        components = "yz";
      }
      if (i == 2 && j == 2) {
        components = "zz";
      }

      string metric_datasetname = "g" + components + ".dat";
      string metric_dt_datasetname = "Dtg" + components + ".dat";
      string metric_dr_datasetname = "Drg" + components + ".dat";

      vector<CCTK_REAL> metric_data(2 * mode_count + 1);
      vector<CCTK_REAL> metric_dt_data(2 * mode_count + 1);
      vector<CCTK_REAL> metric_dr_data(2 * mode_count + 1);
      metric_data.at(0) = cctk_time;
      metric_dt_data.at(0) = cctk_time;
      metric_dr_data.at(0) = cctk_time;
      for (int l = 0; l <= lmax; l++) {
        for (int m = -l; m < l + 1; m++) {
          int mode_index = l_m_to_index(l, m);
          metric_data.at(2 * mode_index + 1) = re_g.at(i).at(j).at(mode_index);
          metric_data.at(2 * mode_index + 2) = im_g.at(i).at(j).at(mode_index);
          metric_dt_data.at(2 * mode_index + 1) =
              re_dt_g.at(i).at(j).at(mode_index);
          metric_dt_data.at(2 * mode_index + 2) =
              im_dt_g.at(i).at(j).at(mode_index);
          metric_dr_data.at(2 * mode_index + 1) =
              re_dr_g.at(i).at(j).at(mode_index);
          metric_dr_data.at(2 * mode_index + 2) =
              im_dr_g.at(i).at(j).at(mode_index);
        }
      }

      Create_Dataset(CCTK_PASS_CTOC, file, metric_datasetname,
                     metric_data.data(), lmax);
      Create_Dataset(CCTK_PASS_CTOC, file, metric_dt_datasetname,
                     metric_dt_data.data(), lmax);
      Create_Dataset(CCTK_PASS_CTOC, file, metric_dr_datasetname,
                     metric_dr_data.data(), lmax);
    }
  }

  // store shift data
  for (int i = 0; i < 3; i++) {
    string component;
    if (i == 0) {
      component = "x";
    }
    if (i == 1) {
      component = "y";
    }
    if (i == 2) {
      component = "z";
    }

    string shift_datasetname = "Shift" + component + ".dat";
    string shift_dt_datasetname = "DtShift" + component + ".dat";
    string shift_dr_datasetname = "DrShift" + component + ".dat";

    vector<CCTK_REAL> shift_data(2 * mode_count + 1);
    vector<CCTK_REAL> shift_dt_data(2 * mode_count + 1);
    vector<CCTK_REAL> shift_dr_data(2 * mode_count + 1);
    shift_data.at(0) = cctk_time;
    shift_dt_data.at(0) = cctk_time;
    shift_dr_data.at(0) = cctk_time;

    for (int l = 0; l <= lmax; l++) {
      for (int m = -l; m < l + 1; m++) {
        int mode_index = l_m_to_index(l, m);
        shift_data.at(2 * mode_index + 1) = re_beta.at(i).at(mode_index);
        shift_data.at(2 * mode_index + 2) = im_beta.at(i).at(mode_index);
        shift_dt_data.at(2 * mode_index + 1) = re_dt_beta.at(i).at(mode_index);
        shift_dt_data.at(2 * mode_index + 2) = im_dt_beta.at(i).at(mode_index);
        shift_dr_data.at(2 * mode_index + 1) = re_dr_beta.at(i).at(mode_index);
        shift_dr_data.at(2 * mode_index + 2) = im_dr_beta.at(i).at(mode_index);
      }
    }

    Create_Dataset(CCTK_PASS_CTOC, file, shift_datasetname, shift_data.data(),
                   lmax);
    Create_Dataset(CCTK_PASS_CTOC, file, shift_dt_datasetname,
                   shift_dt_data.data(), lmax);
    Create_Dataset(CCTK_PASS_CTOC, file, shift_dr_datasetname,
                   shift_dr_data.data(), lmax);
  }

  // store lapse data
  string lapse_datasetname = "Lapse.dat";
  string lapse_dt_datasetname = "DtLapse.dat";
  string lapse_dr_datasetname = "DrLapse.dat";

  vector<CCTK_REAL> lapse_data(2 * mode_count + 1);
  vector<CCTK_REAL> lapse_dt_data(2 * mode_count + 1);
  vector<CCTK_REAL> lapse_dr_data(2 * mode_count + 1);
  lapse_data[0] = cctk_time;
  lapse_dt_data[0] = cctk_time;
  lapse_dr_data[0] = cctk_time;
  for (int l = 0; l <= lmax; l++) {
    for (int m = -l; m < l + 1; m++) {
      int mode_index = l_m_to_index(l, m);
      lapse_data.at(2 * mode_index + 1) = re_alpha.at(mode_index);
      lapse_data.at(2 * mode_index + 2) = im_alpha.at(mode_index);
      lapse_dt_data.at(2 * mode_index + 1) = re_dt_alpha.at(mode_index);
      lapse_dt_data.at(2 * mode_index + 2) = im_dt_alpha.at(mode_index);
      lapse_dr_data.at(2 * mode_index + 1) = re_dr_alpha.at(mode_index);
      lapse_dr_data.at(2 * mode_index + 2) = im_dr_alpha.at(mode_index);
    }
  }

  Create_Dataset(CCTK_PASS_CTOC, file, lapse_datasetname, lapse_data.data(),
                 lmax);
  Create_Dataset(CCTK_PASS_CTOC, file, lapse_dt_datasetname,
                 lapse_dt_data.data(), lmax);
  Create_Dataset(CCTK_PASS_CTOC, file, lapse_dr_datasetname,
                 lapse_dr_data.data(), lmax);

  HDF5_ERROR(H5Fclose(file));
}

} // namespace CCE_export
