#ifndef CCE_EXPORT_HH
#define CCE_EXPORT_HH

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"
#include "util_Table.h"

#include <string>
#include <vector>

namespace CCE_export {

using std::vector;
using std::string;

extern "C" void CCE_Export(CCTK_ARGUMENTS);

} // namespace CCE_export

#endif // CCE_EXPORT_HH
