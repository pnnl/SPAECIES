#include "rainshaft_grid.hpp"
#include "rainshaft_state.hpp"
#include <stdexcept>
#include <cmath>

std::vector<double> get_midpoints(const std::vector<double>& int_vec) {
  std::vector<double> mid_vec;
  int size = int_vec.size() - 1;
  for (int i = 0; i != size; ++i) {
    mid_vec.push_back(0.5 * (int_vec[i] + int_vec[i+1]));
  }
  return mid_vec;
}

RainshaftGrid::RainshaftGrid(const std::vector<double>& p_int_vec)
  : nlev(p_int_vec.size() - 1), p_int(p_int_vec), p_mid(get_midpoints(p_int))
{
  // Currently nothing to do.
}

// "ilev" copied from E3SM output; this is a set of interface pressure levels
// provided in hPa.
const std::size_t num_e3sm_levs = 72;
const double e3sm_ilev[num_e3sm_levs+1] = {
  0.100000003197544, 0.147650822913689, 0.218007648100008,
  0.321890076141867, 0.475273331103897, 0.701744962025609,
  1.03613217805539, 1.52985763845446, 2.25884732035832, 3.3352065502282,
  4.92445975982147, 7.01243897441468, 9.74236978346886, 13.2052046776003,
  17.4626717668335, 22.5300041895007, 28.3593888304225, 34.827113752148,
  41.905452436706, 49.4369434421476, 57.1821794266375, 64.8481839141014,
  72.1045965453237, 78.6060752491983, 85.2864750002711, 92.5346112988506,
  100.398735591023, 108.931198982486, 118.188803860441, 128.233178994491,
  139.131180869721, 150.95535019568, 163.784405484479, 177.703753558685,
  192.806050899617, 209.191828842591, 226.970165223411, 246.25940795157,
  267.187963858948, 289.895152760722, 314.532127483044, 341.262884695297,
  370.265374181464, 401.732673658358, 435.874252100169, 472.917375878188,
  512.019772166866, 551.259290548059, 589.990515658606, 627.297033184562,
  663.342902738581, 697.653186137899, 729.756090502968, 759.193568143422,
  785.532107743287, 808.373391694185, 827.364255260287, 842.826086512935,
  856.496411708255, 869.856439590426, 882.884857077784, 895.560615114759,
  907.863051465681, 919.771957375133, 931.267538919537, 942.330450725871,
  952.941866982882, 963.083665759796, 972.738495041755, 981.889632697215,
  990.52103362492, 996.992878983524, 1000
};

std::vector<double> RainshaftGrid::calc_dz(RainshaftConstants constants,
                                           std::vector<double> t_v) const {
  std::vector<double> dz;
  double rog = constants.rdry / constants.g;
  for (std::size_t i = 0; i != nlev; ++i) {
    double pdel = p_int[i+1] - p_int[i];
    dz.push_back(rog * t_v[i] * pdel / p_mid[i]);
  }
  return dz;
}

RainshaftGrid make_e3sm_like_grid(RainshaftConstants constants,
                                  double model_top, double srf_pres,
                                  double srf_temp, double lapse_rate) {
  // Convert e3sm_ilev to meters.
  std::vector<double> e3sm_ilev_meters;
  double length_scale = srf_temp / lapse_rate;
  double p_exp = constants.rdry * lapse_rate / constants.g;
  for (size_t i = 0; i != num_e3sm_levs+1; ++i) {
    // Ratio of interface level pressure to surface pressure.
    double ilev_rat = 100.*e3sm_ilev[i] / srf_pres;
    e3sm_ilev_meters.push_back(length_scale * (1 - std::pow(ilev_rat, p_exp)));
  }
  // Find the model top for the rainshaft model.
  std::size_t top_level_idx = 0;
  // SPS: Should throw an error when model_top is out of range.
  for (; top_level_idx != num_e3sm_levs; ++top_level_idx) {
    if (e3sm_ilev_meters[top_level_idx] < model_top) {
      --top_level_idx;
      break;
    }
  }
  // Calculate p_int so we can construct and return the grid object.
  std::vector<double> p_int;
  for (size_t i = top_level_idx; i != num_e3sm_levs+1; ++i) {
    p_int.push_back(100. * e3sm_ilev[i]);
  }
  return RainshaftGrid(p_int);
}
