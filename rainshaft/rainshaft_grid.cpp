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


void generate_pcoord(double* e3sm_ilev, int num_e3sm_levs, const RainshaftConstants& constants) {
  if (constants.uniform_grid_flag)
  {
    double dp = (1000.0 - 1.0) / (num_e3sm_levs-1);

    for (std::size_t i = 0; i < num_e3sm_levs; ++i) {
      e3sm_ilev[i] = 1.0 + i*dp;
    }
  } else {
    // "ilev" copied from E3SM output; this is a set of interface pressure levels
    // provided in hPa.
    e3sm_ilev[0] = 0.100000003197544;
    e3sm_ilev[1] = 0.147650822913689;
    e3sm_ilev[2] = 0.218007648100008;
    e3sm_ilev[3] = 0.321890076141867;
    e3sm_ilev[4] = 0.475273331103897;
    e3sm_ilev[5] = 0.701744962025609;
    e3sm_ilev[6] = 1.03613217805539;
    e3sm_ilev[7] = 1.52985763845446;
    e3sm_ilev[8] = 2.25884732035832;
    e3sm_ilev[9] = 3.3352065502282;
    e3sm_ilev[10] = 4.92445975982147;
    e3sm_ilev[11] = 7.01243897441468;
    e3sm_ilev[12] = 9.74236978346886;
    e3sm_ilev[13] = 13.2052046776003;
    e3sm_ilev[14] = 17.4626717668335;
    e3sm_ilev[15] = 22.5300041895007;
    e3sm_ilev[16] = 28.3593888304225;
    e3sm_ilev[17] = 34.827113752148;
    e3sm_ilev[18] = 41.905452436706;
    e3sm_ilev[19] = 49.4369434421476;
    e3sm_ilev[20] = 57.1821794266375;
    e3sm_ilev[21] = 64.8481839141014;
    e3sm_ilev[22] = 72.1045965453237;
    e3sm_ilev[23] = 78.6060752491983;
    e3sm_ilev[24] = 85.2864750002711;
    e3sm_ilev[25] = 92.5346112988506;
    e3sm_ilev[26] = 100.398735591023;
    e3sm_ilev[27] = 108.931198982486;
    e3sm_ilev[28] = 118.188803860441;
    e3sm_ilev[29] = 128.233178994491;
    e3sm_ilev[30] = 139.131180869721;
    e3sm_ilev[31] = 150.95535019568;
    e3sm_ilev[32] = 163.784405484479;
    e3sm_ilev[33] = 177.703753558685;
    e3sm_ilev[34] = 192.806050899617;
    e3sm_ilev[35] = 209.191828842591;
    e3sm_ilev[36] = 226.970165223411;
    e3sm_ilev[37] = 246.25940795157;
    e3sm_ilev[38] = 267.187963858948;
    e3sm_ilev[39] = 289.895152760722;
    e3sm_ilev[40] = 314.532127483044;
    e3sm_ilev[41] = 341.262884695297;
    e3sm_ilev[42] = 370.265374181464;
    e3sm_ilev[43] = 401.732673658358;
    e3sm_ilev[44] = 435.874252100169;
    e3sm_ilev[45] = 472.917375878188;
    e3sm_ilev[46] = 512.019772166866;
    e3sm_ilev[47] = 551.259290548059;
    e3sm_ilev[48] = 589.990515658606;
    e3sm_ilev[49] = 627.297033184562;
    e3sm_ilev[50] = 663.342902738581;
    e3sm_ilev[51] = 697.653186137899;
    e3sm_ilev[52] = 729.756090502968;
    e3sm_ilev[53] = 759.193568143422;
    e3sm_ilev[54] = 785.532107743287;
    e3sm_ilev[55] = 808.373391694185;
    e3sm_ilev[56] = 827.364255260287;
    e3sm_ilev[57] = 842.826086512935;
    e3sm_ilev[58] = 856.496411708255;
    e3sm_ilev[59] = 869.856439590426;
    e3sm_ilev[60] = 882.884857077784;
    e3sm_ilev[61] = 895.560615114759;
    e3sm_ilev[62] = 907.863051465681;
    e3sm_ilev[63] = 919.771957375133;
    e3sm_ilev[64] = 931.267538919537;
    e3sm_ilev[65] = 942.330450725871;
    e3sm_ilev[66] = 952.941866982882;
    e3sm_ilev[67] = 963.083665759796;
    e3sm_ilev[68] = 972.738495041755;
    e3sm_ilev[69] = 981.889632697215;
    e3sm_ilev[70] = 990.52103362492;
    e3sm_ilev[71] = 996.992878983524;
    e3sm_ilev[72] = 1000.;
  }

}



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
  int num_e3sm_levs;

  if (constants.uniform_grid_flag)
  {
    num_e3sm_levs = constants.num_uniform_grid_points;
  } else {
    num_e3sm_levs = 72;
  }

  double e3sm_ilev[num_e3sm_levs+1];
  generate_pcoord(e3sm_ilev, num_e3sm_levs+1, constants);

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
  for (size_t i = top_level_idx; i != num_e3sm_levs; ++i) {
    p_int.push_back(100. * e3sm_ilev[i]);
  }
  return RainshaftGrid(p_int);
}
