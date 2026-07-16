#include "rainshaft_constants.hpp"

namespace {
  VarDescList get_prognostic_variables(spaecies::Domain dom, spaecies::DimensionPtr lev_dim) {
    spaecies::VarDescPtr t_desc = dom.add_var_desc("T", spaecies::Float64Type, {lev_dim}, "K");
    spaecies::VarDescPtr q_desc = dom.add_var_desc("q", spaecies::Float64Type, {lev_dim}, "kg/kg");
    spaecies::VarDescPtr nr_desc = dom.add_var_desc("nr", spaecies::Float64Type, {lev_dim}, "1/kg");
    spaecies::VarDescPtr qr_desc = dom.add_var_desc("qr", spaecies::Float64Type, {lev_dim}, "kg/kg");
    return {t_desc, q_desc, nr_desc, qr_desc};
  }

  VarDescList get_diagnostic_variables(spaecies::Domain dom, spaecies::DimensionPtr lev_dim, const bool budget_diagnostics) {
    if (!budget_diagnostics) {
      return {};
    }

    spaecies::VarDescPtr evap_nr_desc = dom.add_var_desc("evap_nr", spaecies::Float64Type, {lev_dim}, "1/kg");
    spaecies::VarDescPtr evap_qr_desc = dom.add_var_desc("evap_qr", spaecies::Float64Type, {lev_dim}, "kg/kg");
    return {evap_nr_desc, evap_qr_desc};
  }
}

// CJV: can consider whether the remainder of this header should be in the namespace above

namespace defaults
{
  const double epsilon_qsat_fac = 1.e-10;
  const double epsilon_self_coll = 0.0;
  const double qsmall = 1.e-10;
  const double regularize_lambdar = true;
  const double regularize_qsat = true;
  const double use_lookup = false;
  const int use_zero_mur = 0;
};

RainshaftConstants create_RainshaftConstants(const double qsmall,
                                             const double epsilon_qsat_fac,
                                             const double epsilon_self_coll,
                                             const double mur)
{
  // SPS: Choose rho_top in a more principled way?
  return RainshaftConstants{3.14159265358979323846,
                            287.04, 1.00464e3, 461.50, 997., 2.501e6,
                            0.62197, qsmall, 9.80616, 1.e-5, 5.e-3, mur,
                            0.988919555598356, 1.e3, 1.e-4, epsilon_qsat_fac,
                            epsilon_self_coll};
}
