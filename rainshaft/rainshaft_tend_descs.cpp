#include "rainshaft_tend_descs.hpp"

std::vector<spaecies::VarDescPtr> tend_descs_from_state_descs(spaecies::Domain dom,
                                                              std::vector<spaecies::VarDescPtr> state_descs) {
  std::vector<spaecies::VarDescPtr> tend_descs;
  for (spaecies::VarDescPtr var_desc : state_descs) {
    auto name = var_desc->name + "_tend";
    auto type = var_desc->type;
    auto dimensions = var_desc->dimensions;
    auto units = var_desc->units + "/s";
    tend_descs.push_back(dom.add_var_desc(name, type, dimensions, units));
  }
  return tend_descs;
}
