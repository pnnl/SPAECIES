#ifndef RAINSHAFT_TEND_DESCS_HPP
#define RAINSHAFT_TEND_DESCS_HPP
#include <spaecies.hpp>

std::vector<spaecies::VarDescPtr> tend_descs_from_state_descs(spaecies::Domain dom,
                                                              std::vector<spaecies::VarDescPtr> state_descs);

#endif // RAINSHAFT_TEND_DESCS_HPP
