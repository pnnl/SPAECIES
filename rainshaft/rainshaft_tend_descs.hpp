#ifndef RAINSHAFT_TEND_DESCS_HPP
#define RAINSHAFT_TEND_DESCS_HPP
#include <spaecies.hpp>

#include "rainshaft_types.hpp"

VarDescList tend_descs_from_state_descs(spaecies::Domain dom,
                                        VarDescList state_descs);

#endif // RAINSHAFT_TEND_DESCS_HPP
