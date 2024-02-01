#include "rainshaft_tendency_jac.hpp"
#include <stdexcept>

RainshaftTendencyJac::RainshaftTendencyJac(const double* t_jac,
                                     const double* q_jac,
                                     const double* nr_jac,
                                     const double* qr_jac)
  : t_tend_jac(t_jac), q_tend_jac(q_jac), nr_tend_jac(nr_jac), qr_tend_jac(qr_jac) {


  // // SPS: Should actually print out mismatched dimensions on failure.
  // if (t_vec.size() != q_vec.size()) {
  //   throw std::invalid_argument("t_vec and q_vec differ in size");
  // }
  // if (q_vec.size() != nr_vec.size()) {
  //   throw std::invalid_argument("q_vec and nr_vec differ in size");
  // }
  // if (nr_vec.size() != qr_vec.size()) {
  //   throw std::invalid_argument("nr_vec and qr_vec differ in size");
  // }
}
