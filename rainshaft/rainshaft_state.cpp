#include "rainshaft_state.hpp"
#include <stdexcept>

RainshaftState::RainshaftState(const std::vector<double>& t_vec,
                               const std::vector<double>& q_vec,
                               const std::vector<double>& nr_vec,
                               const std::vector<double>& qr_vec)
  : t(t_vec), q(q_vec), nr(nr_vec), qr(qr_vec) {
  // Check vector sizes.
  // SPS: Should actually print out mismatched dimensions on failure.
  if (t_vec.size() != q_vec.size()) {
    throw std::invalid_argument("t_vec and q_vec differ in size");
  }
  if (q_vec.size() != nr_vec.size()) {
    throw std::invalid_argument("q_vec and nr_vec differ in size");
  }
  if (nr_vec.size() != qr_vec.size()) {
    throw std::invalid_argument("nr_vec and qr_vec differ in size");
  }
}
