#pragma once

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace karst {

inline std::string make_key(const std::vector<int>& x) {
  std::string out;
  out.reserve(x.size() * 3);
  for (size_t i = 0; i < x.size(); ++i) {
    if (i) out.push_back('.');
    out += std::to_string(x[i]);
  }
  return out;
}

// Copy numbers are unsigned decimal integers separated by single dots: no
// sign, no space, no empty field (also none at the end, which getline() would
// drop). Leading zeros are read as the number they write; callers that match
// labels as strings first rewrite them with canonical_karyotype_labels().
inline std::vector<int> parse_label(const std::string& s) {
  std::vector<int> out;
  if (s.empty()) throw std::invalid_argument("Malformed karyotype label: empty string.");
  if (s.back() == '.') throw std::invalid_argument("Malformed karyotype label \'" + s + "\': empty copy-number field.");
  std::string token;
  std::stringstream ss(s);
  while (std::getline(ss, token, '.')) {
    if (token.empty()) throw std::invalid_argument("Malformed karyotype label \'" + s + "\': empty copy-number field.");
    for (char c : token) {
      if (c < '0' || c > '9') {
        throw std::invalid_argument("Malformed karyotype label \'" + s + "\': copy numbers must be integers.");
      }
    }
    char* end = nullptr;
    errno = 0;
    long v = std::strtol(token.c_str(), &end, 10);
    if (end == token.c_str() || *end != '\0' || errno == ERANGE) {
      throw std::invalid_argument("Malformed karyotype label \'" + s + "\': copy numbers must be integers.");
    }
    if (v < 0 || v > std::numeric_limits<int>::max()) {
      throw std::invalid_argument("Malformed karyotype label \'" + s + "\': copy number outside integer range.");
    }
    out.push_back(static_cast<int>(v));
  }
  if (out.empty()) throw std::invalid_argument("Malformed karyotype label: empty string.");
  return out;
}

inline double choose_double(int n, int k) {
  if (k < 0 || k > n) return 0.0;
  if (k == 0 || k == n) return 1.0;
  k = std::min(k, n - k);
  double out = 1.0;
  for (int i = 1; i <= k; ++i) {
    out *= static_cast<double>(n - k + i) / i;
  }
  return out;
}

inline double chromosome_transition_probability(int parent_cn, int child_cn, double beta) {
  if (parent_cn < 0 || child_cn < 0) return 0.0;
  int diff = std::abs(parent_cn - child_cn);
  if (diff > parent_cn) return 0.0;
  double total = 0.0;
  for (int m = diff; m <= parent_cn; m += 2) {
    int lost_to_child = (m + parent_cn - child_cn) / 2;
    if (lost_to_child < 0 || lost_to_child > m) continue;
    total += choose_double(parent_cn, m) *
      std::pow(beta, m) *
      std::pow(1.0 - beta, parent_cn - m) *
      std::pow(0.5, m) *
      choose_double(m, lost_to_child);
  }
  return total;
}

} // namespace karst
