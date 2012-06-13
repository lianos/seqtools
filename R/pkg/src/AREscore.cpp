#include "AREscore.h"

SEXP find_au_start_end(SEXP au_, SEXP pstart_, SEXP pend_) {
BEGIN_RCPP
  Rcpp::NumericVector au(au_);
  double pstart = Rcpp::as<double>(pstart_);
  double pend = Rcpp::as<double>(pend_);

  std::vector<int> starts = std::vector<int>();
  std::vector<int> ends = std::vector<int>();

  bool in_block = false;
  double val = 0.0;

  for (int i = 0; i < au.size(); i++) {
    val = au[i];
    if (!in_block && val >= pstart) {
      in_block = true;
      starts.push_back(i + 1);
    } else if (in_block && val < pend) {
      in_block = false;
      ends.push_back(i + 1);
    }
  }

  if (starts.size() == 0) {
    starts.push_back(-1);
    ends.push_back(-1);
  }
  if (ends.size() == starts.size() - 1) {
    ends.push_back(au.size());
  }

  if (starts.size() != ends.size()) {
    throw std::range_error("Size of starts and ends do not match");
  }

  return Rcpp::DataFrame::create(Rcpp::_["start"]=starts, Rcpp::_["end"]=ends);
END_RCPP
}
