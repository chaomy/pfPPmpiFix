#ifndef pfOpt_H_
#define pfOpt_H_

#include "pfHome.h"

class pfOptimizer {
 private:
  int nvars;
  bool gtag;
  pfHome* pfhm;
  nlopt::opt* opt;

 public:
  pfOptimizer();
  pfOptimizer(pfHome* pf) : pfhm(pf) {
    nvars = pfhm->ini.size();
    gtag = false;
  };
  void prepOpt();
  void runOpt();
  double objectiveFuc(const vector<double>& x, vector<double>& g);
  double addConstraint1(const vector<double>& x, vector<double>& g);
  double addConstraint2(const vector<double>& x, vector<double>& g);
};

#endif  // pfOptimizer