/*
 * @Author: chaomingyang
 * @Date:   2017-11-02 22:06:59
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-02 23:44:49
 */

#include "pfOptimizer.h"
#include "pfHome.h"

double objWrapper(const vector<double> &x, vector<double> &g, void *data) {
  pfOptimizer *obj = static_cast<pfOptimizer *>(data);
  return obj->objectiveFuc(x, g);
}

double cnt1Wrapper(const vector<double> &x, vector<double> &g, void *data) {
  pfOptimizer *obj = static_cast<pfOptimizer *>(data);
  return obj->addConstraint1(x, g);
}

double cnt2Wrapper(const vector<double> &x, vector<double> &g, void *data) {
  pfOptimizer *obj = static_cast<pfOptimizer *>(data);
  return obj->addConstraint2(x, g);
}

void pfOptimizer::prepOpt() {
  opt = new nlopt::opt(nlopt::LD_MMA, 2);
  vector<double> lb(2);
  lb[0] = -1e10;
  lb[1] = 0.0;
  opt->set_lower_bounds(lb);
  opt->set_min_objective(objWrapper, this);
  opt->add_inequality_constraint(cnt1Wrapper, this);
  opt->add_inequality_constraint(cnt2Wrapper, this);
  opt->set_xtol_rel(1e-4);
}

void pfOptimizer::runOpt() {
  vector<double> x(2);
  x[0] = 1.2;
  x[1] = 5.6;
  double minf;
  nlopt::result result = opt->optimize(x, minf);
  printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
  delete opt;
}

double pfOptimizer::objectiveFuc(const vector<double> &x,
                                 vector<double> &grad) {
  if (!grad.empty()) {
    grad[0] = 0.0;
    grad[1] = 0.5 / sqrt(x[1]);
  }
  return sqrt(x[1]);
}

double pfOptimizer::addConstraint1(const vector<double> &x,
                                   vector<double> &grad) {
  double a = 2., b = 0.0;
  if (!grad.empty()) {
    grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
    grad[1] = -1.0;
  }
  return ((a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1]);
}

double pfOptimizer::addConstraint2(const vector<double> &x,
                                   vector<double> &grad) {
  double a = -1., b = 1.0;
  if (!grad.empty()) {
    grad[0] = 3 * a * (a * x[0] + b) * (a * x[0] + b);
    grad[1] = -1.0;
  }
  return ((a * x[0] + b) * (a * x[0] + b) * (a * x[0] + b) - x[1]);
}