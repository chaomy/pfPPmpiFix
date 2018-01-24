/*
 * @Author: chaomingyang
 * @Date:   2017-11-02 22:06:59
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-13 21:56:39
 */

#include "pfOptimizer.h"

using std::cout;
using std::endl;

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
  nlopt::algorithm nalg;
  const string atag = pfhm->gsparams()["alg"];
  cout<<"atag"<<atag<<endl;
  if (atag == "LD_MMA") nalg = nlopt::LD_MMA; 
  else if (atag == "LD_LBFGS") nalg = nlopt::LD_LBFGS;
  else if (atag == "GN_ESCH") nalg = nlopt::GN_ESCH; 
  else if (atag == "GN_ISRES") nalg = nlopt::GN_ISRES;
  else if (atag == "LN_SBPLX") nalg = nlopt::LN_SBPLX;
  else nalg = nlopt::GN_CRS2_LM;

  opt = new nlopt::opt(nalg, nvars);
  // bounds 
  opt->set_lower_bounds(pfhm->lob);
  opt->set_upper_bounds(pfhm->hib);

  // targets   
  opt->set_min_objective(objWrapper, this);
  // constrain for rho
  // opt->add_inequality_constraint(cnt1Wrapper, this);
  // opt->add_inequality_constraint(cnt2Wrapper, this);
  cout<<"maxstep "<<pfhm->giparams()["maxstep"];
  cout<<"xtol "<<pfhm->gdparams()["xtol"];

  opt->set_maxeval(pfhm->giparams()["maxstep"]); 
  opt->set_xtol_rel(pfhm->gdparams()["xtol"]);
}

void pfOptimizer::runOpt() {   
  double minf;
  nlopt::result result = opt->optimize(pfhm->ini, minf);
  if (result) printf("found min\n");
  delete opt;
}

double pfOptimizer::objectiveFuc(const vector<double> &x, vector<double> &grad) {
    if (gtag) return pfhm->errFunctGrad(x, grad);
    else return pfhm->errFunct(x);
}

double pfOptimizer::addConstraint1(const vector<double> &x, vector<double> &grad) {
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