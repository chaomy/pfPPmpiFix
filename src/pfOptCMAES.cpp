/*
 * @Xuthor: chaomy
 * @Date:   2018-01-10 20:08:18
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-24 11:03:50
 *
 * Modified from mlpack
 * Implementation of the Covariance Matrix Adaptation Evolution Strategy as
 * proposed by N. Hansen et al. in "Completely Derandomized Self-Adaptation in
 * Evolution Strategies".
 *
 */

#include "pfHome.h"

using arma::accu;
using arma::linspace;
using arma::mat;
using arma::randu;
using arma::vec;
using std::cout;
using std::endl;
using std::to_string;

double pfHome::testFunc(arma::mat& vc) {
  return pow((vc[0] - 2.0), 2) + pow((vc[1] - 3.0), 2) + pow((vc[2] - 0.01), 2);
}

void pfHome::cntcmaes() {
  arma::mat iterate = encodev(ini);  // [0, 10]
  (this->*calobj[sparams["ptype"]])(iterate, 1);
  if (cmm.rank() == PFROOT) {
    double cr = cmaes(iterate);
    writePot();
    sparams["lmpfile"] = "dummy." + sparams["ptype"];
    (this->*write[sparams["ptype"]])();
    (this->*calobj[sparams["ptype"]])(iterate, EXT);
  }
}

void pfHome::loopcmaes() {
  arma::mat iterate(nvars, 1, arma::fill::randu);
  double cr = 1e30, op = 1e30;
  for (int i = 0; i < 1; i++) {
    for (int k = 0; k < nvars; k++)  // random
      iterate[k] = lob[k] + randUniform() * (hib[k] - lob[k]);

    (this->*calobj[sparams["ptype"]])(iterate, 1);
    if (cmm.rank() == PFROOT) {
      op = (cr = cmaes(iterate)) < op ? cr : op;
      writePot("dummy.tmp." + to_string(i) + "_" + to_string(cr));
      sparams["lmpfile"] = "dummy." + sparams["ptype"] + "." + to_string(i);
      (this->*write[sparams["ptype"]])();
      (this->*calobj[sparams["ptype"]])(iterate, EXT);
    }
  }
}

double pfHome::cmaes(arma::mat& iterate) {
  int maxIt = 50000;
  double tolerance = 1e-18;

  // Population size.
  int lambda = (4 + std::round(3 * std::log(iterate.n_elem))) * 10;

  // Parent weights.
  const size_t mu = std::round(lambda / 2);

  cout << "nrows: " << iterate.n_rows << " lambda: " << lambda << " mu: " << mu
       << endl;

  arma::vec w = std::log(mu + 0.5) -
                arma::log(arma::linspace<arma::vec>(0, mu - 1, mu) + 1.0);
  w /= arma::sum(w);

  // Number of effective solutions.
  const double muEffective = 1 / arma::accu(arma::pow(w, 2));

  // Step size control parameters.
  arma::vec sigma(3);

  double upperBound = 10., lowerBound = -10.;
  sigma(0) = 0.1 * (upperBound - lowerBound);

  const double cs = (muEffective + 2) / (iterate.n_elem + muEffective + 5);
  const double ds =
      1 + cs +
      2 * std::max(std::sqrt((muEffective - 1) / (iterate.n_elem + 1)) - 1,
                   0.0);
  const double enn =
      std::sqrt(iterate.n_elem) * (1.0 - 1.0 / (4.0 * iterate.n_elem) +
                                   1.0 / (21 * std::pow(iterate.n_elem, 2)));

  // Covariance update parameters.
  // Cumulation for distribution.
  const double cc = (4 + muEffective / iterate.n_elem) /
                    (4 + iterate.n_elem + 2 * muEffective / iterate.n_elem);
  const double h = (1.4 + 2.0 / (iterate.n_elem + 1.0)) * enn;

  const double c1 = 2 / (std::pow(iterate.n_elem + 1.3, 2) + muEffective);
  const double alphaMu = 2;
  const double cmu = std::min(
      1 - c1,
      alphaMu * (muEffective - 2 + 1 / muEffective) /
          (std::pow(iterate.n_elem + 2, 2) + alphaMu * muEffective / 2));

  arma::cube mps(iterate.n_rows, iterate.n_cols, 3);  // meam
  // mps.slice(0) =
  //     lowerBound +
  //     arma::randu(iterate.n_rows, iterate.n_cols) * (upperBound -
  //     lowerBound);
  mps.slice(0) = iterate;

  arma::mat step = arma::zeros(iterate.n_rows, iterate.n_cols);

  // Calculate the first objective function.
  double currentobj =
      (this->*calobj[sparams["ptype"]])(decodev(mps.slice(0)), 1);
  double overallobj = currentobj;
  double lastobj = 1e30;

  // Population parameters.
  arma::cube pStep(iterate.n_rows, iterate.n_cols, lambda);
  arma::cube pps(iterate.n_rows, iterate.n_cols, lambda);
  arma::vec pobj(lambda);
  arma::cube ps = arma::zeros(iterate.n_rows, iterate.n_cols, 2);
  arma::cube pc = ps;
  arma::cube C(iterate.n_elem, iterate.n_elem, 2);
  C.slice(0).eye();

  // Covariance matrix parameters.
  arma::vec eigval;
  arma::mat eigvec;
  arma::vec eigvalZero = arma::zeros(iterate.n_elem);

  // The current visitation order (sorted by population objectives).
  arma::uvec idx = arma::linspace<arma::uvec>(0, lambda - 1, lambda);

  for (size_t i = 1; i < maxIt; ++i) {
    const size_t idx0 = (i - 1) % 2;
    const size_t idx1 = i % 2;

    const arma::mat covLower = arma::chol(C.slice(idx0), "lower");

    for (size_t j = 0; j < lambda; ++j) {
      if (iterate.n_rows > iterate.n_cols) {
        pStep.slice(idx(j)) =
            covLower * arma::randn(iterate.n_rows, iterate.n_cols);
      } else {
        pStep.slice(idx(j)) =
            arma::randn(iterate.n_rows, iterate.n_cols) * covLower;
      }
      pps.slice(idx(j)) = mps.slice(idx0) + sigma(idx0) * pStep.slice(idx(j));
      pobj(idx(j)) =
          (this->*calobj[sparams["ptype"]])(decodev(pps.slice(idx(j))), 1);
    }

    // Sort population.
    idx = sort_index(pobj);

    step = w(0) * pStep.slice(idx(0));
    for (size_t j = 1; j < mu; ++j) step += w(j) * pStep.slice(idx(j));

    mps.slice(idx1) = mps.slice(idx0) + sigma(idx0) * step;

    currentobj = (this->*calobj[sparams["ptype"]])(decodev(mps.slice(idx1)), 1);

    // Update best parameters.
    if (currentobj < overallobj) {
      overallobj = currentobj;
      iterate = mps.slice(idx1);
      // output
      writePot("dummy.tmp");
      (this->*write[sparams["ptype"]])();
    }

    // Update Step Size.
    if (iterate.n_rows > iterate.n_cols) {
      ps.slice(idx1) =
          (1 - cs) * ps.slice(idx0) +
          std::sqrt(cs * (2 - cs) * muEffective) * covLower.t() * step;
    } else {
      ps.slice(idx1) =
          (1 - cs) * ps.slice(idx0) +
          std::sqrt(cs * (2 - cs) * muEffective) * step * covLower.t();
    }

    const double psNorm = arma::norm(ps.slice(idx1));
    sigma(idx1) =
        sigma(idx0) * std::pow(std::exp(cs / ds * psNorm / enn - 1), 0.3);

    // Update covariance matrix.
    if ((psNorm / sqrt(1 - std::pow(1 - cs, 2 * i))) < h) {
      pc.slice(idx1) = (1 - cc) * pc.slice(idx0) +
                       std::sqrt(cc * (2 - cc) * muEffective) * step;

      if (iterate.n_rows > iterate.n_cols) {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1) * pc.slice(idx1).t());
      } else {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1).t() * pc.slice(idx1));
      }
    } else {
      pc.slice(idx1) = (1 - cc) * pc.slice(idx0);

      if (iterate.n_rows > iterate.n_cols) {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1) * pc.slice(idx1).t() +
                              (cc * (2 - cc)) * C.slice(idx0));
      } else {
        C.slice(idx1) = (1 - c1 - cmu) * C.slice(idx0) +
                        c1 * (pc.slice(idx1).t() * pc.slice(idx1) +
                              (cc * (2 - cc)) * C.slice(idx0));
      }
    }

    if (iterate.n_rows > iterate.n_cols) {
      for (size_t j = 0; j < mu; ++j) {
        C.slice(idx1) = C.slice(idx1) + cmu * w(j) * pStep.slice(idx(j)) *
                                            pStep.slice(idx(j)).t();
      }
    } else {
      for (size_t j = 0; j < mu; ++j) {
        C.slice(idx1) = C.slice(idx1) + cmu * w(j) * pStep.slice(idx(j)).t() *
                                            pStep.slice(idx(j));
      }
    }

    arma::eig_sym(eigval, eigvec, C.slice(idx1));
    const arma::uvec negativeEigval = find(eigval < 0, 1);
    if (!negativeEigval.is_empty()) {
      if (negativeEigval(0) == 0) {
        C.slice(idx1).zeros();
      } else {
        C.slice(idx1) = eigvec.cols(0, negativeEigval(0) - 1) *
                        arma::diagmat(eigval.subvec(0, negativeEigval(0) - 1)) *
                        eigvec.cols(0, negativeEigval(0) - 1).t();
      }
    }

    // Output current objective function.
    cout << "CMA-ES: iteration " << i << ", objective " << overallobj << " "
         << error["frc"] << " " << error["punish"] << " " << error["shift"]
         << " cs " << sigma(idx1) << " " << ominrho << " " << omaxrho << " "
         << funcs[EMF].xx.front() << " " << funcs[EMF].xx.back() << " "
         << (lastobj - overallobj) / lastobj << " " << configs[locstt].fitengy
         << " " << configs[locstt].engy << endl;

    if (std::isnan(overallobj) || std::isinf(overallobj)) {
      cout << "CMA-ES: converged to " << overallobj << "; "
           << "terminating with failure.  Try a smaller step size?"
           << std::endl;
      return overallobj;
    }

    if (std::abs(lastobj - overallobj) < tolerance && i > 5000) {
      cout << "CMA-ES: minimized within tolerance " << tolerance << "; "
           << "terminating optimization." << std::endl;
      return overallobj;
    }

    // if (i % iparams["resfreq"] == 1 && rescaleEMF(iterate)) {
    //   cout << "before rescale " << currentobj << endl;
    //   currentobj = (this->*calobj[sparams["ptype"]])(decodev(iterate), 1);
    //   cout << "after rescale " << currentobj << endl;
    //   overallobj = currentobj;
    // }
    lastobj = overallobj;
  }

  return overallobj;
}