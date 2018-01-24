/*
 * @Xuthor: chaomy
 * @Date:   2018-01-10 20:08:18
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-11 22:40:31
 */
#include "pfHome.h"

#define UPBND 10
#define LWBND -10

using arma::accu;
using arma::linspace;
using arma::mat;
using arma::randu;
using arma::vec;
using std::cout;
using std::endl;

double pfHome::testFunc(arma::mat& vc) {
  // printf("%f %f %f\n", vc[0], vc[1], vc[2]);
  // return -std::exp(-std::abs(vc[0])) + std::pow(vc[1], 2) + std::pow(vc[2],
  // 4) +
  //        3 * std::pow(vc[2], 2);
  return pow((vc[0] - 2.0), 2) + pow((vc[1] - 3.0), 2) + pow((vc[2] - 0.01), 2);
}

void pfHome::cmaes() {
  int nr = 3, nc = 1;
  arma::mat X(nr, nc, arma::fill::randu);

  cout << "X.n_rows: " << X.n_rows << " X.n_cols: " << X.n_cols
       << " X.n_elem: " << X.n_elem << endl;

  int n = X.n_elem;

  /* generation size */
  int lmbd = 4 + std::round(3 * std::log(n)) * 10;  // (48)
  const int mu = std::round(lmbd / 2);

  /* weights */
  vec w = std::log(mu + 0.5) -
          arma::log(linspace<vec>(0, mu - 1, mu) + 1.0);  // (49)
  w /= arma::sum(w);

  const double mue = 1. / accu(arma::pow(w, 2));

  /* step size control */
  vec sgm(3);
  sgm(0) = 0.3 * (UPBND - LWBND);

  const double cs = (mue + 2) / (n + mue + 5);  // (55)
  const double ds =
      1 + cs + 2 * std::max(std::sqrt((mue - 1) / (n + 1)) - 1, 0.0);  // (55)
  const double exnn =
      std::sqrt(n) * (1. - 1. / (4. * n) + 1. / (21 * std::pow(n, 2)));
  const double hh = (1.4 + 2. / (n + 1.)) * exnn;

  /* covariance matrix adaptation */
  const double cc = (4 + mue / n) / (n + 4 + 2 * mue / n);  // (56)
  const double alcov = 2.;
  const double c1 = alcov / (std::pow(n + 1.3, 2) + mue);  // (57)
  const double cmu =
      std::min(1 - c1, alcov * (mue - 2. + 1. / mue) /
                           (std::pow(n + 2, 2) + alcov * mue / 2.));  // (58)

  printf("n = %d lmbda = %d mu = %d\n", n, lmbd, mu);

  arma::cube mn(nr, nc, 3);
  mn.slice(0) = LWBND + arma::randu(nr, nc) * (UPBND - LWBND);
  arma::mat wy = arma::zeros(nr, nc);

  /* initialize */
  double cr = testFunc(mn.slice(0));
  double op = cr;
  double ls = 1e30;

  /* population */
  arma::cube yk(nr, nc, lmbd);
  arma::cube xn(nr, nc, lmbd);
  arma::cube ps = arma::zeros(nr, nc, 2);
  arma::cube pc = arma::zeros(nr, nc, 2);
  arma::cube C(n, n, 2);
  C.slice(0).eye();

  arma::vec err(lmbd);

  /* covariance matrix params */
  arma::vec egval;
  arma::mat egvec;
  arma::vec egvalz = arma::zeros(n);

  /* visitation order (sorted by population objectives). */
  arma::uvec idx = arma::linspace<arma::uvec>(0, lmbd - 1, lmbd);

  int maxIt = 20;
  double tol = 1e-18;

  for (int i = 1; i < maxIt; ++i) {
    const int vr0 = (i - 1) % 2;
    const int vr1 = i % 2;
    const arma::mat covL = arma::chol(C.slice(vr0), "lower");

    /* sample populations */
    for (int j = 0; j < lmbd; ++j) {
      yk.slice(idx(j)) = covL * arma::randn(nr, nc);
      xn.slice(idx(j)) = mn.slice(vr0) + sgm(vr0) * yk.slice(idx(j));
      err(idx(j)) = testFunc(xn.slice(idx(j)));
    }

    /* update mean */
    idx = arma::sort_index(err);
    wy = w(0) * yk.slice(idx(0));
    for (int j = 1; j < mu; ++j) wy += w(j) * yk.slice(idx(j));
    mn.slice(vr1) = mn.slice(vr0) + sgm(vr0) * wy;

    /* update best */
    cr = testFunc(mn.slice(vr1));
    if (cr < op) {
      op = cr;
      X = mn.slice(vr1);
      mn.slice(vr1).print("Mn = ");
    }

    /* control step size */
    ps.slice(vr1) = (1 - cs) * ps.slice(vr0) +
                    std::sqrt(cs * (2 - cs) * mue) * covL.t() * wy;

    const double psnm = arma::norm(ps.slice(vr1));
    sgm(vr1) = sgm(vr0) * std::pow(std::exp(cs / ds * psnm / exnn - 1), 0.3);

    /* update covariance matrix */
    if ((psnm / sqrt(1 - std::pow(1 - cs, 2 * i))) < hh) {
      pc.slice(vr1) =
          (1 - cc) * pc.slice(vr0) + std::sqrt(cc * (2 - cc) * mue) * wy;
      C.slice(vr1) = (1 - c1 - cmu) * C.slice(vr0) +
                     c1 * (pc.slice(vr1) * pc.slice(vr1).t());
    } else {
      pc.slice(vr1) = (1 - cc) * pc.slice(vr0);
      C.slice(vr1) = (1 - c1 - cmu) * C.slice(vr0) +
                     c1 * (pc.slice(vr1) * pc.slice(vr1).t() +
                           (cc * (2 - cc)) * C.slice(vr0));
    }
    for (int j = 0; j < mu; ++j)
      C.slice(vr1) =
          C.slice(vr1) + cmu * w(j) * yk.slice(idx(j)) * yk.slice(idx(j)).t();

    arma::eig_sym(egval, egvec, C.slice(vr1));
    const arma::uvec neg = find(egval < 0, 1);

    if (!neg.is_empty()) {
      if (neg(0) == 0)
        C.slice(vr1).zeros();
      else
        C.slice(vr1) = egvec.cols(0, neg(0) - 1) *
                       arma::diagmat(egval.subvec(0, neg(0) - 1)) *
                       egvec.cols(0, neg(0) - 1).t();
    }

    /* output */
    cout << "iteration " << i << " object " << op << endl;

    if (std::isnan(op) || std::isinf(op)) {
      cout << "converged to " << op << "try a smaller step size" << endl;
      return;
    }

    if (std::abs(ls - op) < tol) {
      cout << "minimized wiithin " << tol << endl;
      X.print();
      return;
    }
    ls = op;
  }
  return;
}