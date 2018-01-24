/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-19 12:15:26
 */

#include "pfHome.h"

double pfHome::forceEAM(const arma::mat& vv) {
  error["frc"] = 0.0;
  error["punish"] = 0.0;
  omaxrho = 1e-10;
  ominrho = 1e10;
  int cnt = 0;
  for (Func& ff : funcs)
    for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];

  if (sparams["spline"] != "nat") {
    funcs[PHI].g1.front() = vv[cnt++];
    funcs[RHO].g1.front() = vv[cnt++];
    funcs[EMF].g1.front() = vv[cnt++];
    funcs[EMF].g1.back() = vv[cnt++];
  }

  for (int i : {0, 1, 2}) splineNe(funcs[i], gradRight[i]);

  for (int i : {0, 1, 2})
    for (int j = 1; j < funcs[i].npts - 1; j++)
      error["punish"] += square11(funcs[i].g1[j]);

  for (Config& cnf : configs) {
    forceEAM(cnf);
    for (pfAtom& atm : cnf.atoms) {
      for (int it : {X, Y, Z}) {
        atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] - atm.frc[it];
        error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }
    }
    error["frc"] += square11(cnf.fitengy - cnf.engy);
    omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
    ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
    if (cnf.rhomx > hirho) error["punish"] += square11(cnf.rhomx - hirho);
    if (cnf.rhomi < lorho) error["punish"] += square11(lorho - cnf.rhomi);
  }
  error["frc"] *= 1e2;
  error["punish"] *= 1e4;
  calLat("bcc");
  return (error["frc"] + error["punish"] + error["lat"]);
}

double pfHome::forceEAM(const vector<double>& vv) {
  int cnt = 0;
  for (Func& ff : funcs)
    for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];

  if (sparams["spline"] != "nat") {
    funcs[PHI].g1.front() = vv[cnt++];
    funcs[RHO].g1.front() = vv[cnt++];
    funcs[EMF].g1.front() = vv[cnt++];
    funcs[EMF].g1.back() = vv[cnt++];
  }

  for (int i : {0, 1, 2}) splineNe(funcs[i], gradRight[i]);
  error["frc"] = 0.0;
  error["punish"] = 0.0;

  for (Config& cnf : configs) {
    forceEAM(cnf);
    for (pfAtom& atm : cnf.atoms) {
      for (int it : {X, Y, Z}) {
        atm.fitfrc[it] = atm.phifrc[it] + atm.rhofrc[it] - atm.frc[it];
        error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }
    }
    error["frc"] += square11(cnf.fitengy - cnf.engy);
    if (cnf.rhomx > hirho) error["punish"] += square11(cnf.rhomx - hirho);
    if (cnf.rhomi < lorho) error["punish"] += square11(lorho - cnf.rhomi);
  }
  return 1e3 * (error["frc"] + error["punish"]);
}

void pfHome::forceEAM(Config& cnf) {
  cnf.phiengy = cnf.emfengy = 0.0;
  cnf.rhomx = -1e3;
  cnf.rhomi = 1e3;
  for (pfAtom& atm : cnf.atoms) { /* reset values */
    atm.crho = 0.0;
    for (int it : {X, Y, Z}) atm.phifrc[it] = atm.rhofrc[it] = 0.0;
  }  // ii

  for (pfAtom& atm : cnf.atoms) { /* atoms pairs densities */
    for (Neigh& ngb : atm.neighs) {
      splint(funcs[PHI], ngb.slots[PHI], ngb.shifts[PHI], ngb.steps[PHI],
             ngb.phi, ngb.phig);
      cnf.phiengy += ngb.phi;

      double tmp[3];
      for (int it : {X, Y, Z}) {
        tmp[it] = ngb.dist2r[it] * ngb.phig;
        atm.phifrc[it] += tmp[it];
        cnf.atoms[ngb.aid].phifrc[it] -= tmp[it];
      }

      splint(funcs[RHO], ngb.slots[RHO], ngb.shifts[RHO], ngb.steps[RHO],
             ngb.rho, ngb.rhog);
      atm.crho += ngb.rho;
      cnf.atoms[ngb.aid].crho += ngb.rho;
    }  // nn

    double embE = 0.0;

    cnf.rhomx = atm.crho > cnf.rhomx ? atm.crho : cnf.rhomx;
    cnf.rhomi = atm.crho < cnf.rhomi ? atm.crho : cnf.rhomi;

    if (atm.crho > funcs[EMF].xx.back())
      spltra(funcs[EMF], atm.crho, embE, atm.gradF);
    else if (atm.crho < funcs[EMF].xx.front())
      spltrai(funcs[EMF], atm.crho, embE, atm.gradF);
    else
      splint(funcs[EMF], atm.crho, embE, atm.gradF);
    cnf.emfengy += embE;
  }  // ii

  for (pfAtom& atm : cnf.atoms) /* eambedding forces */
    for (Neigh& ngb : atm.neighs) {
      double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
      double tmp[3];
      for (int it : {0, 1, 2}) {
        tmp[it] = ngb.dist2r[it] * emf;
        atm.rhofrc[it] += tmp[it];
        cnf.atoms[ngb.aid].rhofrc[it] -= tmp[it];
      }
    }  // nn
  cnf.fitengy = (cnf.phiengy + cnf.emfengy) / cnf.natoms;
}