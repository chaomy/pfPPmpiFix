/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-13 12:17:18
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

double pfHome::forceMEAM(const arma::mat &vv) {
  int cnt = 0;
  for (Func &ff : funcs)
    for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];

  funcs[PHI].g1.front() = vv[cnt++];
  funcs[EMF].g1.front() = vv[cnt++];
  funcs[MEAMG].g1.front() = vv[cnt++];
  funcs[EMF].g1.back() = vv[cnt++];
  funcs[MEAMG].g1.back() = vv[cnt++];

  for (int i : {0, 1, 2, 3, 4}) splineNe(funcs[i], gradRight[i]);

  error["frc"] = 0.0;
  for (Config &cnf : configs) {
    forceMEAM(cnf);
    for (pfAtom &atm : cnf.atoms)
      if (atm.nneighsFull == N5)
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] =
              atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it] - atm.frc[it];
          error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
        }
    error["frc"] += square11(cnf.fitengy - cnf.engy);
  }  // cc

  // calLat("bcc");  // cal phys
  // calElas();
  // calPV();
  return error["frc"];
  // return error["frc"] + error["lat"] + error["ela"];
  // error["pv"];  // err per atom
}

double pfHome::forceMEAM(const vector<double> &vv) {
  int cnt = 0;
  for (Func &ff : funcs)
    for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];

  funcs[PHI].g1.front() = vv[cnt++];
  funcs[EMF].g1.front() = vv[cnt++];
  funcs[MEAMG].g1.front() = vv[cnt++];
  funcs[EMF].g1.back() = vv[cnt++];
  funcs[MEAMG].g1.back() = vv[cnt++];

  for (int i : {0, 1, 2, 3, 4}) splineNe(funcs[i], gradRight[i]);

  error["frc"] = 0.0;
  for (Config &cnf : configs) {
    forceMEAM(cnf);
    for (pfAtom &atm : cnf.atoms)
      if (atm.nneighsFull == N5)
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] =
              atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it] - atm.frc[it];
          error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
        }
    error["frc"] += square11(cnf.fitengy - cnf.engy);
  }  // cc

  // calLat("bcc");  // cal phys
  // calElas();
  // calPV();
  return error["frc"];
  // return error["frc"] + error["lat"] + error["ela"];
  // error["pv"];  // err per atom
}

void pfHome::forceMEAM(Config &cnf) {
  cnf.phiengy = cnf.emfengy = 0.0;
  for (pfAtom &atm : cnf.atoms) { /* loop over atoms to reset values */
    atm.crho = atm.prho = 0.0;
    for (int it : {X, Y, Z})
      atm.phifrc[it] = atm.rhofrc[it] = atm.trifrc[it] = 0.0;
  }  // ii

  double e0 = 0.0, e0g = 0.0;
  splint(funcs[EMF], 0.0, e0, e0g);

  for (pfAtom &atm : cnf.atoms) { /* loop over atoms pairs, densities */
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];

      // pair (phi)
      splint(funcs[PHI], ngbj.slots[PHI], ngbj.shifts[PHI], ngbj.steps[PHI],
             ngbj.phi, ngbj.phig);
      cnf.phiengy += ngbj.phi;

      for (int it : {X, Y, Z}) atm.phifrc[it] += ngbj.dist2r[it] * ngbj.phig;

      // rho
      if (ngbj.r >= rhcut)
        spltra(funcs[RHO], ngbj.r, ngbj.rho, ngbj.rhog);
      else
        splint(funcs[RHO], ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.steps[RHO],
               ngbj.rho, ngbj.rhog);
      atm.crho += ngbj.rho;

      // partial sum
      if (ngbj.r >= rhcut)
        spltra(funcs[MEAMF], ngbj.r, ngbj.fval, ngbj.fgrad);
      else
        splint(funcs[MEAMF], ngbj.slots[MEAMF - 1], ngbj.shifts[MEAMF - 1],
               ngbj.steps[MEAMF - 1], ngbj.fval, ngbj.fgrad);
      for (int kk = 0; kk < jj; kk++) {
        Neigh &ngbk = atm.neighsFull[kk];
        Angle &agl = atm.angMat[jj][kk];
        splint(funcs[MEAMG], agl.slot, agl.shift, agl.step, agl.gval,
               agl.ggrad);
        atm.prho += ngbk.fval * ngbj.fval * agl.gval;
      }  // kk
    }    //

    atm.rho = atm.crho + atm.prho;

    double embE = 0.0;
    splint(funcs[EMF], atm.rho, embE, atm.gradF);
    cnf.emfengy += (embE - e0);

    double forces_i[3] = {0.0, 0.0, 0.0};
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];

      double f_rij = ngbj.fval;
      double f_rij_prime = ngbj.fgrad;

      double forces_j[3] = {0.0, 0.0, 0.0};
      for (int kk = 0; kk < jj; kk++) {
        Neigh &ngbk = atm.neighsFull[kk];

        double gcos = atm.angMat[jj][kk].gcos;
        double gval = atm.angMat[jj][kk].gval;
        double gprime = atm.angMat[jj][kk].ggrad;

        double f_rik = ngbk.fval;
        double f_rik_prime = ngbk.fgrad;

        double fij = -atm.gradF * gval * f_rik * f_rij_prime;
        double fik = -atm.gradF * gval * f_rij * f_rik_prime;

        double prefactor = atm.gradF * f_rij * f_rik * gprime;

        double prefactor_ij = prefactor / ngbj.r;
        double prefactor_ik = prefactor / ngbk.r;

        fij += prefactor_ij * gcos;
        fik += prefactor_ik * gcos;

        double fj[3], fk[3];

        for (int it : {X, Y, Z}) {
          fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
          forces_j[it] += fj[it];

          fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
          forces_i[it] -= fk[it];
          cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
        }
      }  // loop over kk
      for (int it : {X, Y, Z}) {
        atm.trifrc[it] -= forces_j[it];
        cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
      }
    }  // loop over jj
    for (int it : {X, Y, Z}) atm.trifrc[it] += forces_i[it];
  }                             // atm
  for (pfAtom &atm : cnf.atoms) /* eambedding forces */
    for (Neigh &ngb : atm.neighsFull) {
      double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
      for (int it : {X, Y, Z}) atm.rhofrc[it] += ngb.dist2r[it] * emf;
    }  // nn
  cnf.fitengy = (cnf.phiengy / 2. + cnf.emfengy) / cnf.natoms;
}
