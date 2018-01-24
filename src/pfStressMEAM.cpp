/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-19 15:35:59
 */

#include "pfHome.h"

void pfHome::stressMEAM(Config &cnf) {
  double tm[3];
  cnf.phiengy = cnf.emfengy = 0.0;
  for (int i : {0, 1, 2, 3, 4, 5}) cnf.strs[i] = 0;
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

      for (int it : {X, Y, Z}) {
        tm[it] = ngbj.dist2r[it] * ngbj.phig;
        atm.phifrc[it] += tm[it];
      }

      cnf.strs[XX] += 0.5 * ngbj.dist[X] * tm[X];
      cnf.strs[YY] += 0.5 * ngbj.dist[Y] * tm[Y];
      cnf.strs[ZZ] += 0.5 * ngbj.dist[Z] * tm[Z];
      cnf.strs[XY] += 0.5 * ngbj.dist[X] * tm[Y];
      cnf.strs[YZ] += 0.5 * ngbj.dist[Y] * tm[Z];
      cnf.strs[ZX] += 0.5 * ngbj.dist[Z] * tm[X];

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
        /* i-j */
        cnf.strs[XX] -= ngbj.dist[X] * fj[X];
        cnf.strs[YY] -= ngbj.dist[Y] * fj[Y];
        cnf.strs[ZZ] -= ngbj.dist[Z] * fj[Z];
        cnf.strs[XY] -= ngbj.dist[X] * fj[Y];
        cnf.strs[YZ] -= ngbj.dist[Y] * fj[Z];
        cnf.strs[ZX] -= ngbj.dist[Z] * fj[X];

        /* i-k */
        cnf.strs[XX] -= ngbk.dist[X] * fk[X];
        cnf.strs[YY] -= ngbk.dist[Y] * fk[Y];
        cnf.strs[ZZ] -= ngbk.dist[Z] * fk[Z];
        cnf.strs[XY] -= ngbk.dist[X] * fk[Y];
        cnf.strs[YZ] -= ngbk.dist[Y] * fk[Z];
        cnf.strs[ZX] -= ngbk.dist[Z] * fk[X];
      }  // loop over kk
      for (int it : {X, Y, Z}) {
        atm.trifrc[it] -= forces_j[it];
        cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
      }
    }  // loop over jj
    for (int it : {X, Y, Z}) atm.trifrc[it] += forces_i[it];
  }  // atm

  for (pfAtom &atm : cnf.atoms) /* eambedding forces */
    for (Neigh &ngb : atm.neighsFull) {
      double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
      for (int it : {X, Y, Z}) {
        tm[it] = ngb.dist2r[it] * emf;
        atm.rhofrc[it] += tm[it];
      }
      cnf.strs[XX] += 0.5 * ngb.dist[X] * tm[X];
      cnf.strs[YY] += 0.5 * ngb.dist[Y] * tm[Y];
      cnf.strs[ZZ] += 0.5 * ngb.dist[Z] * tm[Z];
      cnf.strs[XY] += 0.5 * ngb.dist[X] * tm[Y];
      cnf.strs[YZ] += 0.5 * ngb.dist[Y] * tm[Z];
      cnf.strs[ZX] += 0.5 * ngb.dist[Z] * tm[X];
    }  // nn

  for (int i = 0; i < 6; i++) cnf.strs[i] /= cnf.vol;
  cnf.fitengy = (cnf.phiengy / 2. + cnf.emfengy) / cnf.natoms;
}