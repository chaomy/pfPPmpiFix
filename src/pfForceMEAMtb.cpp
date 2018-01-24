/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-13 12:17:33
 */

#include "pfHome.h"

double pfHome::forceMEAMtb(const vector<double>& vv, int tag) {
  vector<Func>& ffs = funcs;

  ffs[chid].g1.front() = (ffs[chid].yy[1] - ffs[chid].yy[0]) / ffs[chid].step;

  if (chid == PHI || chid == RHO || chid == MEAMF)
    ffs[chid].g1.back() = 0;
  else
    ffs[chid].g1.back() =
        (ffs[chid].yy.back() - ffs[chid].yy[ffs[chid].npts - 2]) /
        ffs[chid].step;

  splineNe(ffs[chid], gradRight[chid]);

  if (chid == PHI)
    changePHI();
  else if (chid == RHO)
    changeRHO();
  else if (chid == EMF)
    changeEMF();
  else if (chid == MEAMF)
    changeMEAMF();
  else if (chid == MEAMG)
    changeMEAMG();

  double err = 0;
  for (int cc = 0; cc < nconfs; cc++) {
    Config& cnf = configs[cc];
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      for (int it = 0; it < 3; it++) {
        // update force
        atm.fitfrc[it] =
            -atm.frc[it] + atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it];
        err += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }
      cnf.fitengy = cnf.phiengy + cnf.emfengy;
      err += square11(cnf.fitengy / cnf.natoms - cnf.engy);
    }
  }
  printf("%d %f\n", chid, err);
  return err;
}

void pfHome::changePHI() {
  for (Config& cnf : configs) {
    /* clear phi */
    cnf.phiengy = 0.0;
    for (int ii = 0; ii < cnf.natoms; ii++)
      for (int it = 0; it < 3; it++) cnf.atoms[ii].phifrc[it] = 0.0;

    /* add phi */
    for (pfAtom& atm : cnf.atoms) {
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];
        splint(funcs[PHI], ngbj.slots[PHI], ngbj.shifts[PHI], ngbj.steps[PHI],
               ngbj.phi, ngbj.phig);
        cnf.phiengy += ngbj.phi;
        for (int it = 0; it < 3; it++)
          atm.phifrc[it] += ngbj.dist2r[it] * ngbj.phig;
      }  // jj
    }    // ii
  }      // cc
}

void pfHome::changeRHO() {
  for (Config& cnf : configs) {
    cnf.emfengy = 0.0;
    /*------------ clear rhos  -------------*/
    for (pfAtom& atm : cnf.atoms) {
      atm.crho = 0.0;
      for (int it = 0; it < 3; it++) {
        atm.rhofrc[it] = 0.0;
        atm.trifrc[it] = 0.0;
      }
    }  // ii

    /*----------- second loop over atoms densities ---------*/
    for (pfAtom& atm : cnf.atoms) {
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        if (ngbj.r >= funcs[RHO].xx.back())
          spltra(funcs[RHO], ngbj.r, ngbj.rho, ngbj.rhog);
        else
          splint(funcs[RHO], ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.steps[RHO],
                 ngbj.rho, ngbj.rhog);
        atm.crho += ngbj.rho;
      }  // jj
      atm.rho = atm.crho + atm.prho;

      double embE = 0.0;
      splint(funcs[EMF], atm.rho, embE, atm.gradF);
      cnf.emfengy += embE;

      double forces_i[3] = {0.0, 0.0, 0.0};
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        double f_rij = ngbj.fval;
        double f_rij_prime = ngbj.fgrad;

        double forces_j[3] = {0.0, 0.0, 0.0};
        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];

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
          for (int it = 0; it < 3; it++) {
            fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
            forces_j[it] += fj[it];

            fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
            forces_i[it] -= fk[it];
            cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
          }
        }  // loop over kk
        for (int it = 0; it < 3; it++) {
          atm.trifrc[it] -= forces_j[it];
          cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
        }
      }  // loop over jj
      for (int it = 0; it < 3; it++) atm.trifrc[it] += forces_i[it];
    }  // ii

    /*----------- last loop over atoms eambedding forces  ------------*/
    for (pfAtom& atm : cnf.atoms) {
      for (int nn = 0; nn < atm.nneighsFull; nn++) {
        Neigh& ngb = atm.neighsFull[nn];

        double rhoj = ngb.rhog;
        double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF;

        for (int it = 0; it < 3; it++) atm.rhofrc[it] += ngb.dist2r[it] * emf;
      }  // nn
    }    // ii
  }      // cc
}

void pfHome::changeEMF() {
  for (Config& cnf : configs) {
    cnf.emfengy = 0.0;
    for (pfAtom& atm : cnf.atoms) {
      for (int it = 0; it < 3; it++) {
        atm.rhofrc[it] = 0.0;
        atm.trifrc[it] = 0.0;
      }
    }  // ii

    /*----------- second loop over atoms pairs, and densities ---------*/
    for (pfAtom& atm : cnf.atoms) {
      double embE = 0.0;
      splint(funcs[EMF], atm.rho, embE, atm.gradF);
      cnf.emfengy += embE;

      double forces_i[3] = {0.0, 0.0, 0.0};

      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        double f_rij = ngbj.fval;
        double f_rij_prime = ngbj.fgrad;

        double forces_j[3] = {0.0, 0.0, 0.0};

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];

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

          for (int it = 0; it < 3; it++) {
            fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
            forces_j[it] += fj[it];

            fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
            forces_i[it] -= fk[it];

            cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
          }
        }  // loop over kk
        for (int it = 0; it < 3; it++) {
          atm.trifrc[it] -= forces_j[it];
          cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
        }
      }  // loop over jj
      for (int it = 0; it < 3; it++) atm.trifrc[it] += forces_i[it];
    }  // ii

    /*----------- last loop over atoms eambedding forces  ------------*/
    for (pfAtom& atm : cnf.atoms) {
      /* embedding forces */
      for (int nn = 0; nn < atm.nneighsFull; nn++) {
        Neigh& ngb = atm.neighsFull[nn];

        double rhoj = ngb.rhog;
        double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF;

        for (int it = 0; it < 3; it++) atm.rhofrc[it] += ngb.dist2r[it] * emf;
      }  // nn
    }    // ii
  }      // cc
}

void pfHome::changeMEAMF() {
  for (Config& cnf : configs) {
    cnf.emfengy = 0.0;
    for (pfAtom& atm : cnf.atoms) {
      atm.prho = 0.0;
      for (int it = 0; it < 3; it++) {
        atm.rhofrc[it] = 0.0;
        atm.trifrc[it] = 0.0;
      }
    }  // ii

    /*----------- second loop over atoms pairs, and densities ---------*/
    for (pfAtom& atm : cnf.atoms) {
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        // partial sum
        if (ngbj.r >= funcs[MEAMF].xx.back())
          spltra(funcs[MEAMF], ngbj.r, ngbj.fval, ngbj.fgrad);
        else
          splint(funcs[MEAMF], ngbj.slots[MEAMF - 1], ngbj.shifts[MEAMF - 1],
                 ngbj.steps[MEAMF - 1], ngbj.fval, ngbj.fgrad);

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];
          Angle& agl = atm.angMat[jj][kk];

          splint(funcs[MEAMG], agl.slot, agl.shift, agl.step, agl.gval,
                 agl.ggrad);
          atm.prho += ngbk.fval * ngbj.fval * agl.gval;
        }  // kk
      }    //
      atm.rho = atm.crho + atm.prho;

      double embE = 0.0;
      splint(funcs[EMF], atm.rho, embE, atm.gradF);
      cnf.emfengy += embE;

      double forces_i[3] = {0.0, 0.0, 0.0};

      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        double f_rij = ngbj.fval;
        double f_rij_prime = ngbj.fgrad;

        double forces_j[3] = {0.0, 0.0, 0.0};

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];

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

          for (int it = 0; it < 3; it++) {
            fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
            forces_j[it] += fj[it];

            fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
            forces_i[it] -= fk[it];

            cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
          }
        }  // loop over kk
        for (int it = 0; it < 3; it++) {
          atm.trifrc[it] -= forces_j[it];
          cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
        }
      }  // loop over jj
      for (int it = 0; it < 3; it++) atm.trifrc[it] += forces_i[it];
    }  // ii
    /*----------- last loop over atoms eambedding forces  ------------*/
    for (pfAtom& atm : cnf.atoms) {
      /* embedding forces */
      for (int nn = 0; nn < atm.nneighsFull; nn++) {
        Neigh& ngb = atm.neighsFull[nn];

        double rhoj = ngb.rhog;
        double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF;

        for (int it = 0; it < 3; it++) atm.rhofrc[it] += ngb.dist2r[it] * emf;
      }  // nn
    }    // ii
  }      // cc
}

void pfHome::changeMEAMG() {
  for (Config& cnf : configs) {
    cnf.emfengy = 0.0;
    for (pfAtom& atm : cnf.atoms) {
      atm.prho = 0.0;
      for (int it = 0; it < 3; it++) {
        atm.rhofrc[it] = 0.0;
        atm.trifrc[it] = 0.0;
      }
    }  // ii

    /*----------- second loop over atoms pairs, and densities ---------*/
    for (pfAtom& atm : cnf.atoms) {
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];
          Angle& agl = atm.angMat[jj][kk];

          splint(funcs[MEAMG], agl.slot, agl.shift, agl.step, agl.gval,
                 agl.ggrad);
          atm.prho += ngbk.fval * ngbj.fval * agl.gval;
        }  // kk
      }    //
      atm.rho = atm.crho + atm.prho;
      double embE = 0.0;
      splint(funcs[EMF], atm.rho, embE, atm.gradF);
      cnf.emfengy += embE;

      double forces_i[3] = {0.0, 0.0, 0.0};
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        double f_rij = ngbj.fval;
        double f_rij_prime = ngbj.fgrad;

        double forces_j[3] = {0.0, 0.0, 0.0};

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];

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

          for (int it = 0; it < 3; it++) {
            fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
            forces_j[it] += fj[it];

            fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
            forces_i[it] -= fk[it];

            cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
          }
        }  // loop over kk
        for (int it = 0; it < 3; it++) {
          atm.trifrc[it] -= forces_j[it];
          cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
        }
      }  // loop over jj
      for (int it = 0; it < 3; it++) atm.trifrc[it] += forces_i[it];
    }  // ii
    /*----------- last loop over atoms eambedding forces  ------------*/
    for (pfAtom& atm : cnf.atoms) {
      /* embedding forces */
      for (int nn = 0; nn < atm.nneighsFull; nn++) {
        Neigh& ngb = atm.neighsFull[nn];

        double rhoj = ngb.rhog;
        double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF;

        for (int it = 0; it < 3; it++) atm.rhofrc[it] += ngb.dist2r[it] * emf;
      }  // nn
    }    // ii
  }      // cc
}