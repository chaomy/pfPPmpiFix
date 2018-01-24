/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-16 16:44:32
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

double pfHome::forceMEAM(const vector<double>& vv, int tag) {
  vector<Func>& ffs = funcs;
  // vv to funcs
  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    for (int j = 0; j < ffs[i].npts; j++) ffs[i].yy[j] = vv[cnt++];
    ffs[i].g1.front() = vv[cnt++];
    ffs[i].g1.back() = vv[cnt++];
    spline(ffs[i], gradRight[i]);  // update splines
  }

  // FILE* ptfile = fopen("force.meam.txt", "w"); 
  double error = 0.0;
  for (int cc = 0; cc < nconfs; cc++) {
    Config& cnf = configs[cc];
    double fitengy = 0.0;

    /*------------ first loop over atoms to reset values  -------------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      atm.rho = 0.0;
      for (int it = 0; it < 3; it++) atm.fitfrc[it] = -atm.frc[it];   
    }  // ii

    /*----------- second loop over atoms pairs, and densities ---------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];

      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        if (ngbj.r < ffs[PHI].xx.back()) {
          double phi = 0, phigrad = 0;
          splint(ffs[PHI], ngbj.slots[PHI], ngbj.shifts[PHI], ngbj.steps[PHI],
                 phi, phigrad);
          fitengy += phi;
          for (int it = 0; it < 3; it++)
            atm.fitfrc[it] += ngbj.dist2r[it] * phigrad;
        }  // fhi

        if (ngbj.r < ffs[RHO].xx.back()) {
          double rho = 0.0;
          splint(ffs[RHO], ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.steps[RHO],
                 rho, ngbj.rhog);
          atm.rho += rho;
        }  // rho

        if (ngbj.r < ffs[MEAMF].xx.back()) {
          splint(ffs[MEAMF], ngbj.slots[MEAMF - 1], ngbj.shifts[MEAMF - 1],
                 ngbj.steps[MEAMF - 1], ngbj.fval, ngbj.fgrad);
          double partial_sum = 0.0;

          for (int kk = 0; kk < jj; kk++) {
            Neigh& ngbk = atm.neighsFull[kk];
            Angle& agl = atm.angMat[jj][kk];
            splint(ffs[MEAMG], agl.slot, agl.shift, agl.step, agl.gval, 
                   agl.ggrad);
            partial_sum += ngbk.fval * ngbj.fval * agl.gval;
          }  // kk
          atm.rho += partial_sum;
        }  // fij * fik * g
      }    // jj
      
      double embE = 0.0;
      splint(ffs[EMF], atm.rho, embE, atm.gradF);
      fitengy += embE;

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

            cnf.atoms[ngbk.aid].fitfrc[it] += fk[it];
          }
        }  // loop over kk

        for (int it = 0; it < 3; it++) {
          atm.fitfrc[it] -= forces_j[it];
          cnf.atoms[ngbj.aid].fitfrc[it] += forces_j[it];
        }
      }  // loop over jj
      for (int it = 0; it < 3; it++) atm.fitfrc[it] += forces_i[it];
    }  // ii

    /*----------- last loop over atoms eambedding forces  ------------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];

      for (int nn = 0; nn < atm.nneighsFull; nn++) {
        Neigh& ngb = atm.neighsFull[nn];

        if (ngb.r < ffs[RHO].xx.back()) {
          double rhoj = ngb.rhog;
          double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF;

          for (int it = 0; it < 3; it++) atm.fitfrc[it] += ngb.dist2r[it] * emf;
        }  // r < rcut
      }    // nn
      // scaleVec(atm.fitfrc, 1. / (FRCEPS + atm.absfrc));
      error += square11(atm.fitfrc[0] * atm.fweigh[0]);
      error += square11(atm.fitfrc[1] * atm.fweigh[1]);
      error += square11(atm.fitfrc[2] * atm.fweigh[2]);

      /* benchmark use */
      // fprintf(ptfile, "%d 1 %.6g %.6g %.6g %.6g %.6g %.6g\n", ii+1, atm.pst[X],atm.pst[Y],atm.pst[Z],
      //                                                    atm.fitfrc[0], atm.fitfrc[1], atm.fitfrc[2]);
    }  // ii

    fitengy /= cnf.natoms;
    // error += dparams["eweight"] * square11(fitengy - cnf.engy);
    error += square11(fitengy - cnf.engy);
  }  // cc

  punish = 0.0;
  physic = 0.0;
  // fclose(ptfile);

  /* not add yet */
  // if (ffs[PHI].g1.front() > 0.0)
  //   punish += error * square11(ffs[PHI].g1.front());
  // if (ffs[RHO].g1.front() > 0.0)
  //   punish += error * square11(ffs[RHO].g1.front());
  // if (ffs[EMF].g1.front() > 0.0)
  //   punish += error * square11(ffs[EMF].g1.front());

  /* too slow to be included  */
  // writeLMPS(vv);
  // lmpdrv->error["tol"] = 0.0;
  // lmpdrv->calLattice();
  // lmpdrv->calElastic();
  // physic = lmpdrv->error["tol"];
  // punish *= dparams["pratio"];
  // error += physic;  
  error += punish;
  return error;
}