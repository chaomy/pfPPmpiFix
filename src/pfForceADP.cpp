/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-13 12:16:09
 */

#include "pfHome.h"

double pfHome::forceADP(const vector<double>& vv, int tag) {
  vector<Func>& ffs = funcs;
  // vv to funcs
  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    for (int j = 0; j < ffs[i].npts; j++) ffs[i].yy[j] = vv[cnt++];
    ffs[i].g1.front() = vv[cnt++];
    ffs[i].g1.back() = vv[cnt++];
    splineNe(ffs[i], gradRight[i]);  // update splines
  }

  FILE* ptfile = fopen("force.adp.loc", "w");

  double err = 0.0;
  for (int cc = 0; cc < nconfs; cc++) {
    Config& cnf = configs[cc];
    double fitengy = 0.0;

    /*------------ first loop over atoms to reset values  -------------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      atm.rho = 0.0;
      // for (int it = 0; it < 3; it++) atm.fitfrc[it] = -atm.frc[it];
      for (int it = 0; it < 3; it++) atm.fitfrc[it] = 0.0;
      for (int it = 0; it < 3; it++) atm.mu[it] = 0.0;
      for (int it = 0; it < 6; it++) atm.lambda[it] = 0.0;
    }  // ii

    /*----------- second loop over atoms pairs, and densities ---------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];

      for (int nn = 0; nn < atm.nneighs; nn++) {
        Neigh& ngb = atm.neighs[nn];
        if (ngb.r < ffs[PHI].xx.back()) {
          double phi = 0, phigrad = 0;
          splint(ffs[PHI], ngb.slots[PHI], ngb.shifts[PHI], ngb.steps[PHI], phi,
                 phigrad);
          fitengy += phi;
          double tmp[3];
          tmp[0] = ngb.dist2r[0] * phigrad;
          tmp[1] = ngb.dist2r[1] * phigrad;
          tmp[2] = ngb.dist2r[2] * phigrad;

          atm.fitfrc[0] += tmp[0];
          atm.fitfrc[1] += tmp[1];
          atm.fitfrc[2] += tmp[2];

          cnf.atoms[ngb.aid].fitfrc[0] -= tmp[0];
          cnf.atoms[ngb.aid].fitfrc[1] -= tmp[1];
          cnf.atoms[ngb.aid].fitfrc[2] -= tmp[2];
        }  // fhi

        // rho
        if (ngb.r < ffs[RHO].xx.back()) {
          double rho = 0.0;
          splint(ffs[RHO], ngb.slots[RHO], ngb.shifts[RHO], ngb.steps[RHO], rho,
                 ngb.rhog);
          atm.rho += rho;
          cnf.atoms[ngb.aid].rho += rho;
        }  // rho

        //  mu
        if (ngb.r < ffs[ADPU].xx.back()) {
          splint(ffs[ADPU], ngb.slots[ADPU - 1], ngb.shifts[ADPU - 1],
                 ngb.steps[ADPU - 1], ngb.uval, ngb.ugrad);
          for (int it = 0; it < 3; it++) {
            double tmpval = ngb.uval * ngb.dist[it];
            atm.mu[it] += tmpval;
            cnf.atoms[ngb.aid].mu[it] -= tmpval;
          }
        }  // mu

        // lambda
        if (ngb.r < ffs[ADPW].xx.back()) {
          splint(ffs[ADPW], ngb.slots[ADPW - 1], ngb.shifts[ADPW - 1],
                 ngb.steps[ADPW - 1], ngb.wval, ngb.wgrad);
          for (int it = 0; it < 6; it++) {
            double tmpval = ngb.wval * ngb.dis2mat[it];
            atm.lambda[it] += tmpval;
            cnf.atoms[ngb.aid].lambda[it] += tmpval;
          }
        }  // lambda
      }    // nn

      fprintf(ptfile, "%.6f %.6f %.6f %.6f %.6f %.6f\n", atm.pst[X], atm.pst[Y],
              atm.pst[Z], atm.fitfrc[X], atm.fitfrc[Y], atm.fitfrc[Z]);

      double embE = 0.0;
      // double extra = 0.0;
      // if (atm.rho > ffs[EMF].xx.back()) {
      //   extra = atm.rho - ffs[EMF].xx.back();
      //   double xval = 0.5 * (ffs[EMF].xx[ffs[EMF].npts - 1] +
      //                        ffs[EMF].xx[ffs[EMF].npts - 2]);
      //   splint(ffs[EMF], xval, embE, atm.gradF);
      //   embE += atm.gradF * extra;
      // } else if (atm.rho < ffs[EMF].xx.front()) {
      //   extra = atm.rho - ffs[EMF].xx.front();
      //   splint(ffs[EMF], ffs[EMF].xx.front(), embE, atm.gradF);
      //   embE += atm.gradF * extra;
      // } else
      splint(ffs[EMF], atm.rho, embE, atm.gradF);

      // energy
      double adpE = 0.0;
      for (int it = 0; it < 3; it++) adpE += square11(atm.mu[it]);
      atm.nu = atm.lambda[0] + atm.lambda[1] + atm.lambda[2];
      double trc = atm.nu / 3.;
      for (int it = 0; it < 3; it++) adpE += square11(atm.lambda[it] - trc);
      for (int it = 3; it < 6; it++) adpE += square11(atm.lambda[it]) * 2.0;
      fitengy += (embE + 0.5 * adpE);
    }  // ii

    /*----------- third loop over atoms eambedding forces  ------------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];

      for (int nn = 0; nn < atm.nneighs; nn++) {
        Neigh& ngb = atm.neighs[nn];

        if (ngb.r < ffs[RHO].xx.back()) {
          double rhoj = ngb.rhog;
          double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF;
          double tmp[3];

          tmp[0] = ngb.dist2r[0] * emf;
          tmp[1] = ngb.dist2r[1] * emf;
          tmp[2] = ngb.dist2r[2] * emf;

          atm.fitfrc[0] += tmp[0];
          atm.fitfrc[1] += tmp[1];
          atm.fitfrc[2] += tmp[2];

          cnf.atoms[ngb.aid].fitfrc[0] -= tmp[0];
          cnf.atoms[ngb.aid].fitfrc[1] -= tmp[1];
          cnf.atoms[ngb.aid].fitfrc[2] -= tmp[2];
        }  // r < rcut

        // adpu
        if (ngb.r < ffs[ADPU].xx.back()) {
          double uforce[3];
          double adpF[3];
          for (int it = 0; it < 3; it++)
            uforce[it] = atm.mu[it] - cnf.atoms[ngb.aid].mu[it];

          double upxyz =
              ngb.ugrad * (uforce[0] * ngb.dist[0] + uforce[1] * ngb.dist[1] +
                           uforce[2] * ngb.dist[2]);

          for (int it = 0; it < 3; it++) {
            adpF[it] = uforce[it] * ngb.uval + upxyz * ngb.dist2r[it];
            atm.fitfrc[it] += adpF[it];
            cnf.atoms[ngb.aid].fitfrc[it] -= adpF[it];
          }
        }  // adpu

        // adpw
        if (ngb.r < ffs[ADPW].xx.back()) {
          double wforce[6];
          double adpw1[3], adpw2[3], adpF[3];
          double wval = ngb.wval;
          double wgrad = ngb.wgrad;

          for (int it = 0; it < 6; it++)
            wforce[it] = atm.lambda[it] + cnf.atoms[ngb.aid].lambda[it];

          double wpxyz = wgrad * (wforce[XX] * ngb.dis2mat[XX] +
                                  wforce[YY] * ngb.dis2mat[YY] +
                                  wforce[ZZ] * ngb.dis2mat[ZZ]);

          adpw1[X] =
              wval * wforce[XX] * 2.0 * ngb.dist[X] + ngb.dist2r[X] * wpxyz;
          adpw1[Y] =
              wval * wforce[YY] * 2.0 * ngb.dist[Y] + ngb.dist2r[Y] * wpxyz;
          adpw1[Z] =
              wval * wforce[ZZ] * 2.0 * ngb.dist[Z] + ngb.dist2r[Z] * wpxyz;

          // trace
          double nu = 1. / 3. * (atm.nu + cnf.atoms[ngb.aid].nu);
          double tmptrc = nu * (wval * 2.0 + wgrad * ngb.r);
          for (int it = 0; it < 3; it++) adpw1[it] -= tmptrc * ngb.dist[it];

          // xy yz zx
          double ppxyz = wgrad * (wforce[XY] * ngb.dis2mat[XY] +
                                  wforce[YZ] * ngb.dis2mat[YZ] +
                                  wforce[ZX] * ngb.dis2mat[ZX]);
          adpw2[X] =
              wval * (wforce[XY] * ngb.dist[Y] + wforce[ZX] * ngb.dist[Z]) +
              ngb.dist2r[X] * ppxyz;
          adpw2[Y] =
              wval * (wforce[XY] * ngb.dist[X] + wforce[YZ] * ngb.dist[Z]) +
              ngb.dist2r[Y] * ppxyz;
          adpw2[Z] =
              wval * (wforce[ZX] * ngb.dist[X] + wforce[YZ] * ngb.dist[Y]) +
              ngb.dist2r[Z] * ppxyz;

          for (int it = 0; it < 3; it++) {
            adpF[it] = adpw1[it] + 2. * adpw2[it];
            atm.fitfrc[it] += adpF[it];
            cnf.atoms[ngb.aid].fitfrc[it] -= adpF[it];
          }
        }  // adpw
      }    // nn
      // scaleVec(atm.fitfrc, 1. / (FRCEPS + atm.absfrc));
      err += square11(atm.fitfrc[0] * atm.fweigh[0]);
      err += square11(atm.fitfrc[1] * atm.fweigh[1]);
      err += square11(atm.fitfrc[2] * atm.fweigh[2]);
    }  // ii
    fitengy /= cnf.natoms;
    err += dparams["eweight"] * square11(fitengy - cnf.engy);
  }  // cc

  punish = 0.0;
  // punish += err * (square11(ffs[EMF].xx.front() + ffs[EMF].xx.back() -
  // omaxrho - ominrho));
  if (ffs[PHI].g1.front() > 0.0) punish += err * square11(ffs[PHI].g1.front());
  if (ffs[RHO].g1.front() > 0.0) punish += err * square11(ffs[RHO].g1.front());
  if (ffs[EMF].g1.front() > 0.0) punish += err * square11(ffs[EMF].g1.front());
  // if (ffs[EMF].g1.back() < 0.0) punish += err * square11(ffs[EMF].g1.back());
  punish *= dparams["pratio"];
  fclose(ptfile);
  err += punish;
  return err;
}