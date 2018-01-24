/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-14 16:39:06
 */

#include "pfHome.h"

double pfHome::forceEAM(const vector<double>& vv, int tag) {
  vector<Func>& ffs = funcs;
  // vv to funcs 
  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    for (int j = 0; j < ffs[i].npts; j++) ffs[i].yy[j] = vv[cnt++];
    ffs[i].g1.front() = vv[cnt++];
    ffs[i].g1.back() = vv[cnt++];
    spline(ffs[i], gradRight[i]);  // update splines
  }

  double error = 0.0;
  FILE* ptfile = fopen("force.pf.txt", "w");

  for (int cc = 0; cc < nconfs; cc++) {
    Config& cnf = configs[cc];
    double fitengy = 0.0;

    /*------------ first loop over atoms to reset values  -------------*/
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      atm.rho = 0.0;
      for (int it = 0; it < 3; it++) atm.fitfrc[it] = 0.0; 
      // atm.fitfrc[0] = -atm.frc[0];
      // atm.fitfrc[1] = -atm.frc[1];
      // atm.fitfrc[2] = -atm.frc[2];
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
        if (ngb.r < ffs[RHO].xx.back()) {
          double rho = 0.0;
          splint(ffs[RHO], ngb.slots[RHO], ngb.shifts[RHO], ngb.steps[RHO], rho,
                 ngb.rhog);
          atm.rho += rho;
          cnf.atoms[ngb.aid].rho += rho;
        }  // rho
      }    // nn

      double embE = 0.0;
      double extra = 0.0;
      if (atm.rho > ffs[EMF].xx.back()) {
        extra = atm.rho - ffs[EMF].xx.back();
        double xval = 0.5 * (ffs[EMF].xx[ffs[EMF].npts-1] 
                           + ffs[EMF].xx[ffs[EMF].npts-2]);
        splint(ffs[EMF], xval, embE, atm.gradF);
        embE += atm.gradF * extra;
      }
      else if (atm.rho < ffs[EMF].xx.front()) {
        extra = atm.rho - ffs[EMF].xx.front();
        splint(ffs[EMF], ffs[EMF].xx.front(), embE,
               atm.gradF);
        embE += atm.gradF * extra;
      } else 
        splint(ffs[EMF], atm.rho, embE, atm.gradF);

      fitengy += embE;
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
        } // r < rcut
      }  // nn
      // scaleVec(atm.fitfrc, 1. / (FRCEPS + atm.absfrc));
      error += square11(atm.fitfrc[0] * atm.fweigh[0]);
      error += square11(atm.fitfrc[1] * atm.fweigh[1]);
      error += square11(atm.fitfrc[2] * atm.fweigh[2]);

      fprintf(ptfile, "%.6f %.6f %.6f %.6f %.6f %.6f\n", atm.pst[X],atm.pst[Y],atm.pst[Z],
                                                         atm.fitfrc[0], atm.fitfrc[1], atm.fitfrc[2]);
    }  // ii
    fitengy /= cnf.natoms; 
    // error += dparams["eweight"] * square11(fitengy - cnf.engy);
    error += square11(fitengy - cnf.engy);
  }  // cc

  fclose(ptfile);
  punish = 0.0; 
  // punish += error * (square11(ffs[EMF].xx.front() + ffs[EMF].xx.back() - omaxrho - ominrho));   
  if (ffs[PHI].g1.front() > 0.0) punish += error * square11(ffs[PHI].g1.front());
  if (ffs[RHO].g1.front() > 0.0) punish += error * square11(ffs[RHO].g1.front());
  if (ffs[EMF].g1.front() > 0.0) punish += error * square11(ffs[EMF].g1.front());
  if (ffs[EMF].g1.back() < 0.0) punish += error * square11(ffs[EMF].g1.back()); 

  punish *= dparams["pratio"];
  error += punish;
  return error;
}