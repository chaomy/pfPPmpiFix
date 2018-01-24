/*
 * @Author: chaomy
 * @Date:   2017-10-30 21:34:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-18 14:11:56
 */

#include "pfHome.h"


void pfHome::updaterhoMEAM(vector<double>& vv) {
  oaverho = 0.0;
  omaxrho = -1e3;
  ominrho = 1e3;
  vector<Func>& ffs = funcs;
  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    for (int j = 0; j < ffs[i].npts; j++) ffs[i].yy[j] = vv[cnt++];
    ffs[i].g1.front() = vv[cnt++];
    ffs[i].g1.back() = vv[cnt++];
    spline(ffs[i], gradRight[i]);  // update splines
  }

  for (int cc = 0; cc < nconfs; cc++) {
    double tmpsum = 0.0;
    Config& cnf = configs[cc];
    for (int ii = 0; ii < cnf.natoms; ii++) cnf.atoms[ii].rho = 0.0;
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      for (int jj = 0; jj < atm.nneighsFull; jj++){
        Neigh& ngbj = atm.neighsFull[jj];

        if (ngbj.r < ffs[RHO].xx.back()){
          double rho = 0;
          splint(ffs[RHO], ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.steps[RHO], rho,
                 ngbj.rhog);
          atm.rho += rho; 
        } // rho

        if (ngbj.r < ffs[MEAMF].xx.back()){
          splint(ffs[MEAMF], ngbj.slots[MEAMF - 1], ngbj.shifts[MEAMF - 1],
                 ngbj.steps[MEAMF - 1], ngbj.fval, ngbj.fgrad); 
          double partial_sum = 0.0; 

          for (int kk = 0; kk < jj; kk++){
            Neigh& ngbk = atm.neighsFull[kk];
            splint(ffs[MEAMF], ngbk.slots[MEAMF - 1], ngbk.shifts[MEAMF - 1],
                   ngbk.steps[MEAMF - 1], ngbk.fval, ngbk.fgrad);
            Angle& agl = atm.angMat[jj][kk];
            splint(ffs[MEAMG], agl.slot, agl.shift, agl.step, agl.gval, 
                   agl.ggrad);
            partial_sum += ngbk.fval * ngbj.fval * agl.gval; 
          } // kk  
          atm.rho += partial_sum; 
        } // partial_sum 
      } // jj 
    } // ii

    for (int ii = 0; ii < cnf.natoms; ii++){
      pfAtom& atm = cnf.atoms[ii];
      omaxrho = fmax(omaxrho, atm.rho);
      ominrho = fmin(ominrho, atm.rho);
      tmpsum += atm.rho;
    }
    oaverho += tmpsum / cnf.natoms;
  }
  oaverho /= nconfs; 
  printf("min = %f ; max = %f ; ave = %f \n", ominrho, omaxrho, oaverho);
  return;
}

void pfHome::updaterho(vector<double>& vv) {
  oaverho = 0.0;
  omaxrho = -1e3;
  ominrho = 1e3;
  vector<Func>& ffs = funcs;
  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    for (int j = 0; j < ffs[i].npts; j++) ffs[i].yy[j] = vv[cnt++];
    ffs[i].g1.front() = vv[cnt++];
    ffs[i].g1.back() = vv[cnt++];
    spline(ffs[i], gradRight[i]);  // update splines
  }
  for (int cc = 0; cc < nconfs; cc++) {
    double tmpsum = 0.0;
    Config& cnf = configs[cc];
    for (int ii = 0; ii < cnf.natoms; ii++) cnf.atoms[ii].rho = 0.0;
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      for (int nn = 0; nn < atm.nneighs; nn++) {
        Neigh& ngb = atm.neighs[nn];
        if (ngb.r < ffs[RHO].xx.back()) {
          double rho = 0.0;
          splint(ffs[RHO], ngb.slots[RHO], ngb.shifts[RHO], ngb.steps[RHO], rho,
                 ngb.rhog);
          atm.rho += rho;
          cnf.atoms[ngb.aid].rho += rho;
        }  // rho
      }
    }
    for (int ii = 0; ii < cnf.natoms; ii++) {
      pfAtom& atm = cnf.atoms[ii];
      omaxrho = fmax(omaxrho, atm.rho);
      ominrho = fmin(ominrho, atm.rho);
      tmpsum += atm.rho;
    }
    oaverho += tmpsum / cnf.natoms;
  }
  oaverho /= nconfs;
  printf("min = %f ; max = %f ; ave = %f \n", ominrho, omaxrho, oaverho);
}
