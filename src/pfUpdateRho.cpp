
/*
 * @Author: chaomy
 * @Date:   2017-10-30 21:34:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-17 15:09:17
 */

#include "pfHome.h"

void pfHome::updaterhoMEAM(vector<double>& vv) {
  oaverho = 0.0;
  omaxrho = -1e3;
  ominrho = 1e3;
  int cnt = 0;
  for (Func& ff : funcs)
    for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];

  // front
  funcs[PHI].g1.front() = vv[cnt++];
  funcs[EMF].g1.front() = vv[cnt++];
  funcs[MEAMG].g1.front() = vv[cnt++];
  // back
  funcs[EMF].g1.back() = vv[cnt++];
  funcs[MEAMG].g1.back() = vv[cnt++];

  for (int i = 0; i < nfuncs; i++) splineNe(funcs[i], gradRight[i]);

  for (int cc = 0; cc < nconfs; cc++) {
    double tmpsum = 0.0;
    Config& cnf = configs[cc];
    for (pfAtom& atm : cnf.atoms) atm.rho = 0.0;
    for (pfAtom& atm : cnf.atoms) {
      for (int jj = 0; jj < atm.nneighsFull; jj++) {
        Neigh& ngbj = atm.neighsFull[jj];

        double rho = 0;
        if (ngbj.r >= rhcut)
          spltra(funcs[RHO], ngbj.r, rho, ngbj.rhog);
        else
          splint(funcs[RHO], ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.steps[RHO],
                 rho, ngbj.rhog);

        atm.rho += rho;
        if (ngbj.r >= rhcut)
          spltra(funcs[MEAMF], ngbj.r, ngbj.fval, ngbj.fgrad);
        else
          splint(funcs[MEAMF], ngbj.slots[MEAMF - 1], ngbj.shifts[MEAMF - 1],
                 ngbj.steps[MEAMF - 1], ngbj.fval, ngbj.fgrad);

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atm.neighsFull[kk];
          Angle& agl = atm.angMat[jj][kk];
          splint(funcs[MEAMG], agl.slot, agl.shift, agl.step, agl.gval,
                 agl.ggrad);
          atm.rho += ngbk.fval * ngbj.fval * agl.gval;
        }  // kk
      }    // jj
    }      // ii
    for (pfAtom& atm : cnf.atoms) {
      omaxrho = fmax(omaxrho, atm.rho);
      ominrho = fmin(ominrho, atm.rho);
      tmpsum += atm.rho;
    }
    oaverho += tmpsum / cnf.natoms;
  }
  oaverho /= nconfs;
}

void pfHome::updaterho(vector<double>& vv) {
  oaverho = 0.0;
  omaxrho = -1e3;
  ominrho = 1e3;

  int cnt = 0;
  // careful !!!
  for (int i = 0; i < nfuncs; i++) {
    Func& ff = funcs[i];
    if (i == PHI || i == RHO || i == MEAMF) {
      for (int j = 0; j < ff.npts - 1; j++) ff.yy[j] = vv[cnt++];
    } else {
      for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];
    }
    ff.s.set_points(ff.xx, ff.yy);
  }

  for (Config& cnf : configs) {
    double tmpsum = 0.0;
    for (pfAtom& atm : cnf.atoms) atm.crho = 0.0;
    for (pfAtom& atm : cnf.atoms) {
      for (Neigh& ngb : atm.neighs) {
        funcs[RHO].s.deriv(ngb.slots[RHO], ngb.shifts[RHO], ngb.rho, ngb.rhog);
        atm.crho += ngb.rho;
        cnf.atoms[ngb.aid].crho += ngb.rho;
      }
    }
    for (pfAtom& atm : cnf.atoms) {
      omaxrho = fmax(omaxrho, atm.crho);
      ominrho = fmin(ominrho, atm.crho);
      tmpsum += atm.crho;
    }
    oaverho += tmpsum / cnf.natoms;
  }
  oaverho /= nconfs;
  printf("fx0 = %f min = %f ; fxn = %f max = %f ; ave = %f \n",
         funcs[EMF].xx[0], ominrho, funcs[EMF].xx.back(), omaxrho, oaverho);
}
