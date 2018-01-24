/*
 * @Author: chaomy
 * @Date:   2017-12-16 22:11:52
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-18 13:40:19
 */

#include "pfHome.h"

#define SQ3B2 0.8660254037844386

void pfHome::buildhcp(const string& kk, const double& gs, const double& dl) {
  double l = gs - 5 * dl;
  mpvc[kk].clear();
  mpcf[kk].clear();
  for (int i = 0; i < 11; i++) {  // 4.22
    double covera = 1.82699686083;
    double a = l + dl * i;
    mpvc[kk].push_back(a);
    mpcf[kk].push_back(buildhcp(a, a * covera));
  }
}

Config pfHome::buildhcp(const double& la, const double& lc) {
  Config cc;
  cc.bvx[X] = la, cc.bvx[Y] = 0.0, cc.bvx[Z] = 0.0;
  cc.bvy[X] = -0.5 * la, cc.bvy[Y] = SQ3B2 * la, cc.bvy[Z] = 0.0;
  cc.bvz[X] = 0.0, cc.bvz[Y] = 0.0, cc.bvz[Z] = lc;
  int cnt = 0;
  for (int ix : {0}) {
    for (int iy : {0}) {
      for (int iz : {0}) {
        pfAtom atm1(cnt++);  //
        atm1.pst[X] = ix * cc.bvx[X] + iy * cc.bvy[X] + iz * cc.bvz[X];
        atm1.pst[Y] = ix * cc.bvx[Y] + iy * cc.bvy[Y] + iz * cc.bvz[Y];
        atm1.pst[Z] = ix * cc.bvx[Z] + iy * cc.bvy[Z] + iz * cc.bvz[Z];
        cc.atoms.push_back(atm1);

        pfAtom atm2(cnt++);
        atm2.pst[X] = (ix + 1. / 3.) * la + (iy + 2. / 3.) * -0.5 * la;
        atm2.pst[Y] = (iy + 2. / 3.) * SQ3B2 * la;
        atm2.pst[Z] = (iz + 0.5) * lc;
        cc.atoms.push_back(atm2);
      }
    }
  }
  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  writePOSCAR(cc);
  return cc;
}