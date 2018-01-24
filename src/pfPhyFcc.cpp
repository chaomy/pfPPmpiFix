/*
 * @Author: chaomy
 * @Date:   2017-12-16 22:11:52
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-16 13:55:55
 */

#include "pfHome.h"

void pfHome::buildfcc(const string& kk, const double& gs, const double& dl) {
  double l = gs - 5 * dl;
  mpvc[kk].clear();
  mpcf[kk].clear();
  for (int i = 0; i < 11; i++) {
    mpvc[kk].push_back(l + i * dl);
    mpcf[kk].push_back(buildfccPrim(mpvc[kk].back()));
  }
}

Config pfHome::buildfccPrim(const double& lat) {
  Config cc;
  double ll = 0.5 * lat;
  cc.bvx[X] = 0., cc.bvx[Y] = ll, cc.bvx[Z] = ll;
  cc.bvy[X] = ll, cc.bvy[Y] = 0., cc.bvy[Z] = ll;
  cc.bvz[X] = ll, cc.bvz[Y] = ll, cc.bvz[Z] = 0.;
  pfAtom atm(0);
  atm.prl[X] = atm.prl[Y] = atm.prl[Z] = 0.0;
  atm.pst[X] = atm.pst[Y] = atm.pst[Z] = 0.0;
  cc.atoms.push_back(atm);
  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}

Config pfHome::buildfccConv(const double& lat) {
  Config cc;
  cc.bvx[X] = lat, cc.bvx[Y] = 0.0, cc.bvx[Z] = 0.0;
  cc.bvy[X] = 0.0, cc.bvy[Y] = lat, cc.bvy[Z] = 0.0;
  cc.bvz[X] = 0.0, cc.bvz[Y] = 0.0, cc.bvz[Z] = lat;
  int cn = 0;
  for (int ix : {0}) {
    for (int iy : {0}) {
      for (int iz : {0}) {
        pfAtom atm1(cn++);
        atm1.prl[X] = atm1.prl[Y] = atm1.prl[Z] = 0.0;
        atm1.pst[X] = ix * lat;
        atm1.pst[Y] = iy * lat;
        atm1.pst[Z] = iz * lat;
        cc.atoms.push_back(atm1);

        pfAtom atm2(cn++);
        atm2.prl[X] = 0.5;
        atm2.prl[Y] = 0.5;
        atm2.prl[Z] = 0.0;
        //
        atm2.pst[X] = (ix + 0.5) * lat;
        atm2.pst[Y] = (iy + 0.5) * lat;
        atm2.pst[Z] = (iz)*lat;
        cc.atoms.push_back(atm2);

        pfAtom atm3(cn++);
        atm3.prl[X] = 0.0;
        atm3.prl[Y] = 0.5;
        atm3.prl[Z] = 0.5;
        //
        atm3.pst[X] = (ix)*lat;
        atm3.pst[Y] = (iy + 0.5) * lat;
        atm3.pst[Z] = (iz + 0.5) * lat;
        cc.atoms.push_back(atm3);

        pfAtom atm4(cn++);
        atm4.prl[X] = 0.5;
        atm4.prl[Y] = 0.0;
        atm4.prl[Z] = 0.5;

        atm4.pst[X] = (ix + 0.5) * lat;
        atm4.pst[Y] = (iy)*lat;
        atm4.pst[Z] = (iz + 0.5) * lat;
        cc.atoms.push_back(atm4);
      }
    }
  }
  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}