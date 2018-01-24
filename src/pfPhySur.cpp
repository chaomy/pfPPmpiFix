/*
 * @Author: chaomy
 * @Date:   2017-12-19 16:23:55
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-20 19:18:55
 */

#include "pfHome.h"

void pfHome::calSurf() {
  double la = exprs["lat"];
  printf("la = %f\n", la);
  Config s1 = buildsur100(la, "sur");
  writePOSCAR(s1);
  forceMEAM(s1);
  printf("%f\n", (s1.fitengy - ubcc.fitengy) * 0.5 * s1.natoms / (la * la));
  // printf("surf100 = %f\n", (s1.fitengy - ubcc.fitengy * s1.natoms) /
  //                              (exprs["lat"] * exprs["lat"]));
}

Config pfHome::buildsur100(const double& lat, const string& tag) { /* 100 */
  Config cc;
  int h = 24;
  int l = (!tag.compare("sur")) ? h - 10 : h;
  double lz = h * lat;

  cc.bvx[X] = lat, cc.bvx[Y] = 0.0, cc.bvx[Z] = 0.0;
  cc.bvy[X] = 0.0, cc.bvy[Y] = lat, cc.bvy[Z] = 0.0;
  cc.bvz[X] = 0.0, cc.bvz[Y] = 0.0, cc.bvz[Z] = lz;

  int cn = 0;
  for (int iz = 0; iz < l; iz++) {
    pfAtom atm1(cn++);
    atm1.prl[X] = atm1.prl[Y] = 0.0;
    atm1.prl[Z] = double(iz) / double(h);

    atm1.pst[X] = atm1.prl[X] * lat;
    atm1.pst[Y] = atm1.prl[Y] * lat;
    atm1.pst[Z] = atm1.prl[Z] * lz;
    cc.atoms.push_back(atm1);

    pfAtom atm2(cn++);
    atm2.prl[X] = atm2.prl[Y] = 0.5;
    atm2.prl[Z] = double(iz + 0.5) / double(h);

    atm2.pst[X] = atm2.prl[X] * lat;
    atm2.pst[Y] = atm2.prl[Y] * lat;
    atm2.pst[Z] = atm2.prl[Z] * lz;
    cc.atoms.push_back(atm2);
  }
  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}

Config pfHome::buildsur110(const double& lat, const string& tag) {
  Config cc;
  return cc;
}

Config pfHome::buildsur211(const double& lat, const string& tag) {
  Config cc;
  return cc;
}