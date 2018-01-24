/*
 * @Author: chaomy
 * @Date:   2017-12-17 00:25:20
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-18 20:21:32
 */

#include "pfHome.h"

Config pfHome::addvolm(const Config& cc, const double& dl) {
  vector<vector<double>> stmt({{1 + dl, 0, 0}, {0, 1 + dl, 0}, {0, 0, 1 + dl}});
  Config nc = addstrain(cc, stmt);
  return nc;
}

Config pfHome::addotho(const Config& cc, const double& dl) {
  vector<vector<double>> stmt(
      {{1 + dl, 0, 0}, {0, 1 - dl, 0}, {0, 0, 1. / (1 - dl * dl)}});
  Config nc = addstrain(cc, stmt);
  return nc;
}

Config pfHome::addmono(const Config& cc, const double& dl) {
  vector<vector<double>> stmt(
      {{1, dl, 0}, {dl, 1, 0}, {0, 0, 1. / (1 - dl * dl)}});
  Config nc = addstrain(cc, stmt);
  return nc;
}

Config pfHome::addstrain(Config cc, const vector<vector<double>>& ss) {
  double n1[DIM];
  double n2[DIM];
  double n3[DIM];
  
  for (int i : {0, 1, 2}) {  // supercell => strain * basis
    n1[i] = ss[i][X] * cc.bvx[X] + ss[i][Y] * cc.bvx[Y] + ss[i][Z] * cc.bvx[Z];
    n2[i] = ss[i][X] * cc.bvy[X] + ss[i][Y] * cc.bvy[Y] + ss[i][Z] * cc.bvy[Z];
    n3[i] = ss[i][X] * cc.bvz[X] + ss[i][Y] * cc.bvz[Y] + ss[i][Z] * cc.bvz[Z];
  }

  for (int i : {X, Y, Z}) {
    cc.bvx[i] = n1[i];
    cc.bvy[i] = n2[i];
    cc.bvz[i] = n3[i];
  }

  for (pfAtom& atm : cc.atoms)  // positions
    for (int i : {X, Y, Z})
      atm.pst[i] = atm.prl[X] * cc.bvx[i] + atm.prl[Y] * cc.bvy[i] +
                   atm.prl[Z] * cc.bvz[i];
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}