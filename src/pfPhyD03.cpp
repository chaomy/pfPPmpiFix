/*
 * @Author: chaomy
 * @Date:   2018-01-16 13:35:16
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-16 15:04:09
 */

#include "pfHome.h"

using std::vector;

void pfHome::buildD03(const string& kk, const double& gs, const double& dl) {
  int pn = 10;
  double lo = gs - pn * dl;
  Config cc;
  for (int i = 0; i < 2 * pn + 1; i++) {
    cc = buildD03(lo + i * dl);
    char buff[10];
    sprintf(buff, "poscar.%03d", i);
    writePOSCAR(cc, string(buff));
  }
}

Config pfHome::buildD03(const double& la) {
  Config cc;
  cc.bvx[X] = la, cc.bvx[Y] = 0.0, cc.bvx[Z] = 0.0;
  cc.bvy[X] = 0.0, cc.bvy[Y] = la, cc.bvy[Z] = 0.0;
  cc.bvz[X] = 0.0, cc.bvz[Y] = 0.0, cc.bvz[Z] = la;

  int cn = 0;
  vector<vector<double>> bs({{0.0, 0.0, 0.0},
                             {0.0, 0.5, 0.5},
                             {0.5, 0.0, 0.5},
                             {0.5, 0.5, 0.0},
                             {0.25, 0.25, 0.75},
                             {0.25, 0.75, 0.25},
                             {0.75, 0.25, 0.25},
                             {0.25, 0.25, 0.25},
                             {0.75, 0.75, 0.25},
                             {0.75, 0.25, 0.75},
                             {0.25, 0.75, 0.75},
                             {0.75, 0.75, 0.75},
                             {0.5, 0.0, 0.0},
                             {0.0, 0.5, 0.0},
                             {0.0, 0.0, 0.5},
                             {0.5, 0.5, 0.5}});
  cc.natomsv = vector<int>({12, 4});
  cc.nelemsv = vector<string>({"Mg", "Nd"});
  for (int ix : {0}) {
    for (int iy : {0}) {
      for (int iz : {0}) {
        for (auto vv : bs) {
          pfAtom atm(cn++);
          atm.pst[X] = (ix + vv[X]) * cc.bvx[X] + (iy + vv[X]) * cc.bvy[X] +
                       (iz + vv[X]) * cc.bvz[X];
          atm.pst[Y] = (ix + vv[Y]) * cc.bvx[Y] + (iy + vv[Y]) * cc.bvy[Y] +
                       (iz + vv[Y]) * cc.bvz[Y];
          atm.pst[Z] = (ix + vv[Z]) * cc.bvx[Z] + (iy + vv[Z]) * cc.bvy[Z] +
                       (iz + vv[Z]) * cc.bvz[Z];
          atm.tp = cn <= 12 ? 1 : 2;
          cc.atoms.push_back(atm);
        }
      }
    }
  }

  cc.natoms = cc.atoms.size();
  initBox(cc);
  initNeighsFull(cc);
  initAngles(cc);
  return cc;
}