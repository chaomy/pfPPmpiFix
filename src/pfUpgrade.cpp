/*
 * @Author: chaomy
 * @Date:   2017-11-16 17:13:01
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-13 12:19:22
 */

#include "pfHome.h"

using std::vector;

// id are PHI, RHO, F
void pfHome::upgrade(int id) {
  int cc = 0;
  for (int i = 0; i < nfuncs; i++)
    for (int j = 0; j < funcs[i].npts; j++) funcs[i].yy[j] = ini[cc++];

  // update gradients
  funcs[PHI].g1.front() = ini[cc++];
  funcs[EMF].g1.front() = ini[cc++];
  funcs[MEAMG].g1.front() = ini[cc++];
  funcs[EMF].g1.back() = ini[cc++];
  funcs[MEAMG].g1.back() = ini[cc++];

  for (int i = 0; i < nfuncs; i++) splineNe(funcs[i], gradRight[i]);

  Func& fnc = funcs[id];

  fnc.xx.push_back(fnc.xx.back());
  fnc.yy.push_back(fnc.yy.back());
  fnc.g1.push_back(fnc.g1.back());
  fnc.g2.push_back(fnc.g2.back());
  fnc.npts = fnc.xx.size();

  double delt = (fnc.xx.back() - fnc.xx.front()) / (fnc.npts - 1);

  /* update y values */
  for (int i = 0; i < fnc.npts - 1; i++)
    splint(fnc, fnc.xx.front() + i * delt, fnc.yy[i]);

  /* update x values */
  for (int i = 0; i < fnc.npts; i++) fnc.xx[i] = fnc.xx.front() + i * delt;

  /* ---------------------------------------------- */
  nvars = 0;
  ini.clear();
  lob.clear();
  hib.clear();
  startps.clear();

  for (int i = 0; i < nfuncs; i++)
    funcs[i].step = funcs[i].xx[1] - funcs[i].xx.front();

  // set up ini
  for (int i = 0; i < nfuncs; i++) {
    startps.push_back(nvars);
    for (int j = 0; j < funcs[i].npts; j++) {
      ini.push_back(funcs[i].yy[j]);
      lob.push_back(lol[i]);
      hib.push_back(hil[i]);
      nvars++;
    }
    // g1 front
    ini.push_back(funcs[i].g1.front());
    lob.push_back(lol.back());
    hib.push_back(hil.back());
    // g1 end
    ini.push_back(funcs[i].g1.back());
    lob.push_back(lol.back());
    hib.push_back(hil.back());
    nvars += 2;
  }

  // update neigh slots
  if (id == PHI || id == RHO || id == MEAMF) {
    for (int cc = 0; cc < nconfs; cc++) {
      Config& tmpc = configs[cc];
      for (int ii = 0; ii < tmpc.natoms; ii++) {
        pfAtom& atmii = tmpc.atoms[ii];
        for (int jj = 0; jj < atmii.nneighsFull; jj++) {
          Neigh& ngbj = atmii.neighsFull[jj];
          if (id == PHI || id == RHO)
            updateNeighslot(ngbj, funcs[id], ngbj.r, id);
          else
            updateNeighslot(ngbj, funcs[id], ngbj.r, id - 1);
        }  // jj
      }    // ii
    }      // cc
  }        // PHI || RHO || MEAMF

  if (id == MEAMG) {
    for (int cc = 0; cc < nconfs; cc++) {
      Config& tmpc = configs[cc];
      for (int ii = 0; ii < tmpc.natoms; ii++) {
        pfAtom& atmii = tmpc.atoms[ii];
        for (int jj = 0; jj < atmii.nneighsFull; jj++) {
          for (int kk = 0; kk < jj; kk++) {
            Angle& tmpang = atmii.angMat[jj][kk];
            setAngleslot(tmpang, funcs[MEAMG], tmpang.gcos);
          }  // kk
        }    // jj
      }      // ii
    }        // cc
  }          // MEAMF
}
