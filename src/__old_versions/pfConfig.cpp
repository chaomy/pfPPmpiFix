/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-04 09:35:55
 */

#include "pfHome.h"

using namespace MMatrix;
using std::cout;
using std::endl;
using std::vector;

inline void printVec(double v[3]) {
  for (int i = 0; i < 3; i++) cout << v[i] << " ";
  cout << endl;
}

void pfHome::initBox() {
  for (int cc = 0; cc < nconfs; cc++) {
    Config& tmpc = configs[cc];

    Mmatrix<double>::crossProd33(tmpc.bvy, tmpc.bvz, tmpc.tvx);
    Mmatrix<double>::crossProd33(tmpc.bvz, tmpc.bvx, tmpc.tvy);
    Mmatrix<double>::crossProd33(tmpc.bvx, tmpc.bvy, tmpc.tvz);

    tmpc.vol = Mmatrix<double>::vecInnProd33(tmpc.bvx, tmpc.tvx);
    double inv = 1. / tmpc.vol;

    // normalize
    scaleVec(tmpc.tvx, inv);
    scaleVec(tmpc.tvy, inv);
    scaleVec(tmpc.tvz, inv);

    double iheight[3];
    iheight[0] = sqrt(square33(tmpc.tvx));
    iheight[1] = sqrt(square33(tmpc.tvy));
    iheight[2] = sqrt(square33(tmpc.tvz));

    tmpc.scale[0] = (int)ceil(rocut * iheight[0]);
    tmpc.scale[1] = (int)ceil(rocut * iheight[1]);
    tmpc.scale[2] = (int)ceil(rocut * iheight[2]);
  }
}

void pfHome::initNeighs() {
  for (int cc = 0; cc < nconfs; cc++) {
    Config& tmpc = configs[cc];
    double d0[3];
    double dij[3];
    for (int ii = 0; ii < tmpc.natoms; ii++) {
      pfAtom& atmii = tmpc.atoms[ii];

      for (int jj = ii; jj < tmpc.natoms; jj++) {
        pfAtom& atmjj = tmpc.atoms[jj];

        d0[0] = atmjj.pst[0] - atmii.pst[0];
        d0[1] = atmjj.pst[1] - atmii.pst[1];
        d0[2] = atmjj.pst[2] - atmii.pst[2];

        for (int ix = -tmpc.scale[0]; ix <= tmpc.scale[0]; ix++) {
          for (int iy = -tmpc.scale[1]; iy <= tmpc.scale[1]; iy++) {
            for (int iz = -tmpc.scale[2]; iz <= tmpc.scale[2]; iz++) {
              if ((ii == jj) && (ix == 0) && (iy == 0) && (iz == 0)) continue;

              dij[0] = d0[0] + ix * tmpc.bvx[0] + iy * tmpc.bvy[0] +
                       iz * tmpc.bvz[0];
              dij[1] = d0[1] + ix * tmpc.bvx[1] + iy * tmpc.bvy[1] +
                       iz * tmpc.bvz[1];
              dij[2] = d0[2] + ix * tmpc.bvx[2] + iy * tmpc.bvy[2] +
                       iz * tmpc.bvz[2];

              double r = sqrt(square33(dij));

              if (r <= rocut) {
                Neigh tmpn;
                double invr = 1. / r;

                tmpn.r = r;
                tmpn.invr = invr;
                tmpn.aid = atmjj.id;

                tmpn.dist[0] = dij[0];
                tmpn.dist[1] = dij[1];
                tmpn.dist[2] = dij[2];

                tmpn.dist2r[0] = dij[0] * invr;
                tmpn.dist2r[1] = dij[1] * invr;
                tmpn.dist2r[2] = dij[2] * invr;

                // set neigh slot
                for (int i = 0; i < 2; i++) setNeighslot(tmpn, funcs[i], r);

                // ADP
                if (!sparams["ptype"].compare("ADP")) {
                  tmpn.dis2mat[XX] = tmpn.dist[X] * tmpn.dist[X];
                  tmpn.dis2mat[YY] = tmpn.dist[Y] * tmpn.dist[Y];
                  tmpn.dis2mat[ZZ] = tmpn.dist[Z] * tmpn.dist[Z];
                  tmpn.dis2mat[XY] = tmpn.dist[X] * tmpn.dist[Y];
                  tmpn.dis2mat[YZ] = tmpn.dist[Y] * tmpn.dist[Z];
                  tmpn.dis2mat[ZX] = tmpn.dist[Z] * tmpn.dist[X];
                  setNeighslot(tmpn, funcs[ADPU], r);
                  setNeighslot(tmpn, funcs[ADPW], r);
                }
                atmii.neighs.push_back(tmpn);
              }  // rcut
            }    // iz
          }      // iy
        }        // ix
      }          // jj
      atmii.nneighs = atmii.neighs.size();
    }  // ii
  }    // cc
}

void pfHome::initAngles() {
  // FILE* ptfile = fopen("angle.txt", "w");
  for (int cc = 0; cc < nconfs; cc++) {
    Config& tmpc = configs[cc];

    for (int ii = 0; ii < tmpc.natoms; ii++) {
      pfAtom& atmii = tmpc.atoms[ii];

      for (int jj = 0; jj < atmii.nneighsFull; jj++) {
        Neigh& ngbj = atmii.neighsFull[jj];
        vector<Angle> angvec;

        for (int kk = 0; kk < jj; kk++) {
          Neigh& ngbk = atmii.neighsFull[kk];
          Angle tmpang;
          tmpang.gcos = ngbj.dist2r[X] * ngbk.dist2r[X] +
                        ngbj.dist2r[Y] * ngbk.dist2r[Y] +
                        ngbj.dist2r[Z] * ngbk.dist2r[Z];
          setAngleslot(tmpang, funcs[MEAMG], tmpang.gcos);
          angvec.push_back(tmpang);
        }  // kk

        atmii.angMat.push_back(angvec);
      }  // jj
    }    // ii
  }      // cc
  // fclose(ptfile);
}

void pfHome::initNeighsFull() {
  for (int cc = 0; cc < nconfs; cc++) {
    Config& tmpc = configs[cc];
    double d0[3];
    double dij[3];
    for (int ii = 0; ii < tmpc.natoms; ii++) {
      pfAtom& atmii = tmpc.atoms[ii];

      for (int jj = 0; jj < tmpc.natoms; jj++) {
        pfAtom& atmjj = tmpc.atoms[jj];

        d0[0] = atmjj.pst[0] - atmii.pst[0];
        d0[1] = atmjj.pst[1] - atmii.pst[1];
        d0[2] = atmjj.pst[2] - atmii.pst[2];

        for (int ix = -tmpc.scale[0]; ix <= tmpc.scale[0]; ix++) {
          for (int iy = -tmpc.scale[1]; iy <= tmpc.scale[1]; iy++) {
            for (int iz = -tmpc.scale[2]; iz <= tmpc.scale[2]; iz++) {
              if ((ii == jj) && (ix == 0) && (iy == 0) && (iz == 0)) continue;

              dij[0] = d0[0] + ix * tmpc.bvx[0] + iy * tmpc.bvy[0] +
                       iz * tmpc.bvz[0];
              dij[1] = d0[1] + ix * tmpc.bvx[1] + iy * tmpc.bvy[1] +
                       iz * tmpc.bvz[1];
              dij[2] = d0[2] + ix * tmpc.bvx[2] + iy * tmpc.bvy[2] +
                       iz * tmpc.bvz[2];

              double r = sqrt(square33(dij));

              if (r < rocut) {
                Neigh tmpn;
                double invr = 1. / r;

                tmpn.r = r;
                tmpn.invr = invr;
                tmpn.aid = atmjj.id;

                tmpn.dist[X] = dij[X];
                tmpn.dist[Y] = dij[Y];
                tmpn.dist[Z] = dij[Z];

                tmpn.dist2r[X] = dij[X] * invr;
                tmpn.dist2r[Y] = dij[Y] * invr;
                tmpn.dist2r[Z] = dij[Z] * invr;

                // set neigh slot
                for (int i = 0; i < 2; i++) setNeighslot(tmpn, funcs[i], r);

                // MEAM
                if (!sparams["ptype"].compare("MEAM"))
                  setNeighslot(tmpn, funcs[MEAMF], r);

                atmii.neighsFull.push_back(tmpn);
              }  // rcut
            }    // iz
          }      // iy
        }        // ix
      }          // jj
      atmii.nneighsFull = atmii.neighsFull.size();
    }  // ii
  }    // cc
}

void pfHome::setNeighslot(Neigh& refn, Func func, double r) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  int slot = 0;

  while ((hi - lo) > 1) {
    slot = (hi + lo) >> 1;
    if (func.xx[slot] > r)
      hi = slot;
    else
      lo = slot;
  }

  slot = lo;
  double step = func.xx[hi] - func.xx[lo];
  double shift = (r - func.xx[lo]) / step;

  if (slot >= (func.npts - 1)) {
    slot--;
    shift += 1;
  }

  refn.slots.push_back(slot);
  refn.steps.push_back(step);
  refn.shifts.push_back(shift);
}

void pfHome::updateNeighslot(Neigh& refn, Func func, double r, int id) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  int slot = 0;

  while ((hi - lo) > 1) {
    slot = (hi + lo) >> 1;
    if (func.xx[slot] > r)
      hi = slot;
    else
      lo = slot;
  }

  slot = lo;
  double step = func.xx[hi] - func.xx[lo];
  double shift = (r - func.xx[lo]) / step;

  if (slot >= (func.npts - 1)) {
    slot--;
    shift += 1;
  }

  refn.slots[id] = slot;
  refn.steps[id] = step;
  refn.shifts[id] = shift;
}

void pfHome::setAngleslot(Angle& refang, Func func, double r) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  int slot = 0;

  while ((hi - lo) > 1) {
    slot = (hi + lo) >> 1;
    if (func.xx[slot] > r)
      hi = slot;
    else
      lo = slot;
  }

  slot = lo;
  double step = func.xx[hi] - func.xx[lo];
  double shift = (r - func.xx[lo]) / step;

  if (slot >= (func.npts - 1)) {
    slot--;
    shift += 1;
  }

  refang.slot = slot;
  refang.step = step;
  refang.shift = shift;
}