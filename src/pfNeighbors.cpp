/*
 * @Author: chaomy
 * @Date:   2017-12-05 10:49:18
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-23 12:06:00
 */

#include "pfHome.h"

using std::cout;
using std::endl;
using std::vector;

void pfHome::initAnglesSameCutOff() {
  for (int cc = 0; cc < nconfs; cc++) {
    Config& tmpc = configs[cc];

    for (int ii = 0; ii < tmpc.natoms; ii++) {
      pfAtom& atmii = tmpc.atoms[ii];

      for (int jj = 0; jj < atmii.nneighsFull; jj++) {
        Neigh& ngbj = atmii.neighsFull[jj];
        vector<Angle> angvec;
        if (ngbj.r < funcs[RHO].xx.back()) {
          for (int kk = 0; kk < jj; kk++) {
            Neigh& ngbk = atmii.neighsFull[kk];
            Angle tmpang;
            if (ngbk.r < funcs[RHO].xx.back()) {
              tmpang.gcos = ngbj.dist2r[X] * ngbk.dist2r[X] +
                            ngbj.dist2r[Y] * ngbk.dist2r[Y] +
                            ngbj.dist2r[Z] * ngbk.dist2r[Z];
              setAngleslot(tmpang, funcs[MEAMG], tmpang.gcos);
            }  // rk < rcut
            angvec.push_back(tmpang);
          }  // kk
        }    // rj < rcut
        atmii.angMat.push_back(angvec);
      }  // jj
    }    // ii
  }      // cc
}

void pfHome::initNeighs() {
  ricut = rocut;
  ftn = 0;
  for (Config& cc : configs) {
    initBox(cc);
    initNeighs(cc);
  }
  cout << "ricut = " << ricut << endl;
}

void pfHome::initNeighs(Config& tmpc) {
  vector<double> d0(3);
  vector<double> dij(3);
  for (int ii = 0; ii < tmpc.natoms; ii++) {
    pfAtom& atmii = tmpc.atoms[ii];
    atmii.nneighs = 0;
    if (!atmii.neighs.empty()) atmii.neighs.clear();

    for (int jj = ii; jj < tmpc.natoms; jj++) {
      pfAtom& atmjj = tmpc.atoms[jj];

      for (int it : {0, 1, 2}) d0[it] = atmjj.pst[it] - atmii.pst[it];

      for (int ix = -tmpc.scale[0]; ix <= tmpc.scale[0]; ix++) {
        for (int iy = -tmpc.scale[1]; iy <= tmpc.scale[1]; iy++) {
          for (int iz = -tmpc.scale[2]; iz <= tmpc.scale[2]; iz++) {
            if ((ii == jj) && (ix == 0) && (iy == 0) && (iz == 0)) continue;

            for (int k : {0, 1, 2})
              dij[k] = d0[k] + ix * tmpc.bvx[k] + iy * tmpc.bvy[k] +
                       iz * tmpc.bvz[k];

            double r = sqrt(square33(dij));

            if (r < rocut) {
              ricut = r < ricut ? r : ricut;

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

              for (int i = 0; i < 2; i++) setNeighslotStd(tmpn, funcs[i], r);

              if (!sparams["ptype"].compare("ADP")) {  // ADP
                tmpn.dis2mat[XX] = tmpn.dist[X] * tmpn.dist[X];
                tmpn.dis2mat[YY] = tmpn.dist[Y] * tmpn.dist[Y];
                tmpn.dis2mat[ZZ] = tmpn.dist[Z] * tmpn.dist[Z];
                tmpn.dis2mat[XY] = tmpn.dist[X] * tmpn.dist[Y];
                tmpn.dis2mat[YZ] = tmpn.dist[Y] * tmpn.dist[Z];
                tmpn.dis2mat[ZX] = tmpn.dist[Z] * tmpn.dist[X];
                setNeighslotStd(tmpn, funcs[ADPU], r);
                setNeighslotStd(tmpn, funcs[ADPW], r);
              }
              atmii.neighs.push_back(tmpn);
            }  // rcut
          }    // iz
        }      // iy
      }        // ix
    }          // jj
    atmii.nneighs = atmii.neighs.size();
  }  // ii
}