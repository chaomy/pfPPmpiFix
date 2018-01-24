/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 15:52:29
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-23 19:29:44
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
using std::cout;
using std::endl;

double pfHome::forceMEAM(const arma::mat &vv, int tg) {
  while (true) {
    broadcast(cmm, tg, PFROOT);
    if (tg == EXT) break;

    int cnt = 0;
    for (int i = 0; i < nfuncs; i++) { /* interpolates */
      Func &ff = funcs[i];
      double mxf = -1e10, mif = 1e10;
      int nt = (i == PHI || i == RHO || i == MEAMF) ? ff.npts - 1 : ff.npts;
      for (int j = 0; j < nt; j++) {
        ff.yy[j] = vv[cnt++];
        mxf = ff.yy[j] > mxf ? ff.yy[j] : mxf;
        mif = ff.yy[j] < mif ? ff.yy[j] : mif;
      }
      ff.rng = mxf - mif;
    }

    for (int i = 0; i < nfuncs; i++) {  // broadcast functions
      broadcast(cmm, funcs[i].xx, PFROOT);
      broadcast(cmm, funcs[i].yy, PFROOT);
      broadcast(cmm, funcs[i].rng, PFROOT);
    }

    for (Func &ff : funcs) ff.s.set_points(ff.xx, ff.yy);

    double efrc = 0.0, epsh = 0.0;
    error["frc"] = 0.0, error["punish"] = 0.0, error["shift"] = 0.0;
    omaxrho = -1e10, ominrho = 1e10;

    int ls[] = {PHI, RHO, MEAMF};
    for (int it : ls) {
      double invrg = 1. / square11(funcs[it].rng);
      double tm = 0.0;
      for (int i = 0; i < funcs[it].s.m_b.size() - 1; i++)
        tm += (square11(funcs[it].s.m_b[i]) +
               0.5 * funcs[it].s.m_b[i] * funcs[it].s.m_b[i + 1]);
      tm += square11(funcs[it].s.m_b.back());
      epsh += tm * invrg;
    }

    for (int i = locstt; i < locend; i++) {
      Config &cnf = configs[i];
      forceMEAM(cnf);
      for (pfAtom &atm : cnf.atoms) {
        for (int it : {X, Y, Z}) {
          atm.fitfrc[it] =
              atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it] - atm.frc[it];
          efrc += square11(atm.fitfrc[it] * atm.fweigh[it]);
        }
      }
      efrc += square11(cnf.fitengy - cnf.engy);
      omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
      ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
    }

    epsh *= dparams["pweight"];
    reduce(cmm, epsh, error["punish"], std::plus<double>(), PFROOT);
    reduce(cmm, efrc, error["frc"], std::plus<double>(), PFROOT);
    if (cmm.rank() == PFROOT) break;
  }
  return error["frc"] * (1. + error["punish"]);
}

double pfHome::forceMEAM(const arma::mat &vv) {
  error["frc"] = 0.0, error["punish"] = 0.0, error["shift"] = 0.0;
  omaxrho = -1e10, ominrho = 1e10;

  int cnt = 0;
  for (int i = 0; i < nfuncs; i++) {
    Func &ff = funcs[i];
    if (i == PHI || i == RHO || i == MEAMF) {
      for (int j = 0; j < ff.npts - 1; j++) ff.yy[j] = vv[cnt++];
    } else {
      for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];
    }
    ff.s.set_points(ff.xx, ff.yy);
  }

  int ls[] = {PHI, RHO, MEAMF};
  for (int it : ls) {
    for (double ee : funcs[it].s.m_b) error["punish"] += square11(ee);
    for (int i = 0; i < funcs[it].s.m_b.size() - 1; i++)
      error["punish"] += 0.5 * funcs[it].s.m_b[i] * funcs[it].s.m_b[i + 1];
  }

  error["frc"] = 0.0;
  for (Config &cnf : configs) {
    forceMEAM(cnf);
    for (pfAtom &atm : cnf.atoms)
      // if (atm.nneighsFull == N5)
      for (int it : {X, Y, Z}) {
        atm.fitfrc[it] =
            atm.phifrc[it] + atm.rhofrc[it] + atm.trifrc[it] - atm.frc[it];
        error["frc"] += square11(atm.fitfrc[it] * atm.fweigh[it]);
      }
    ominrho = cnf.rhomi < ominrho ? cnf.rhomi : ominrho;
    omaxrho = cnf.rhomx > omaxrho ? cnf.rhomx : omaxrho;
    error["frc"] += square11(cnf.fitengy - cnf.engy);
  }  // cc
  error["punish"] *= dparams["pweight"];
  error["frc"] *= 1e2;
  return error["frc"] + error["punish"];  // + error["shift"];
}

double pfHome::forceMEAM(const vector<double> &vv) {
  int cnt = 0;
  for (Func &ff : funcs)  // left examples of how to use it
    for (int j = 0; j < ff.npts; j++) ff.yy[j] = vv[cnt++];
  return error["frc"];
}

void pfHome::forceMEAM(Config &cnf) {
  cnf.phiengy = cnf.emfengy = 0.0;
  for (pfAtom &atm : cnf.atoms) { /* loop over atoms to reset values */
    atm.crho = 0.0;
    for (int it : {X, Y, Z})
      atm.phifrc[it] = atm.rhofrc[it] = atm.trifrc[it] = 0.0;
  }  // ii
  double e0 = funcs[EMF].s(0.0);
  for (pfAtom &atm : cnf.atoms) { /* loop over atoms pairs, densities */
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];

      // pair (phi)
      funcs[PHI].s.deriv(ngbj.slots[PHI], ngbj.shifts[PHI], ngbj.phi,
                         ngbj.phig);
      cnf.phiengy += ngbj.phi;
      for (int it : {X, Y, Z}) atm.phifrc[it] += ngbj.dist2r[it] * ngbj.phig;

      // rho
      if (ngbj.r < funcs[RHO].xx.back()) {
        funcs[RHO].s.deriv(ngbj.slots[RHO], ngbj.shifts[RHO], ngbj.rho,
                           ngbj.rhog);
        atm.crho += ngbj.rho;
        // partial sum
        double psum = 0.0;
        funcs[MEAMF].s.deriv(ngbj.slots[2], ngbj.shifts[2], ngbj.fval,
                             ngbj.fgrad);
        for (int kk = 0; kk < jj; kk++) {
          Neigh &ngbk = atm.neighsFull[kk];
          if (ngbk.r < funcs[RHO].xx.back()) {
            Angle &agl = atm.angMat[jj][kk];
            funcs[MEAMG].s.deriv(agl.slot, agl.shift, agl.gval, agl.ggrad);
            psum += ngbk.fval * agl.gval;
          }  // rhocut
        }    // kk
        atm.crho += psum * ngbj.fval;
      }  // r < rcut(rho)
    }    // jj

    double embE;
    funcs[EMF].s.deriv(atm.crho, embE, atm.gradF);

    cnf.emfengy += (embE - e0);

    cnf.rhomx = atm.crho > cnf.rhomx ? atm.crho : cnf.rhomx;
    cnf.rhomi = atm.crho < cnf.rhomi ? atm.crho : cnf.rhomi;

    double forces_i[3] = {0.0, 0.0, 0.0};
    for (int jj = 0; jj < atm.nneighsFull; jj++) {
      Neigh &ngbj = atm.neighsFull[jj];
      if (ngbj.r < funcs[RHO].xx.back()) {
        double f_rij = ngbj.fval;
        double f_rij_prime = ngbj.fgrad;

        double forces_j[3] = {0.0, 0.0, 0.0};
        for (int kk = 0; kk < jj; kk++) {
          Neigh &ngbk = atm.neighsFull[kk];
          if (ngbk.r < funcs[RHO].xx.back()) {
            double gcos = atm.angMat[jj][kk].gcos;
            double gval = atm.angMat[jj][kk].gval;
            double gprime = atm.angMat[jj][kk].ggrad;

            double f_rik = ngbk.fval;
            double f_rik_prime = ngbk.fgrad;

            double fij = -atm.gradF * gval * f_rik * f_rij_prime;
            double fik = -atm.gradF * gval * f_rij * f_rik_prime;

            double prefactor = atm.gradF * f_rij * f_rik * gprime;

            double prefactor_ij = prefactor / ngbj.r;
            double prefactor_ik = prefactor / ngbk.r;

            fij += prefactor_ij * gcos;
            fik += prefactor_ik * gcos;

            double fj[3], fk[3];

            for (int it : {X, Y, Z}) {
              fj[it] = ngbj.dist2r[it] * fij - ngbk.dist2r[it] * prefactor_ij;
              forces_j[it] += fj[it];

              fk[it] = ngbk.dist2r[it] * fik - ngbj.dist2r[it] * prefactor_ik;
              forces_i[it] -= fk[it];
              cnf.atoms[ngbk.aid].trifrc[it] += fk[it];
            }
          }  // r < cutoff(rho)
        }    // loop over kk
        for (int it : {X, Y, Z}) {
          atm.trifrc[it] -= forces_j[it];
          cnf.atoms[ngbj.aid].trifrc[it] += forces_j[it];
        }
      }  // r < cutoff(rho)
    }    // loop over jj
    for (int it : {X, Y, Z}) atm.trifrc[it] += forces_i[it];
  }                             // atm
  for (pfAtom &atm : cnf.atoms) /* eambedding forces */
    for (Neigh &ngb : atm.neighsFull) {
      if (ngb.r < funcs[RHO].xx.back()) {
        double emf = ngb.rhog * (atm.gradF + cnf.atoms[ngb.aid].gradF);
        for (int it : {X, Y, Z}) atm.rhofrc[it] += ngb.dist2r[it] * emf;
      }  // cutoff(rho)
    }    // nn
  cnf.fitengy = (cnf.phiengy / 2. + cnf.emfengy) / cnf.natoms;
}
