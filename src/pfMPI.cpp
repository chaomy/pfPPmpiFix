/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:11:45
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-23 22:55:51
 */

#include "pfHome.h"

using std::cout;
using std::endl;

void pfHome::bcdata() {
  /* params */
  broadcast(cmm, dparams, PFROOT);
  broadcast(cmm, iparams, PFROOT);
  broadcast(cmm, sparams, PFROOT);

  /* funcstions */
  broadcast(cmm, nfuncs, PFROOT);
  broadcast(cmm, ini, PFROOT);
  broadcast(cmm, lob, PFROOT);
  broadcast(cmm, hib, PFROOT);
  broadcast(cmm, deb, PFROOT);
  broadcast(cmm, startps, PFROOT);
  broadcast(cmm, endps, PFROOT);

  if (cmm.rank() != PFROOT) {
    for (int i = 0; i < nfuncs; i++) {
      Func tmp;
      funcs.push_back(tmp);
    }
  }

  broadcast(cmm, nvars, PFROOT);
  for (int i = 0; i < nfuncs; i++) {
    broadcast(cmm, funcs[i].npts, PFROOT);
    broadcast(cmm, funcs[i].xx, PFROOT);
    broadcast(cmm, funcs[i].yy, PFROOT);
    broadcast(cmm, funcs[i].g1, PFROOT);
    broadcast(cmm, funcs[i].g2, PFROOT);
  }

  // assign cutoffs
  ricut = funcs[PHI].xx.front();
  rocut = funcs[PHI].xx.back();
  rhcut = funcs[RHO].xx.back();

  if (!sparams["ptype"].compare("MEAM")) {  // boundary
    funcs[MEAMF].s.set_boundary(tk::spline::second_deriv, 0.0,
                                tk::spline::first_deriv, 0.0, true);
    funcs[MEAMG].s.set_boundary(tk::spline::second_deriv, 0.0,
                                tk::spline::second_deriv, 0.0, true);
  }
  funcs[PHI].s.set_boundary(tk::spline::second_deriv, 0.0,
                            tk::spline::first_deriv, 0.0, true);
  funcs[RHO].s.set_boundary(tk::spline::second_deriv, 0.0,
                            tk::spline::first_deriv, 0.0, true);
  funcs[EMF].s.set_boundary(tk::spline::second_deriv, 0.0,
                            tk::spline::second_deriv, 0.0, true);

  /* configurations */
  broadcast(cmm, nconfs, PFROOT);
  broadcast(cmm, configs, PFROOT);
}

// void pfHome::bfunc() {}