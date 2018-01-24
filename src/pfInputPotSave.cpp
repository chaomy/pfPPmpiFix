/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-13 16:35:33
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

// read dummy.pot
void pfHome::readPot() {
  funcs.clear();
  ini.clear();

  ifstream fid;
  pfUtil pfu;

  char tmp[MAXLEN];
  fid.open(sparams["potfile"].c_str());
  if (!fid.is_open()) cerr << "error opening " + sparams["potfile"] << endl;

  string buff;
  vector<string> segs(1, " ");
  vector<int> bnds;
  vector<int> npts;

  while (getline(fid, buff)) {
    segs.clear();
    pfu.split(buff, " ", segs);
    if (!segs[0].compare("#F"))
      sscanf(buff.c_str(), "%s %s %d", tmp, tmp, &nfuncs);
    else if (!segs[0].compare("#T"))
      sparams["ptype"] = segs[1];
    else if (!segs[0].compare("#G"))
      for (unsigned int i = 1; i < segs.size(); i++)
        bnds.push_back(stoi(segs[i]));
    else if (!segs[0].compare("#E"))
      break;
  }

  for (int i = 0; i < nfuncs; i++) {
    getline(fid, buff);
    npts.push_back(stoi(buff));
  }

  double v2[2];
  double b2[2];

  for (int i = 0; i < nfuncs; i++) {
    Func tmp;
    tmp.bnd = bnds[i];
    tmp.npts = npts[i];

    getline(fid, buff);
    getline(fid, buff);
    sscanf(buff.c_str(), "%lf %lf", &b2[0], &b2[1]);

    for (int j = 0; j < npts[i]; j++) {
      getline(fid, buff);
      sscanf(buff.c_str(), "%lf %lf", &v2[0], &v2[1]);
      tmp.xx.push_back(v2[0]);
      tmp.yy.push_back(v2[1]);
      tmp.g1.push_back(0.0);
      tmp.g2.push_back(0.0);
    }
    for (int it = 0; it < 2; it++) {
      if (b2[it] < lol.back()) b2[it] = 0.5 * lol.back();
      if (b2[it] > hil.back()) b2[it] = 0.5 * hil.back();
    }
    tmp.g1.front() = b2[0];
    tmp.g1.back() = b2[1];
    funcs.push_back(tmp);
  }
  fid.close();

  // assign cutoffs
  ricut = funcs[PHI].xx.front();
  rocut = funcs[PHI].xx.back();
  rhcut = funcs[RHO].xx.back();

  // func -> ini
  nvars = 0;
  for (Func& ff : funcs) ff.step = ff.xx[1] - ff.xx[0];
  for (int i = 0; i < nfuncs; i++) {
    Func& ff = funcs[i];
    startps.push_back(nvars);
    endps.push_back(nvars + ff.npts);
    for (int j = 0; j < ff.npts; j++) {
      ini.push_back(ff.yy[j]);
      lob.push_back(lol[i]);
      hib.push_back(hil[i]);
      nvars++;
    }
  }  // i

  if (!sparams["ptype"].compare("MEAM")) {  // boundary
    ini.push_back(funcs[PHI].g1.front());
    ini.push_back(funcs[EMF].g1.front());
    ini.push_back(funcs[MEAMG].g1.front());
    ini.push_back(funcs[EMF].g1.back());
    ini.push_back(funcs[MEAMG].g1.back());
    nvars += 5;
    for (int i : {0, 1, 3}) funcs[i].g1.back() = 0.0;
    for (int i : {1, 3}) funcs[i].g1.front() = 0.0;
  } else if (!sparams["ptype"].compare("EAM")) {
    hirho = funcs[EMF].xx.back();
    lorho = funcs[EMF].xx.front();
    if (sparams["spline"] != "nat") {
      funcs[PHI].g1.front() = 1e30;
      funcs[RHO].g1.front() = 1e30;
      funcs[EMF].g1.front() = 1e30;
      funcs[EMF].g1.back() = 1e30;
    } else {  // optimize  spline
      ini.push_back(funcs[PHI].g1.front());
      ini.push_back(funcs[RHO].g1.front());
      ini.push_back(funcs[EMF].g1.front());
      ini.push_back(funcs[EMF].g1.back());
      nvars += 4;
    }
    funcs[PHI].g1.back() = 0.0;
    funcs[RHO].g1.back() = 0.0;
  }
  funcs[PHI].s.set_boundary(tk::spline::second_deriv, 0.0,
                            tk::spline::first_deriv, 0.0, false);
  funcs[RHO].s.set_boundary(tk::spline::second_deriv, 0.0,
                            tk::spline::first_deriv, 0.0, false);
  funcs[EMF].s.set_boundary(tk::spline::second_deriv, 0.0,
                            tk::spline::second_deriv, 0.0, false);
}