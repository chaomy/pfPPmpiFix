/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-24 11:04:41
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

  // func -> ini
  nvars = 0;
  for (Func& ff : funcs) ff.step = ff.xx[1] - ff.xx[0];
  for (int i = 0; i < nfuncs; i++) {
    Func& ff = funcs[i];
    startps.push_back(nvars);
    int nt = (i == PHI || i == RHO || i == MEAMF) ? ff.npts - 1 : ff.npts;
    for (int j = 0; j < nt; j++) {
      ini.push_back(ff.yy[j]);
      lob.push_back(lol[i]);
      hib.push_back(hil[i]);
      deb.push_back(hil[i] - lol[i]);
      nvars++;
    }
    endps.push_back(nvars);
  }  // i

  if (!sparams["ptype"].compare("MEAM"))  // boundary
    funcs[MEAMF].yy.back() = 0.0;
  funcs[PHI].yy.back() = 0.0;
  funcs[RHO].yy.back() = 0.0;
}