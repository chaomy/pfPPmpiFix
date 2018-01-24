/*
 * @Author: yangchaoming
 * @Date:   2017-10-23 14:04:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-25 14:34:52
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

enum {LO=0, HI=1};
void pfHome::readPot() {
  funcs.clear();
  ini.clear();
  // read dummy.pot
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
    for (int it = 0; it < 2; it++){
      if (b2[it] < lol.back()) b2[it] = 0.5 * lol.back();
      if (b2[it] > hil.back()) b2[it] = 0.5 * hil.back();
    }
    tmp.g1.front() = b2[0];
    tmp.g1.back() = b2[1];
    funcs.push_back(tmp);
  }
  fid.close();

  // func -> vv  
  nvars = 0;
  // cout << "num of funcs = " << nfuncs << " " << funcs.size() << endl;
  for (int i = 0; i < nfuncs; i++) {
    startps.push_back(nvars);
    for (int j = 0; j < funcs[i].npts; j++) {
      ini.push_back(funcs[i].yy[j]);
      lob.push_back(lol[i]);
      hib.push_back(hil[i]);
      nvars++;  // total nvars;
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
  }  // i
  rocut = funcs[PHI].xx.back();
}