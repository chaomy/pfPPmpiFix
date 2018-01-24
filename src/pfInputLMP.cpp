/*
 * @Author: chaomy
 * @Date:   2017-12-17 14:00:51
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-18 11:58:07
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::vector;

void pfHome::readLmpMEAM() {
  funcs.clear();
  ini.clear();

  ifstream fid;
  fid.open(sparams["lmppot"].c_str());
  if (!fid.is_open()) cerr << "error open " << sparams["lmppot"] << endl;

  string buff;
  vector<string> segs(1, " ");

  getline(fid, buff);  // read head line
  int cnt = nfuncs = 5;
  while (--cnt >= 0) {
    Func tm;
    getline(fid, buff);
    tm.npts = stoi(buff);
    cout << tm.npts << endl;
    tm.g1 = vector<double>(tm.npts, 0);
    tm.g2 = vector<double>(tm.npts, 0);
    tm.xx = vector<double>(tm.npts, 0);
    tm.yy = vector<double>(tm.npts, 0);
    getline(fid, buff);
    sscanf(buff.c_str(), "%lf %lf", &tm.g1.front(), &tm.g1.back());
    getline(fid, buff);
    for (int j = 0; j < tm.npts; j++) {
      getline(fid, buff);
      sscanf(buff.c_str(), "%lf %lf %lf", &tm.xx[j], &tm.yy[j], &tm.g2[j]);
    }
    funcs.push_back(tm);
  }
  fid.close();
  ricut = funcs[PHI].xx.front();
  rocut = funcs[PHI].xx.back();
  rhcut = funcs[RHO].xx.back();

  // for (auto ff : funcs) {
  //   cout << "npts = " << ff.npts << endl;
  //   for (int i = 0; i < ff.npts; i++)
  //     cout << ff.xx[i] << " " << ff.yy[i] << " " << ff.g2[i] << endl;
  // }
}