/*
 * @Author: chaomy
 * @Date:   2017-12-19 08:58:15
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-19 15:57:43
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::pow;
using std::string;
using std::to_string;
using std::vector;

void pfHome::calPV() {
  string ptg(sparams["elem"] + "p");
  string vtg(sparams["elem"] + "v");
  string etg(sparams["elem"] + "e");

  const vector<double>& pp = mele.pvm[ptg];
  const vector<double>& vv = mele.pvm[vtg];

  mpcf["pv"].clear();
  for (double dl : vv)
    mpcf["pv"].push_back(addvolm(ubcc, pow(dl, 1. / 3.) - 1.0));

  error["pv"] = 0.0;
  for (int i = 0; i < mpcf["pv"].size(); i++) {
    auto& ee = mpcf["pv"][i];
    stressMEAM(ee);
    for (int i = 0; i < 6; i++) ee.strs[i] *= EVA3_GPA;
    error["pv"] += square11(pp[i] + ee.strs[XX]);
  }
}