/*
 * @Author: chaomy
 * @Date:   2017-11-05 22:29:46
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-23 15:09:36
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::string;
using std::unordered_map;

void pfHome::initParam() {
  sparams["elem"] = string("Nb");
  sparams["ptype"] = string("EAM");
  sparams["pairstyle"] = string("meam/spline");
  sparams["alg"] = string("LN_SBPLX");
  sparams["opt"] = string("nlopt");

  dparams["temp"] = 0.01;
  dparams["istep"] = 0.1;
  dparams["pratio"] = 0.5;
  dparams["eweight"] = 1.0;
  dparams["xtol"] = 1e-4;
  dparams["bwidth"] = 15.0;

  iparams["maxstep"] = 10000;
  iparams["resfreq"] = 10;
  iparams["lmpfreq"] = 15;
  iparams["kmax"] = 8;  // number of outer loop in simulated annealing
  readParam();
  sparams["lmpfile"] = string("dummy.lammps.") + sparams["ptype"];
}

void pfHome::readParam() {
  ifstream fid;
  pfUtil pfu;

  fid.open(sparams["parfile"].c_str());
  if (!fid.is_open()) cerr << " error opening " + sparams["parfile"] << endl;
  vector<string> segs(1, " ");
  string buff;

  while (getline(fid, buff)) {
    segs.clear();
    pfu.split(buff, " ", segs);
    cout << segs[0] <<" "<< segs[1] << endl;
    if (!segs[0].compare("alg"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("elem"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("opt"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("ptype"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("pairstyle"))
      sparams[segs[0]] = segs[1];
    else if (!segs[0].compare("maxstep"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("resfreq"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("lmpfreq"))
      iparams[segs[0]] = stoi(segs[1]);
    else if (!segs[0].compare("kmax"))
      iparams[segs[0]] = stoi(segs[1]);
    else
      dparams[segs[0]] = stof(segs[1]);
  }
  fid.close();
}

void pfHome::parseArgs(int argc, char *argv[]) {
  /**** parse the command line ****/
  for (int i = 0; i < argc; i++) {
    if (!strcmp(argv[i], "--p") || !strcmp(argv[i], "-p"))
      sparams["parfile"] = string(argv[++i]);
    else if (!strcmp(argv[i], "--c") || !strcmp(argv[i], "-c"))
      sparams["cnffile"] = string(argv[++i]);
    else if (!strcmp(argv[i], "--f") || !strcmp(argv[i], "-f"))
      sparams["potfile"] = string(argv[++i]); 
  }
  return;
}