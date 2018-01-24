/*
 * @Author: chaomy
 * @Date:   2018-01-20 16:53:38
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-23 19:26:41
 */

#include "pfHome.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::vector;

void pfHome::readConfig() { /* read atomic force file */
  configs.clear();          // clear

  ifstream fid;
  pfUtil pfu;

  char tmp[MAXLEN];
  fid.open(sparams["cnffile"].c_str());

  if (!fid.is_open()) cerr << "error open" + sparams["cnffile"] << endl;

  int cnt = 0;
  string buff;
  vector<string> segs(1, " ");
  Config config;

  while (getline(fid, buff)) {
    segs.clear();
    if (!config.atoms.empty()) config.atoms.clear();
    pfu.split(buff, " ", segs);
    if (!segs[0].compare("#N")) {
      sscanf(buff.c_str(), "%s %d %s", tmp, &config.natoms, tmp);
    } else if (!segs[0].compare("#X")) {
      sscanf(buff.c_str(), "%s %lf %lf %lf", tmp, &config.bvx[0],
             &config.bvx[1], &config.bvx[2]);
    } else if (segs[0] == "#Y") {
      sscanf(buff.c_str(), "%s %lf %lf %lf", tmp, &config.bvy[0],
             &config.bvy[1], &config.bvy[2]);
    } else if (segs[0] == "#Z") {
      sscanf(buff.c_str(), "%s %lf %lf %lf", tmp, &config.bvz[0],
             &config.bvz[1], &config.bvz[2]);
    } else if (!segs[0].compare("#W")) {
      sscanf(buff.c_str(), "%s %lf", tmp, &config.weigh);
    } else if (!segs[0].compare("#E")) {
      sscanf(buff.c_str(), "%s %lf", tmp, &config.engy);
    } else if (!segs[0].compare("#F")) {
      for (int i = 0; i < config.natoms; i++) {
        getline(fid, buff);
        pfAtom atom(i);
        sscanf(buff.c_str(), "%s %lf %lf %lf %lf %lf %lf", tmp, &atom.pst[0],
               &atom.pst[1], &atom.pst[2], &atom.frc[0], &atom.frc[1],
               &atom.frc[2]);
        atom.absfrc = sqrt(square33(atom.frc));
        for (int ii : {0, 1, 2}) {
          atom.fweigh[ii] = 1. / (fabs(atom.frc[ii]) + FRCEPS);
          // atom.fweigh[ii] =
          //     1e4 * exp(-dparams["bwidth"] * atom.frc[ii] * atom.frc[ii]);
          // mfrc[ii] += fabs(atom.frc[ii]);
        }
        config.atoms.push_back(atom);
      }
      tln += config.natoms;
      config.cfgid = cnt++;
      configs.push_back(config);
    }  // #F
  }    // while

  for (int ii : {0, 1, 2}) mfrc[ii] /= tln;
  printf("%f %f %f\n", mfrc[0], mfrc[1], mfrc[2]);

  fid.close();
  cout << "finish reading " << configs.size() << " configs" << endl;
  nconfs = configs.size();
}