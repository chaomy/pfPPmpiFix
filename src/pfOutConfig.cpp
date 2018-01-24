/*
 * @Author: chaomy
 * @Date:   2017-12-16 16:21:26
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-16 14:38:49
 */

#include "pfHome.h"

using std::string;

void pfHome::writePOSCAR(const Config &cc, string fnm) {
  FILE *fid = fopen(fnm.c_str(), "w");
  fprintf(fid, "%s\n", sparams["elem"].c_str());
  fprintf(fid, "1 \n");
  fprintf(fid, "%f %f %f\n", cc.bvx[X], cc.bvx[Y], cc.bvx[Z]);
  fprintf(fid, "%f %f %f\n", cc.bvy[X], cc.bvy[Y], cc.bvy[Z]);
  fprintf(fid, "%f %f %f\n", cc.bvz[X], cc.bvz[Y], cc.bvz[Z]);
  for (auto ee : cc.nelemsv) fprintf(fid, "%s ", ee.c_str());
  fprintf(fid, "\n");
  for (auto ee : cc.natomsv) fprintf(fid, "%d ", ee);
  fprintf(fid, "\n");
  fprintf(fid, "Cartesian\n");
  for (const pfAtom &atm : cc.atoms)
    fprintf(fid, "%f %f %f\n", atm.pst[X], atm.pst[Y], atm.pst[Z]);
  fclose(fid);
}
