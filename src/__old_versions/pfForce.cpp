/*
 * @Author: chaomy
 * @Date:   2017-10-30 15:31:59
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-26 15:53:37
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
#include "pfOptimizer.h"

using std::cout;
using std::endl;

void pfHome::increAnneal() {
  scnt = 0;
  // int jobl[] = {PHI, PHI, MEAMF, MEAMG,
  //               PHI, PHI, MEAMF, MEAMG, EMF, PHI};
  int jobl[] = {PHI, PHI, RHO, MEAMF, MEAMG, PHI, RHO, PHI, RHO, MEAMF, MEAMG};
  FILE *ptfile;
  for (int i = 0; i < 10; i++) {
    simAnneal();
    upgrade(jobl[i]);  // enlarge the PHI
    iparams["kmax"]++;
    ptfile = fopen("record.txt", "a");
    fprintf(ptfile, "%d %f\n", i, recorderr[i]);
    fclose(ptfile);
  }
}

void pfHome::run(int argc, char *argv[]) {
  if (!sparams["opt"].compare("lmp"))
    lmpdrv->calPhy();
  else if (!sparams["opt"].compare("make") || !sparams["opt"].compare("err"))
    calErr();
  else if (!sparams["opt"].compare("anneal")) {
    simAnneal();
    nloptGlobal();
  } else if (!sparams["opt"].compare("incre"))
    increAnneal();
  else if (!sparams["opt"].compare("nlopt"))
    nloptGlobal();
  else if (!sparams["opt"].compare("shift"))
    doShift();
}

void pfHome::calErr() {
  // make potential
  double err = 0.0;
  if (!sparams["ptype"].compare("EAM"))
    err = forceEAM(ini, 0);
  else if (!sparams["ptype"].compare("ADP"))
    err = forceADP(ini, 0);
  else if (!sparams["ptype"].compare("MEAM"))
    err = forceMEAM(ini, 0);
  printf("force err = %.10f\n", err);
  writeLMPS(ini);
  // run lammps
  lmpdrv->calPhy();
}

void pfHome::nloptGlobal() {
  gcnt = 0;
  optdrv->prepOpt();
  optdrv->runOpt();
  writePot(ini);
}

void pfHome::doShift() {
  cout << "before err = " << forceMEAM(ini, 0) << endl;

  sparams["tmpfile"] = "dummy.tmp.0";
  writePot(ini);

  updaterhoMEAM(ini);
  double rho1 = oaverho;
  shiftRHO(ini);

  updaterhoMEAM(ini);
  double rho2 = oaverho;

  shiftEMF(rho2 - rho1);

  sparams["tmpfile"] = "dummy.tmp.1";
  writePot(ini);

  cout << "after err = " << forceMEAM(ini, 0) << endl;
}

double pfHome::errFunct(const vector<double> &x) {
  double err = 0.0;
  int pfreq = 500;
  if (sparams["ptype"] == "EAM") {
    gcnt++;
    err = forceEAM(x, 0);
    if (gcnt % pfreq == 0) printf("cnt %d err = %f \n", gcnt++, err);
  }
  return err;
}

double pfHome::errFunctGrad(const vector<double> &x, vector<double> &grad) {
  double y1, y2;
  vector<double> tmpx = x;
  for (unsigned int i = 0; i < x.size(); i++) {
    double dx = 1e-4 * x[i];
    tmpx[i] += dx;
    y1 = errFunct(tmpx);
    tmpx[i] -= (2 * dx);
    y2 = errFunct(tmpx);
    grad[i] = (y1 - y2) / (2 * dx);
  }
  return errFunct(x);
}