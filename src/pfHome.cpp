/*
 * @Author: chaomy
 * @Date:   2018-01-15 00:24:43
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-24 00:05:21
 */

#include "pfHome.h"
#include "pfLmpDrv.h"
#include "pfOptimizer.h"
namespace mpi = boost::mpi;

pfHome::pfHome(int argc, char* argv[])
    : ricut(2.08),
      rocut(5.05),
      mfrc(3),
      hil({1.5, 0.5, 0.0, 1.0, 1.0}),
      lol({-1.0, -0.1, -10.0, -0.5, -0.5}),
      gradRight(5, 0) {
  printf("I am %d of %d \n", cmm.rank(), cmm.size());

  calfrc["MEAM"] = &pfHome::forceMEAM;
  calfrc["EAM"] = &pfHome::forceEAM;
  calobj["MEAM"] = &pfHome::forceMEAM;
  calobj["EAM"] = &pfHome::forceEAM;

  write["MEAM"] = &pfHome::writeMEAM;
  write["EAM"] = &pfHome::writeLMPS;
  write["TMP"] = &pfHome::writePot;

  if (cmm.rank() == PFROOT) {
    ftn = tln = 0;
    for (int ii : {0, 1, 2}) mfrc[ii] = 0.0;
    iparams["atomicNum"] = 40;
    dparams["mass"] = 92.906400;
    sparams["tmpdir"] = string("dirtmp");
    sparams["lmpdir"] = string("dirlmp");
    outMkdir(sparams["tmpdir"]);
    outMkdir(sparams["lmpdir"]);

    sparams["tmpfile"] = string("dummy.tmp");
    sparams["parfile"] = string("dummy.param");
    sparams["cnffile"] = string("dummy.config");
    sparams["potfile"] = string("dummy.pot");
    sparams["lmppot"] = string("dummy.lmp");

    lorho = 0.4, hirho = 1.0;
    gradRight[EMF] = 1, gradRight[MEAMG] = 1;
    parseArgs(argc, argv);
    pfInit();
  }

  (cmm.barrier)();  //  important!
  bcdata();
  if (!sparams["ptype"].compare("MEAM")) {
    initNeighsFull();
    initAngles();
  } else
    initNeighs();
  // forceDis();

  int del = nconfs / cmm.size();
  locstt = del * cmm.rank();
  locend = del * (cmm.rank() + 1);
  cout << "report " << del << " " << locstt << " " << locend << endl;
  (cmm.barrier)();  //  important!
  // temporarily close these functionalities
  // optdrv = new pfOptimizer(this);
  // lmpdrv = new pfLMPdrv(argc, argv, this);
};

pfHome::~pfHome() {
  funcs.clear();
  ini.clear();
  lob.clear();
  hib.clear();
  configs.clear();
  startps.clear();
  if (lmpdrv) delete lmpdrv;
  if (optdrv) delete optdrv;
}

void pfHome::pfInit() {
  initParam();
  initTargs();
  readConfig();
  readPot();
}
