/*
 * @Author: chaomy
 * @Date:   2017-11-13 15:58:23
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-31 15:14:18
 */

#include "pfLmpDrv.h"

using namespace LAMMPS_NS;
using std::cout;
using std::endl;

pfLMPdrv::pfLMPdrv(int argc, char* argv[]) : pfhm(NULL) {
  paraInit(argc, argv);
}
pfLMPdrv::pfLMPdrv(int argc, char* argv[], pfHome* pt) : pfhm(pt) {
  paraInit(argc, argv);
}
pfLMPdrv::~pfLMPdrv() {
  if (lmp) delete lmp;
  MPI_Finalize();
}

void pfLMPdrv::calPhy() {
  if (label["bcc"]) calLatticeBCC();
  if (label["fcc"]) calLatticeFCC();
  if (label["hcp"]) calLatticeHCP();
  if (label["cij"]) calElastic();
  if (label["suf"]) calSurface();
  if (label["gsf"]) calGSF();
  if (label["pv"]) calPV();
  if (label["vac"]) calVac();
  if (label["iten"]) calIten();
  calPhyErr();
  // remove("no");
  // remove("log.lammps");
}

void pfLMPdrv::calPhyErr() {
  exprs["bcc2hcp"] = exprs["ehcp"] - exprs["ebcc"];
  exprs["bcc2fcc"] = exprs["efcc"] - exprs["ebcc"];

  vector<string> aa({"lat", "c11", "c12", "c44", "suf110", "suf100", "suf111",
                     "bcc2fcc", "bcc2hcp"});
  for (string ee : aa)
    cout << ee << " " << exprs[ee] << " " << targs[ee] << endl;
  for (string ee : aa) error[ee] = 100 * relerr(exprs[ee], targs[ee]);
  for (string ee : aa) cout << "error " << ee << " " << error[ee] << endl;

  error["gsf110"] = error["gsf211"] = 0.0;

  for (double ee : lgsf["111e110"]) error["gsf110"] += ee;
  error["gsf110"] /= lgsf["111e110"].size();

  for (double ee : lgsf["111e211"]) error["gsf211"] += ee;
  error["gsf211"] /= lgsf["111e211"].size();

  error["gsf110"] *= 100;
  error["gsf211"] *= 100;

  // usf
  exprs["usf211"] = exprs["usf110"] = 0.0;
  for (double ee : lgsf["111z110"]) exprs["usf110"] = fmax(exprs["usf110"], ee);
  for (double ee : lgsf["111z211"]) exprs["usf211"] = fmax(exprs["usf211"], ee);

  // pv
  error["pv"] = 0.0;
  for (double ee : lmpv["Nbe"]) error["pv"] += ee;
  error["pv"] /= lmpv["Nbe"].size();
  error["pv"] *= 100;

  // output
  string fmt(
      "%3d %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f %06.2f "
      "%06.2f\n");
  FILE* fid = fopen("perr.log", "a");
  fprintf(fid, fmt.c_str(), cnt++, exprs["lat"], exprs["c11"], exprs["c12"],
          exprs["c44"], exprs["suf110"], exprs["suf100"], exprs["suf111"],
          exprs["usf110"], exprs["usf211"], lmpv["Nbp"].back());
  fprintf(fid, fmt.c_str(), 0, error["lat"], error["c11"], error["c12"],
          error["c44"], error["suf110"], error["suf100"], error["suf111"],
          error["gsf110"], error["gsf211"], error["pv"]);
  fclose(fid);
}

void pfLMPdrv::lmpTargets() {
  sttag["elem"] = pfhm->gsparams()["elem"];
  sttag["pairstyle"] = pfhm->gsparams()["pairstyle"];
  sttag["lmpfile"] = pfhm->gsparams()["lmpfile"];
  error = pfhm->error;
  targs = pfhm->targs;
  weigh = pfhm->weigh;
}

void pfLMPdrv::paraInit(int argc, char* argv[]) {
  label["bcc"] = 1;
  label["fcc"] = 1;
  label["hcp"] = 1;
  label["cij"] = 1;
  label["suf"] = 1;
  label["gsf"] = 1;
  label["pv"] = 1;
  // label["vac"] = 1;
  // label["iten"] = 1;

  // mpi
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  lmpTargets();  // targets

  char** lmparg = new char*[3];
  lmparg[0] = NULL;  // placeholder for program name
  lmparg[1] = (char*)"-screen";
  lmparg[2] = (char*)"no";

  lmp = new LAMMPS(3, lmparg, MPI_COMM_WORLD);
  delete[] lmparg;
}
