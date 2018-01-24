/*
 * @Author: chaomy
 * @Date:   2017-11-10 14:44:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-23 06:38:46
 */

#include "pfLmpDrv.h"

using std::string;
using std::vector;

void pfLMPdrv::calItenRun(const string& tag, double& engy, vector<double>& stsv) {
  char lmpcmds[100][MAXLEN];
  int i = 0;

  //  --------------------- INITIALIZAITION ---------------------
  sprintf(lmpcmds[i++], "clear");
  sprintf(lmpcmds[i++], "units  metal");
  sprintf(lmpcmds[i++], "dimension  3");
  sprintf(lmpcmds[i++], "boundary p p p");

  if (!tag.compare("tpath")) {
    sprintf(lmpcmds[i++],
            "lattice custom %f "
            "a1 %8.6f %8.6f %8.6f a2 %8.6f %8.6f %8.6f a3 %8.6f %8.6f "
            "%8.6f "
            "basis 0 0 0 basis 0.5 0.5 0.5",
            exprs["lat"], bsm[0][0], bsm[0][1], bsm[0][2], bsm[1][0], bsm[1][1],
            bsm[1][2], bsm[2][0], bsm[2][1], bsm[2][2]);
  } else if (!tag.compare("opath")) {
    sprintf(lmpcmds[i++],
            "lattice custom %f "
            "a1 %8.6f %8.6f %8.6f a2 %8.6f %8.6f %8.6f a3 %8.6f %8.6f "
            "%8.6f "
            "basis 0.0 0.0 0.0 basis 0.5 0.5 0.0 "
            "basis 0.5 0.0 0.5 basis 0.0 0.5 0.5 ",
            exprs["lat"], bsm[0][0], bsm[0][1], bsm[0][2], bsm[1][0], bsm[1][1],
            bsm[1][2], bsm[2][0], bsm[2][1], bsm[2][2]);
  }

  /*** use units lattice ! ***/
  sprintf(lmpcmds[i++],
          "region    box   block  0  1  0   1   0   1  units lattice");
  sprintf(lmpcmds[i++], "create_box    1   box");
  sprintf(lmpcmds[i++], "create_atoms  1   region   box");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(lmpcmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(lmpcmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(lmpcmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(lmpcmds[i++], "neighbor 1.0 bin");
  sprintf(lmpcmds[i++], "neigh_modify  every 1  delay  0 check yes");

  sprintf(lmpcmds[i++], "compute  mystress  all stress/atom  NULL");

  sprintf(lmpcmds[i++], "compute p1 all reduce sum c_mystress[1]");
  sprintf(lmpcmds[i++], "compute p2 all reduce sum c_mystress[2]");
  sprintf(lmpcmds[i++], "compute p3 all reduce sum c_mystress[3]");
  sprintf(lmpcmds[i++], "compute p4 all reduce sum c_mystress[4]");
  sprintf(lmpcmds[i++], "compute p5 all reduce sum c_mystress[5]");
  sprintf(lmpcmds[i++], "compute p6 all reduce sum c_mystress[6]");

  /*** stress unts is changed from bars to Mpa x 0.1 to Gpa x 1e-4 ***/
  sprintf(lmpcmds[i++], "variable  Sxx equal 1e-4*(c_p1)/vol");  // xx
  sprintf(lmpcmds[i++], "variable  Syy equal 1e-4*(c_p2)/vol");  // yy
  sprintf(lmpcmds[i++], "variable  Szz equal 1e-4*(c_p3)/vol");  // zz
  sprintf(lmpcmds[i++], "variable  Sxy equal 1e-4*(c_p4)/vol");  // zz
  sprintf(lmpcmds[i++], "variable  Sxz equal 1e-4*(c_p5)/vol");  // zz
  sprintf(lmpcmds[i++], "variable  Syz equal 1e-4*(c_p6)/vol");  // zz

  // ---------------------- THERMO  --------------------------
  sprintf(lmpcmds[i++], "thermo  10000");
  sprintf(lmpcmds[i++],
          "thermo_style custom  step  pe  etotal  lx   ly   lz   v_Sxx  "
          "v_Syy  v_Szz  v_Sxy  v_Sxz  v_Syz");

  sprintf(lmpcmds[i++], "min_style  cg");
  sprintf(lmpcmds[i++], "minimize  1e-12  1e-12  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, lmpcmds[iter]);

  engy = lammps_get_thermo(lmp, ETOTAL);

  stsv[0] = *(double*)lammps_extract_variable(lmp, SXX, NULL);
  stsv[1] = *(double*)lammps_extract_variable(lmp, SYY, NULL);
  stsv[2] = *(double*)lammps_extract_variable(lmp, SZZ, NULL);
  stsv[3] = *(double*)lammps_extract_variable(lmp, SXY, NULL);
  stsv[4] = *(double*)lammps_extract_variable(lmp, SXZ, NULL);
  stsv[5] = *(double*)lammps_extract_variable(lmp, SYZ, NULL);
}