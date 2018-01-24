/*
 * @Author: chaomy
 * @Date:   2017-11-10 14:28:37
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-25 11:25:11
 */

#include "pfLmpDrv.h"

void pfLMPdrv::calLatticeBCC() {
  char cmds[100][MAXLEN];
  int i = 0;

  //  --------------------- INITIALIZAITION ---------------------
  sprintf(cmds[i++], "clear");
  // sprintf(cmds[i++], "print logfile screen no");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "atom_style atomic");
  sprintf(cmds[i++], "variable  a equal  %.7f", targs["lat"]);

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++], "region whole block 0 ${a} 0 ${a} 0 ${a}  units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 0  orient  z 0 0 1");
  sprintf(cmds[i++], "create_atoms 1 region whole");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(cmds[i++], "thermo 10000");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(cmds[i++],
          "fix  1  all box/relax  x  0.0  y  0.0  z  0.0  couple none  vmax "
          "0.001");
  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-12  1e-12  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  /* extract the lattice constant */
  exprs["lat"] = exprs["abcc"] = lammps_get_thermo(lmp, LX);
  exprs["ebcc"] =
      lammps_get_thermo(lmp, ETOTAL) / double(lammps_get_thermo(lmp, NATOM));

  lammps_command(lmp, CLEAR);

  error["lat"] = abs(exprs["abcc"] - targs["abcc"]) / targs["abcc"];
}