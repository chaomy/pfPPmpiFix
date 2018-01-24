/*
 * @Author: chaomy
 * @Date:   2017-11-23 09:40:32
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-23 11:56:37
 */

#include "pfLmpDrv.h"

void pfLMPdrv::calVac() {
  char lmpcmds[100][MAXLEN];
  int i = 0;

  // bulk
  sprintf(lmpcmds[i++], "clear");
  sprintf(lmpcmds[i++], "units  metal");
  sprintf(lmpcmds[i++], "dimension  3");
  sprintf(lmpcmds[i++], "boundary p p p");
  sprintf(lmpcmds[i++], "atom_style atomic");
  sprintf(lmpcmds[i++], "variable  a equal  %.7f", targs["lat"]);
  sprintf(lmpcmds[i++], "variable  sz equal 10*${a}");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(lmpcmds[i++], "lattice   bcc  ${a}");
  sprintf(lmpcmds[i++],
          "region whole block 0 ${sz} 0 ${sz} 0 ${sz}  units box");
  sprintf(lmpcmds[i++], "create_box  1  whole");
  sprintf(lmpcmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 0  orient  z 0 0 1");
  sprintf(lmpcmds[i++], "create_atoms 1 region whole");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(lmpcmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(lmpcmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(lmpcmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(lmpcmds[i++], "neighbor 1.0 bin");
  sprintf(lmpcmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(lmpcmds[i++], "thermo 10000");
  sprintf(lmpcmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(lmpcmds[i++],
          "fix  1  all box/relax  x  0.0  y  0.0  z  0.0  couple none  vmax "
          "0.001");
  sprintf(lmpcmds[i++], "min_style  cg");
  sprintf(lmpcmds[i++], "minimize  1e-12  1e-12  100000  100000");

  for (int mm = 0; mm < i; mm++) lammps_command(lmp, lmpcmds[mm]);

  double vac =
      lammps_get_thermo(lmp, ETOTAL) / double(lammps_get_thermo(lmp, NATOM));
  lammps_command(lmp, CLEAR);

  i = 0.0;
  // vac
  sprintf(lmpcmds[i++], "units  metal");
  sprintf(lmpcmds[i++], "dimension  3");
  sprintf(lmpcmds[i++], "boundary p p p");
  sprintf(lmpcmds[i++], "atom_style atomic");
  sprintf(lmpcmds[i++], "variable  a equal  %.7f", targs["lat"]);
  sprintf(lmpcmds[i++], "variable  sz equal 10*${a}");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(lmpcmds[i++], "lattice   bcc  ${a}");
  sprintf(lmpcmds[i++],
          "region whole block 0 ${sz} 0 ${sz} 0 ${sz}  units box");
  sprintf(lmpcmds[i++], "create_box  1  whole");
  sprintf(lmpcmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 0  orient  z 0 0 1");
  sprintf(lmpcmds[i++], "create_atoms 1 region whole");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(lmpcmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(lmpcmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(lmpcmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(lmpcmds[i++], "neighbor 1.0 bin");
  sprintf(lmpcmds[i++], "neigh_modify  every 1  delay  0 check yes");

  sprintf(lmpcmds[i++], "group del id 1000");
  sprintf(lmpcmds[i++], "delete_atoms group del");

  // ---------------------- THERMO  --------------------------
  sprintf(lmpcmds[i++], "thermo 10000");
  sprintf(lmpcmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(lmpcmds[i++],
          "fix  1  all box/relax  x  0.0  y  0.0  z  0.0  couple none  vmax "
          "0.001");
  sprintf(lmpcmds[i++], "min_style  cg");
  sprintf(lmpcmds[i++], "minimize  1e-12  1e-12  100000  100000");

  for (int mm = 0; mm < i; mm++) lammps_command(lmp, lmpcmds[mm]);

  vac = lammps_get_thermo(lmp, ETOTAL) - vac * lammps_get_thermo(lmp, NATOM);
  // vac /= double(lammps_get_thermo(lmp, NATOM));
  lammps_command(lmp, CLEAR);

  exprs["vac"] = vac;
}