/*
 * @Author: chaomy
 * @Date:   2017-11-10 14:36:04
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-26 13:57:39
 */

#include "pfLmpDrv.h"

void pfLMPdrv::calSurface() {
  double eB100 = 0.0, eS100 = 0.0;
  double eB110 = 0.0, eS110 = 0.0;
  double eB111 = 0.0, eS111 = 0.0;

  double areaS100 = 0.0, areaS110 = 0.0, areaS111 = 0.0;

  int i = 0;
  char cmds[100][MAXLEN];

  /*******************************************************************
   * calculate the  100  with surface
   * *****************************************************************/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);

  sprintf(cmds[i++], "variable  width equal 2*${a}");
  sprintf(cmds[i++], "variable  zlow  equal 10*${a}-0.00001");
  sprintf(cmds[i++], "variable  zhigh equal 20*${a}");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++],
          "region whole block 0 ${width} 0 ${width} 0 ${zhigh}  "
          "units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 0  orient  z "
          "0 0 1");
  sprintf(cmds[i++],
          "region rdn  block INF INF INF INF  0.0000  ${zlow}  units box");
  sprintf(cmds[i++],
          "region rup  block INF INF INF INF ${zlow} ${zhigh} units box");
  sprintf(cmds[i++], "create_atoms 1 region rdn");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(cmds[i++], "thermo 100");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-10  1e-10  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  eS100 = lammps_get_thermo(lmp, ETOTAL);

  /*******************************************************************
   * calculate the bulk energy 100 pure bulk
   * *****************************************************************/
  i = 0;
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);
  sprintf(cmds[i++], "variable  width equal 2*${a}");
  sprintf(cmds[i++], "variable  zhigh equal 20*${a}");

  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++],
          "region whole block 0 ${width} 0 ${width} 0 ${zhigh}  "
          "units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 0  orient  z "
          "0 0 1");
  sprintf(cmds[i++], "create_atoms 1 region whole");

  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  sprintf(cmds[i++], "thermo 100");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-10  1e-10  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  eB100 = lammps_get_thermo(lmp, ETOTAL);
  double lx = lammps_get_thermo(lmp, LX);
  double ly = lammps_get_thermo(lmp, LY);
  areaS100 = lx * ly;
  exprs["suf100"] = 0.5 * (eS100 - 0.5 * eB100) / areaS100 * EV_A2toJ_M2;

  /*******************************************************************
   * calculate the energy 110 with surface
   * *****************************************************************/
  i = 0;
  //  --------------------- INITIALIZAITION ---------------------
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);
  sprintf(cmds[i++], "variable  widthx equal 4*${a}");
  sprintf(cmds[i++], "variable  widthy equal 4*${a}*sqrt(2)");
  sprintf(cmds[i++], "variable  zlow  equal  5*${a}*sqrt(2)-0.00001");
  sprintf(cmds[i++], "variable  zhigh equal 10*${a}");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++],
          "region whole block 0 ${widthx} 0 ${widthy} 0 ${zhigh}  "
          "units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 -1  orient  z "
          "0 1 1");
  sprintf(cmds[i++],
          "region rdn  block INF INF INF INF  0.0000  ${zlow} units box");
  sprintf(cmds[i++],
          "region rup  block INF INF INF INF ${zlow} ${zhigh} units box");
  sprintf(cmds[i++], "create_atoms 1 region rdn");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(cmds[i++], "thermo 100");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-10  1e-10  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  eS110 = lammps_get_thermo(lmp, ETOTAL);

  /*******************************************************************
   * calculate the energy 110 bulk
   * *****************************************************************/
  i = 0;
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);
  sprintf(cmds[i++], "variable  widthx equal 4*${a}");
  sprintf(cmds[i++], "variable  widthy equal 4*${a}*sqrt(2)");
  sprintf(cmds[i++], "variable  zhigh equal  5*${a}*sqrt(2)");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++],
          "region whole block 0 ${widthx} 0 ${widthy} 0 ${zhigh}  "
          "units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 0 0 orient y 0 1 -1  orient  z "
          "0 1 1");
  sprintf(cmds[i++], "create_atoms 1 region whole");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(cmds[i++], "thermo 100");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-10  1e-10  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  eB110 = lammps_get_thermo(lmp, ETOTAL);
  lx = lammps_get_thermo(lmp, LX);
  ly = lammps_get_thermo(lmp, LY);
  areaS110 = lx * ly;
  exprs["suf110"] = 0.5 * (eS110 - eB110) / areaS110 * EV_A2toJ_M2;

  /*******************************************************************
   * calculate the energy 111 with surface
   * *****************************************************************/
  i = 0;
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);
  sprintf(cmds[i++], "variable  widthx equal 1*${a}*sqrt(2)");
  sprintf(cmds[i++], "variable  widthy equal 2*${a}*sqrt(6)/2");
  sprintf(cmds[i++], "variable  zlow  equal  5*${a}*sqrt(3)-0.00001");
  sprintf(cmds[i++], "variable  zhigh equal  10*${a}*sqrt(3)");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++],
          "region whole block 0 ${widthx} 0 ${widthy} 0 ${zhigh}  "
          "units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 -1 0 orient y  1 1 -2 orient  "
          "z 1 1 1");
  sprintf(cmds[i++],
          "region rdn  block INF INF INF INF  0.0000  ${zlow} units box");
  sprintf(cmds[i++],
          "region rup  block INF INF INF INF ${zlow} ${zhigh} units box");
  sprintf(cmds[i++], "create_atoms 1 region rdn");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(cmds[i++], "thermo 100");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-10  1e-10  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  eS111 = lammps_get_thermo(lmp, ETOTAL);

  /*******************************************************************
   * calculate the energy 111 bulk
   * *****************************************************************/
  i = 0;
  //  --------------------- INITIALIZAITION ---------------------
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "units  metal");
  sprintf(cmds[i++], "dimension  3");
  sprintf(cmds[i++], "boundary p p p");
  sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);
  sprintf(cmds[i++], "variable  widthx equal 1*${a}*sqrt(2)");
  sprintf(cmds[i++], "variable  widthy equal 2*${a}*sqrt(6)/2");
  sprintf(cmds[i++], "variable  zhigh equal  5*${a}*sqrt(3)");

  // --------------------- ATOM DEFINITION ---------------------
  sprintf(cmds[i++], "lattice   bcc  ${a}");
  sprintf(cmds[i++],
          "region whole block 0 ${widthx} 0 ${widthy} 0 ${zhigh}  "
          "units box");
  sprintf(cmds[i++], "create_box  1  whole");
  sprintf(cmds[i++],
          "lattice   bcc  ${a} orient x 1 -1 0 orient y  1 1 -2  orient  "
          "z 1 1 1");
  sprintf(cmds[i++], "create_atoms 1 region whole");

  // --------------------- FORCE FIELDS ---------------------
  sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
  sprintf(cmds[i++], "neighbor 1.0 bin");
  sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

  // ---------------------- THERMO  --------------------------
  sprintf(cmds[i++], "thermo 100");
  sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal temp");

  // ---------------------- RELAX   --------------------------
  sprintf(cmds[i++], "min_style  cg");
  sprintf(cmds[i++], "minimize  1e-10  1e-10  100000  100000");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  eB111 = lammps_get_thermo(lmp, ETOTAL);
  lx = lammps_get_thermo(lmp, LX);
  ly = lammps_get_thermo(lmp, LY);
  areaS111 = lx * ly;
  exprs["suf111"] = 0.5 * (eS111 - eB111) / areaS111 * EV_A2toJ_M2;

  lammps_command(lmp, CLEAR);
}
