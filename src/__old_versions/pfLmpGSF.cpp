/*
 * @Author: chaomy
 * @Date:   2017-11-10 14:40:55
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-26 10:42:52
 */

#include "pfLmpDrv.h"

/*******************************************************************
 * gsf curve  [111](110); [111](211); [111](123) serial version
 * *****************************************************************/
using std::vector;

void pfLMPdrv::calGSF() {
  char cmds[100][MAXLEN];

  const vector<double>& dft110 = mele.gsf[sttag["elem"] + "111z110"];
  const vector<double>& dft211 = mele.gsf[sttag["elem"] + "111z211"];
  const int npts = dft110.size();
  const double dlt = 1. / (npts - 1);

  vector<double>& gz110 = lgsf["111z110"] = vector<double>(npts, 0.);
  vector<double>& gz211 = lgsf["111z211"] = vector<double>(npts, 0.);

  vector<double>& ez110 = lgsf["111e110"] = vector<double>(npts, 0.);
  vector<double>& ez211 = lgsf["111e211"] = vector<double>(npts, 0.);

  double area = 1.0, equi = 0.0, disp = 0.0;
  int i;

  /*******************************************************************
   * gsf curve [111](110);
   * x [-1  1  1]
   * y [1, -1  2]
   * z [1,  1, 0]
   * *****************************************************************/
  for (int it = 0; it < npts; it++, disp += dlt) {
    i = 0;

    //  --------------------- INITIALIZAITION ---------------------
    sprintf(cmds[i++], "clear");
    sprintf(cmds[i++], "units  metal");
    sprintf(cmds[i++], "dimension  3");
    sprintf(cmds[i++], "boundary p p p");
    sprintf(cmds[i++], "atom_style atomic");
    sprintf(cmds[i++], "variable  a equal  %.7f", exprs["lat"]);
    sprintf(cmds[i++], "variable  wx  equal ${a}*sqrt(3)/2");
    sprintf(cmds[i++], "variable  wy  equal ${a}*sqrt(6)");
    sprintf(cmds[i++], "variable  zlo  equal  ${a}*sqrt(2)*10");
    sprintf(cmds[i++], "variable  zhi  equal  ${a}*sqrt(2)*20");
    sprintf(cmds[i++], "variable  ztl  equal ${a}*sqrt(2)*23");
    sprintf(cmds[i++], "variable  xdis equal  %f*${wx}", disp);

    // --------------------- ATOM DEFINITION ---------------------
    sprintf(cmds[i++], "lattice bcc  ${a}");
    sprintf(cmds[i++],
            "region whole  block  0  ${wx}  0  ${wy}  0  ${ztl} units box");
    sprintf(cmds[i++], "create_box  1   whole");
    sprintf(
        cmds[i++],
        "lattice   bcc  ${a} orient x -1 1 1 orient y 1 -1  2 orient  z 1 1 0");
    sprintf(cmds[i++],
            "region mybody  block  0  ${wx}  0  ${wy}  0  ${zhi}  units box");
    sprintf(cmds[i++], "create_atoms   1  region   mybody");
    sprintf(cmds[i++],
            "region      top   block   0  ${wx}  0  ${wy}  ${zlo}  ${zhi} "
            "units box");
    sprintf(cmds[i++], "group    gtop  region   top");

    // --------------------- FORCE FIELDS ---------------------
    sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
    sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
            sttag["elem"].c_str());
    sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
    sprintf(cmds[i++], "neighbor 1.0 bin");
    sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

    // ---------------------- THERMO  --------------------------
    sprintf(cmds[i++], "thermo 5000");
    sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal");
    sprintf(cmds[i++],
            "displace_atoms  gtop    move   ${xdis}  0  0 units box");

    // ---------------------- RELAX   --------------------------
    sprintf(cmds[i++], "fix  1  all  setforce  0  0  NULL");
    sprintf(cmds[i++], "min_style  cg");
    sprintf(cmds[i++], "minimize  1e-12  1e-12  100000  100000");
    sprintf(cmds[i++], "unfix  1");

    for (int it = 0; it < i; it++) lammps_command(lmp, cmds[it]);

    gz110[it] = lammps_get_thermo(lmp, ETOTAL);

    if (it == 0) {
      area = lammps_get_thermo(lmp, LX) * lammps_get_thermo(lmp, LY);
      equi = gz110[0];
    }

    gz110[it] = (gz110[it] - equi) / area;
    lammps_command(lmp, CLEAR);
  }

  /*******************************************************************
   * gsf curve [111](211) 20 points
   * x [-1  1  1]
   * y [ 0 -1  1
   * z [ 2  1  1]
   * *****************************************************************/
  disp = 0.0;
  for (int it = 0; it < npts; it++, disp += dlt) {
    i = 0;

    //  --------------------- INITIALIZAITION ---------------------
    sprintf(cmds[i++], "clear");
    sprintf(cmds[i++], "units  metal");
    sprintf(cmds[i++], "dimension  3");
    sprintf(cmds[i++], "boundary p p p");
    sprintf(cmds[i++], "atom_style atomic");
    sprintf(cmds[i++], "variable  a equal    %.7f", exprs["lat"]);
    sprintf(cmds[i++], "variable  wx  equal  ${a}*sqrt(3)/2");
    sprintf(cmds[i++], "variable  wy  equal  ${a}*sqrt(2)");
    sprintf(cmds[i++], "variable  zlo equal  ${a}*sqrt(6)*8");
    sprintf(cmds[i++], "variable  zhi equal  ${a}*sqrt(6)*16");
    sprintf(cmds[i++], "variable  ztl  equal  ${a}*sqrt(6)*22");
    sprintf(cmds[i++], "variable  xdis equal %f*${wx}", disp);

    // --------------------- ATOM DEFINITION ---------------------
    sprintf(cmds[i++], "lattice bcc  ${a}");
    sprintf(cmds[i++],
            "region  whole  block  0  ${wx}  0  ${wy}  0  ${ztl}  units box");
    sprintf(cmds[i++], "create_box  1  whole");
    sprintf(cmds[i++],
            "lattice  bcc  ${a} orient x -1 1 1  orient y 0 -1 1  orient  z "
            "2 1 1");
    sprintf(cmds[i++],
            "region   mybody  block  0  ${wx}  0  ${wy}  0  ${zhi}  units box");
    sprintf(cmds[i++], "create_atoms   1  region   mybody");
    sprintf(cmds[i++],
            "region  top   block   0  ${wx}  0  ${wy}  ${zlo}  ${zhi}  units "
            "box");
    sprintf(cmds[i++], "group  gtop  region   top");

    // --------------------- FORCE FIELDS ---------------------
    sprintf(cmds[i++], "pair_style  %s", sttag["pairstyle"].c_str());
    sprintf(cmds[i++], "pair_coeff  *  *  %s %s", sttag["lmpfile"].c_str(),
            sttag["elem"].c_str());
    sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);
    sprintf(cmds[i++], "neighbor 1.0 bin");
    sprintf(cmds[i++], "neigh_modify  every 1  delay  0 check yes");

    // ---------------------- THERMO  --------------------------
    sprintf(cmds[i++], "thermo 5000");
    sprintf(cmds[i++], "thermo_style  custom  lx  ly  lz step  etotal");
    sprintf(cmds[i++],
            "displace_atoms  gtop    move   ${xdis}  0  0  units box");

    // ---------------------- RELAX   --------------------------
    sprintf(cmds[i++], "fix  1  all  setforce  0  0  NULL");
    sprintf(cmds[i++], "min_style  cg");
    sprintf(cmds[i++], "minimize  1e-12  1e-12  100000  100000");
    sprintf(cmds[i++], "unfix  1");

    for (int it = 0; it < i; it++) lammps_command(lmp, cmds[it]);

    gz211[it] = lammps_get_thermo(lmp, ETOTAL);

    if (it == 0) {
      area = lammps_get_thermo(lmp, LX) * lammps_get_thermo(lmp, LY);
      equi = gz211[0];
    }

    gz211[it] = (gz211[it] - equi) / area;
    lammps_command(lmp, CLEAR);
  }

  for (int j = 0; j < npts; j++) {
    if (dft110[j] != 0.0) ez110[j] = fabs(gz110[j] - dft110[j]) / dft110[j];
    if (dft211[j] != 0.0) ez211[j] = fabs(gz211[j] - dft211[j]) / dft110[j];

    printf("%.4f %.4f %.8f\n", gz110[j], dft110[j], ez110[j]);
    printf("%.4f %.4f %.8f\n", gz211[j], dft211[j], ez211[j]);
  }
}
