/*
 * @Author: chaomy
 * @Date:   2017-11-14 14:19:02
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-25 17:19:01
 */

#include "pfLmpDrv.h"

void pfLMPdrv::calElastic() {
  char cmds[1000][MAXLEN];

  double *ptc11 = (double *)NULL, *ptc12 = (double *)NULL,
         *ptc44 = (double *)NULL;

  int i = 0;
  //  --------------------- INITIALIZAITION ---------------------
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "variable up equal 1.0e-6");
  sprintf(cmds[i++], "variable atomjiggle equal 1.0e-5");

  sprintf(cmds[i++], "units   metal");
  sprintf(cmds[i++], "variable  cfac equal 1.0e-4");
  sprintf(cmds[i++], "variable  cunits string GPa");

  sprintf(cmds[i++], "variable  etol equal 0.0");
  sprintf(cmds[i++], "variable  ftol    equal 1.0e-10");
  sprintf(cmds[i++], "variable  maxiter equal 100");
  sprintf(cmds[i++], "variable  maxeval equal 1000");
  sprintf(cmds[i++], "variable  dmax    equal 1.0e-2");

  sprintf(cmds[i++], "variable a  equal  %.5f", exprs["lat"]);
  sprintf(cmds[i++], "boundary p p p");

  sprintf(cmds[i++], "lattice  bcc $a");
  sprintf(cmds[i++],
          "region   box prism  0  2.0  0  3.0  0  4.0  0.0  0.0  0.0");
  sprintf(cmds[i++], "create_box   1  box");
  sprintf(cmds[i++], "create_atoms 1  box");

  sprintf(cmds[i++], "mass  *  %f", pfhm->gdparams()["mass"]);

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff  * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*---------------------------------------------------------------
   * compute initial state
   * -------------------------------------------------------------*/
  sprintf(cmds[i++], "fix 3 all box/relax  aniso 0.0");
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal lx");
  sprintf(cmds[i++], "variable lx0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal ly");
  sprintf(cmds[i++], "variable ly0 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal lz");
  sprintf(cmds[i++], "variable lz0 equal ${tmp}");

  // These formulas define the derivatives w.r.t. strain components
  // Constants uses $, variables use v_
  sprintf(cmds[i++],
          "variable d1 equal -(v_pxx1-${pxx0})/(v_delta/v_len0)*${cfac}");
  sprintf(cmds[i++],
          "variable d2 equal -(v_pyy1-${pyy0})/(v_delta/v_len0)*${cfac}");
  sprintf(cmds[i++],
          "variable d3 equal -(v_pzz1-${pzz0})/(v_delta/v_len0)*${cfac}");
  sprintf(cmds[i++],
          "variable d4 equal -(v_pyz1-${pyz0})/(v_delta/v_len0)*${cfac}");
  sprintf(cmds[i++],
          "variable d5 equal -(v_pxz1-${pxz0})/(v_delta/v_len0)*${cfac}");
  sprintf(cmds[i++],
          "variable d6 equal -(v_pxy1-${pxy0})/(v_delta/v_len0)*${cfac}");

  sprintf(cmds[i++],
          "displace_atoms all random ${atomjiggle} ${atomjiggle} "
          "${atomjiggle} 87287 units box");

  sprintf(cmds[i++], "unfix 3");
  sprintf(cmds[i++], "write_restart restart.equil");

  /* ###########################################################################
   * uxx Perturbation
   * ###########################################################################*/
  sprintf(cmds[i++], "variable dir equal 1");
  sprintf(cmds[i++], "variable len0 equal ${lx0}");

  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor   1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify norm no");

  /*--- Negative deformation ---*/
  sprintf(cmds[i++], "variable delta   equal -${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal -${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal -${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal -${up}*yz");

  /*--- change box ---*/
  sprintf(cmds[i++],
          "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz "
          "delta ${deltaxz} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1neg equal ${d1}");
  sprintf(cmds[i++], "variable C2neg equal ${d2}");
  sprintf(cmds[i++], "variable C3neg equal ${d3}");
  sprintf(cmds[i++], "variable C4neg equal ${d4}");
  sprintf(cmds[i++], "variable C5neg equal ${d5}");
  sprintf(cmds[i++], "variable C6neg equal ${d6}");

  /*--- Reset box and simulation parameters ---*/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*--- Positive deformation ---*/
  sprintf(cmds[i++], "variable delta   equal ${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal ${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal ${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal ${up}*yz");

  sprintf(cmds[i++],
          "change_box all x delta 0 ${delta} xy delta ${deltaxy} xz "
          "delta ${deltaxz} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pe");
  sprintf(cmds[i++], "variable e1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal press");
  sprintf(cmds[i++], "variable p1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1pos equal ${d1}");
  sprintf(cmds[i++], "variable C2pos equal ${d2}");
  sprintf(cmds[i++], "variable C3pos equal ${d3}");
  sprintf(cmds[i++], "variable C4pos equal ${d4}");
  sprintf(cmds[i++], "variable C5pos equal ${d5}");
  sprintf(cmds[i++], "variable C6pos equal ${d6}");

  /*--- Combine positive and negative ---*/
  sprintf(cmds[i++], "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})");
  sprintf(cmds[i++], "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})");
  sprintf(cmds[i++], "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})");
  sprintf(cmds[i++], "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})");
  sprintf(cmds[i++], "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})");
  sprintf(cmds[i++], "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})");

  /*--- Delete dir to make sure it is not reused ---*/
  sprintf(cmds[i++], "variable dir delete");

  /* ###########################################################################
   * uyy Perturbation
   * ###########################################################################*/
  sprintf(cmds[i++], "variable dir equal 2");
  sprintf(cmds[i++], "variable len0 equal ${ly0}");
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor   1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify norm no");

  /*--- Negative deformation ---*/
  sprintf(cmds[i++], "variable delta   equal -${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal -${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal -${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal -${up}*yz");

  /*--- change box ---*/
  sprintf(cmds[i++],
          "change_box all y delta 0 ${delta} yz delta ${deltayz} remap "
          "units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1neg equal ${d1}");
  sprintf(cmds[i++], "variable C2neg equal ${d2}");
  sprintf(cmds[i++], "variable C3neg equal ${d3}");
  sprintf(cmds[i++], "variable C4neg equal ${d4}");
  sprintf(cmds[i++], "variable C5neg equal ${d5}");
  sprintf(cmds[i++], "variable C6neg equal ${d6}");

  /*--- Reset box and simulation parameters ---*/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*--- Positive deformation ---*/
  sprintf(cmds[i++], "variable delta   equal ${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal ${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal ${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal ${up}*yz");

  sprintf(cmds[i++],
          "change_box all y delta 0 ${delta} yz delta ${deltayz} remap "
          "units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pe");
  sprintf(cmds[i++], "variable e1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal press");
  sprintf(cmds[i++], "variable p1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1pos equal ${d1}");
  sprintf(cmds[i++], "variable C2pos equal ${d2}");
  sprintf(cmds[i++], "variable C3pos equal ${d3}");
  sprintf(cmds[i++], "variable C4pos equal ${d4}");
  sprintf(cmds[i++], "variable C5pos equal ${d5}");
  sprintf(cmds[i++], "variable C6pos equal ${d6}");

  /*--- Combine positive and negative ---*/
  sprintf(cmds[i++], "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})");
  sprintf(cmds[i++], "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})");
  sprintf(cmds[i++], "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})");
  sprintf(cmds[i++], "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})");
  sprintf(cmds[i++], "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})");
  sprintf(cmds[i++], "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})");

  /*--- Delete dir to make sure it is not reused ---*/
  sprintf(cmds[i++], "variable dir delete");

  /* ###########################################################################
   * uzz Perturbation
   * ###########################################################################*/
  sprintf(cmds[i++], "variable dir equal 3");
  sprintf(cmds[i++], "variable len0 equal ${lz0}");
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor   1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify norm no");

  /*--- Negative deformation ---*/
  sprintf(cmds[i++], "variable delta   equal -${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal -${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal -${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal -${up}*yz");

  /*--- change box ---*/
  sprintf(cmds[i++], "change_box all z delta 0 ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1neg equal ${d1}");
  sprintf(cmds[i++], "variable C2neg equal ${d2}");
  sprintf(cmds[i++], "variable C3neg equal ${d3}");
  sprintf(cmds[i++], "variable C4neg equal ${d4}");
  sprintf(cmds[i++], "variable C5neg equal ${d5}");
  sprintf(cmds[i++], "variable C6neg equal ${d6}");

  /*--- Reset box and simulation parameters ---*/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*--- Positive deformation ---*/
  sprintf(cmds[i++], "variable delta   equal ${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal ${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal ${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal ${up}*yz");

  sprintf(cmds[i++], "change_box all z delta 0 ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pe");
  sprintf(cmds[i++], "variable e1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal press");
  sprintf(cmds[i++], "variable p1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1pos equal ${d1}");
  sprintf(cmds[i++], "variable C2pos equal ${d2}");
  sprintf(cmds[i++], "variable C3pos equal ${d3}");
  sprintf(cmds[i++], "variable C4pos equal ${d4}");
  sprintf(cmds[i++], "variable C5pos equal ${d5}");
  sprintf(cmds[i++], "variable C6pos equal ${d6}");

  /*--- Combine positive and negative ---*/
  sprintf(cmds[i++], "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})");
  sprintf(cmds[i++], "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})");
  sprintf(cmds[i++], "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})");
  sprintf(cmds[i++], "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})");
  sprintf(cmds[i++], "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})");
  sprintf(cmds[i++], "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})");

  /*--- Delete dir to make sure it is not reused ---*/
  sprintf(cmds[i++], "variable dir delete");

  /* ###########################################################################
   * uyz Perturbation
   * ###########################################################################*/
  sprintf(cmds[i++], "variable dir equal 4");
  sprintf(cmds[i++], "variable len0 equal ${lz0}");
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor   1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify norm no");

  /*--- Negative deformation ---*/
  sprintf(cmds[i++], "variable delta   equal -${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal -${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal -${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal -${up}*yz");

  /*--- change box ---*/
  sprintf(cmds[i++], "change_box all yz delta ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1neg equal ${d1}");
  sprintf(cmds[i++], "variable C2neg equal ${d2}");
  sprintf(cmds[i++], "variable C3neg equal ${d3}");
  sprintf(cmds[i++], "variable C4neg equal ${d4}");
  sprintf(cmds[i++], "variable C5neg equal ${d5}");
  sprintf(cmds[i++], "variable C6neg equal ${d6}");

  /*--- Reset box and simulation parameters ---*/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*--- Positive deformation ---*/
  sprintf(cmds[i++], "variable delta   equal ${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal ${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal ${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal ${up}*yz");

  sprintf(cmds[i++], "change_box all yz delta ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pe");
  sprintf(cmds[i++], "variable e1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal press");
  sprintf(cmds[i++], "variable p1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1pos equal ${d1}");
  sprintf(cmds[i++], "variable C2pos equal ${d2}");
  sprintf(cmds[i++], "variable C3pos equal ${d3}");
  sprintf(cmds[i++], "variable C4pos equal ${d4}");
  sprintf(cmds[i++], "variable C5pos equal ${d5}");
  sprintf(cmds[i++], "variable C6pos equal ${d6}");

  /*--- Combine positive and negative ---*/
  sprintf(cmds[i++], "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})");
  sprintf(cmds[i++], "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})");
  sprintf(cmds[i++], "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})");
  sprintf(cmds[i++], "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})");
  sprintf(cmds[i++], "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})");
  sprintf(cmds[i++], "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})");

  /*--- Delete dir to make sure it is not reused ---*/
  sprintf(cmds[i++], "variable dir delete");

  /* ###########################################################################
   * uxz Perturbation
   * ###########################################################################*/
  sprintf(cmds[i++], "variable dir equal 5");
  sprintf(cmds[i++], "variable len0 equal ${lz0}");
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor   1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify norm no");

  /*--- Negative deformation ---*/
  sprintf(cmds[i++], "variable delta   equal -${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal -${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal -${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal -${up}*yz");

  /*--- change box ---*/
  sprintf(cmds[i++], "change_box all xz delta ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1neg equal ${d1}");
  sprintf(cmds[i++], "variable C2neg equal ${d2}");
  sprintf(cmds[i++], "variable C3neg equal ${d3}");
  sprintf(cmds[i++], "variable C4neg equal ${d4}");
  sprintf(cmds[i++], "variable C5neg equal ${d5}");
  sprintf(cmds[i++], "variable C6neg equal ${d6}");

  /*--- Reset box and simulation parameters ---*/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*--- Positive deformation ---*/
  sprintf(cmds[i++], "variable delta   equal ${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal ${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal ${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal ${up}*yz");

  sprintf(cmds[i++], "change_box all xz delta ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pe");
  sprintf(cmds[i++], "variable e1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal press");
  sprintf(cmds[i++], "variable p1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1pos equal ${d1}");
  sprintf(cmds[i++], "variable C2pos equal ${d2}");
  sprintf(cmds[i++], "variable C3pos equal ${d3}");
  sprintf(cmds[i++], "variable C4pos equal ${d4}");
  sprintf(cmds[i++], "variable C5pos equal ${d5}");
  sprintf(cmds[i++], "variable C6pos equal ${d6}");

  /*--- Combine positive and negative ---*/
  sprintf(cmds[i++], "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})");
  sprintf(cmds[i++], "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})");
  sprintf(cmds[i++], "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})");
  sprintf(cmds[i++], "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})");
  sprintf(cmds[i++], "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})");
  sprintf(cmds[i++], "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})");

  /*--- Delete dir to make sure it is not reused ---*/
  sprintf(cmds[i++], "variable dir delete");

  /* ###########################################################################
   * uxy Perturbation
   * ###########################################################################*/
  sprintf(cmds[i++], "variable dir equal 6");
  sprintf(cmds[i++], "variable len0 equal ${ly0}");
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor   1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify norm no");

  /*--- Negative deformation ---*/
  sprintf(cmds[i++], "variable delta   equal -${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal -${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal -${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal -${up}*yz");

  /*--- change box ---*/
  sprintf(cmds[i++], "change_box all xy delta ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1neg equal ${d1}");
  sprintf(cmds[i++], "variable C2neg equal ${d2}");
  sprintf(cmds[i++], "variable C3neg equal ${d3}");
  sprintf(cmds[i++], "variable C4neg equal ${d4}");
  sprintf(cmds[i++], "variable C5neg equal ${d5}");
  sprintf(cmds[i++], "variable C6neg equal ${d6}");

  /*--- Reset box and simulation parameters ---*/
  sprintf(cmds[i++], "clear");
  sprintf(cmds[i++], "box tilt large");
  sprintf(cmds[i++], "read_restart restart.equil");

  /****** potential ******/
  sprintf(cmds[i++], "pair_style   %s", sttag["pairstyle"].c_str());
  sprintf(cmds[i++], "pair_coeff   * *  ./%s  %s", sttag["lmpfile"].c_str(),
          sttag["elem"].c_str());
  /****** neighbor ******/
  sprintf(cmds[i++], "neighbor  1.0  nsq");
  sprintf(cmds[i++],
          "neigh_modify   once  no   every  1  delay  0  check  yes");
  /****** minimization ******/
  sprintf(cmds[i++], "min_style    cg");
  sprintf(cmds[i++], "min_modify  dmax ${dmax} line quadratic");
  /****** output ******/
  sprintf(cmds[i++], "thermo  1");
  sprintf(cmds[i++],
          "thermo_style  custom step temp pe press pxx pyy pzz pxy pxz "
          "pyz lx ly lz vol");
  sprintf(cmds[i++], "thermo_modify  norm no");

  /*--- Positive deformation ---*/
  sprintf(cmds[i++], "variable delta   equal ${up}*${len0}");
  sprintf(cmds[i++], "variable deltaxy equal ${up}*xy");
  sprintf(cmds[i++], "variable deltaxz equal ${up}*xz");
  sprintf(cmds[i++], "variable deltayz equal ${up}*yz");

  sprintf(cmds[i++], "change_box all xy delta ${delta} remap units box");

  /*--- relaxation ---*/
  sprintf(cmds[i++], "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");

  /*--- new stress tensor ---*/
  sprintf(cmds[i++], "variable tmp equal pe");
  sprintf(cmds[i++], "variable e1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal press");
  sprintf(cmds[i++], "variable p1 equal ${tmp}");

  sprintf(cmds[i++], "variable tmp equal pxx");
  sprintf(cmds[i++], "variable pxx1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyy");
  sprintf(cmds[i++], "variable pyy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pzz");
  sprintf(cmds[i++], "variable pzz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxy");
  sprintf(cmds[i++], "variable pxy1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pxz");
  sprintf(cmds[i++], "variable pxz1 equal ${tmp}");
  sprintf(cmds[i++], "variable tmp equal pyz");
  sprintf(cmds[i++], "variable pyz1 equal ${tmp}");

  /*--- Compute elastic constant from pressure tensor ---*/
  sprintf(cmds[i++], "variable C1pos equal ${d1}");
  sprintf(cmds[i++], "variable C2pos equal ${d2}");
  sprintf(cmds[i++], "variable C3pos equal ${d3}");
  sprintf(cmds[i++], "variable C4pos equal ${d4}");
  sprintf(cmds[i++], "variable C5pos equal ${d5}");
  sprintf(cmds[i++], "variable C6pos equal ${d6}");

  /*--- Combine positive and negative ---*/
  sprintf(cmds[i++], "variable C1${dir} equal 0.5*(${C1neg}+${C1pos})");
  sprintf(cmds[i++], "variable C2${dir} equal 0.5*(${C2neg}+${C2pos})");
  sprintf(cmds[i++], "variable C3${dir} equal 0.5*(${C3neg}+${C3pos})");
  sprintf(cmds[i++], "variable C4${dir} equal 0.5*(${C4neg}+${C4pos})");
  sprintf(cmds[i++], "variable C5${dir} equal 0.5*(${C5neg}+${C5pos})");
  sprintf(cmds[i++], "variable C6${dir} equal 0.5*(${C6neg}+${C6pos})");

  /*--- Delete dir to make sure it is not reused ---*/
  sprintf(cmds[i++], "variable dir delete");

  /*--- output final values ---*/

  sprintf(cmds[i++], "variable C11all equal ${C11}");
  sprintf(cmds[i++], "variable C22all equal ${C22}");
  sprintf(cmds[i++], "variable C33all equal ${C33}");

  sprintf(cmds[i++], "variable C12all equal 0.5*(${C12}+${C21})");
  sprintf(cmds[i++], "variable C13all equal 0.5*(${C13}+${C31})");
  sprintf(cmds[i++], "variable C23all equal 0.5*(${C23}+${C32})");

  sprintf(cmds[i++], "variable C44all equal ${C44}");
  sprintf(cmds[i++], "variable C55all equal ${C55}");
  sprintf(cmds[i++], "variable C66all equal ${C66}");

  sprintf(cmds[i++], "variable C14all equal 0.5*(${C14}+${C41})");
  sprintf(cmds[i++], "variable C15all equal 0.5*(${C15}+${C51})");
  sprintf(cmds[i++], "variable C16all equal 0.5*(${C16}+${C61})");

  sprintf(cmds[i++], "variable C24all equal 0.5*(${C24}+${C42})");
  sprintf(cmds[i++], "variable C25all equal 0.5*(${C25}+${C52})");
  sprintf(cmds[i++], "variable C26all equal 0.5*(${C26}+${C62})");

  sprintf(cmds[i++], "variable C34all equal 0.5*(${C34}+${C43})");
  sprintf(cmds[i++], "variable C35all equal 0.5*(${C35}+${C53})");
  sprintf(cmds[i++], "variable C36all equal 0.5*(${C36}+${C63})");

  sprintf(cmds[i++], "variable C45all equal 0.5*(${C45}+${C54})");
  sprintf(cmds[i++], "variable C46all equal 0.5*(${C46}+${C64})");
  sprintf(cmds[i++], "variable C56all equal 0.5*(${C56}+${C65})");

  sprintf(cmds[i++],
          "variable C11cubic equal (${C11all}+${C22all}+${C33all})/3.0");
  sprintf(cmds[i++],
          "variable C12cubic equal (${C12all}+${C13all}+${C23all})/3.0");
  sprintf(cmds[i++],
          "variable C44cubic equal (${C44all}+${C55all}+${C66all})/3.0");

  sprintf(cmds[i++],
          "variable bulkmodulus equal (${C11cubic}+2*${C12cubic})/3.0");
  sprintf(cmds[i++], "variable shearmodulus1 equal ${C44cubic}");
  sprintf(cmds[i++],
          "variable shearmodulus2 equal (${C11cubic}-${C12cubic})/2.0");
  sprintf(cmds[i++],
          "variable poissonratio equal 1.0/(1.0+${C11cubic}/${C12cubic})");

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, cmds[iter]);

  sprintf(cmds[0], "C11cubic");
  sprintf(cmds[1], "C12cubic");
  sprintf(cmds[2], "C44cubic");

  ptc11 = (double *)lammps_extract_variable(lmp, cmds[0], NULL);
  ptc12 = (double *)lammps_extract_variable(lmp, cmds[1], NULL);
  ptc44 = (double *)lammps_extract_variable(lmp, cmds[2], NULL);

  exprs["c11"] = *ptc11;
  exprs["c12"] = *ptc12;
  exprs["c44"] = *ptc44;

  double de1 = abs(exprs["c11"] - targs["c11"]);
  double de2 = abs(exprs["c12"] - targs["c12"]);
  double de3 = abs(exprs["c44"] - targs["c44"]);

  error["c11"] = de1 / targs["c11"];
  error["c12"] = de2 / targs["c12"];
  error["c44"] = de3 / targs["c44"];

  free(ptc11);
  free(ptc12);
  free(ptc44);

  sprintf(cmds[0], "clear");
  lammps_command(lmp, cmds[0]);
}
