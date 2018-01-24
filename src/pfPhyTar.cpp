/*
 * @Author: chaomy
 * @Date:   2017-11-14 14:24:20
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-26 23:19:05
 */

#include "pfHome.h"

void pfHome::initTargs() {
  error["tol"] = 0.0;

  // lattice
  targs["lat"] = 3.308;
  targs["abcc"] = 3.308;
  targs["afcc"] = 4.220;
  targs["ahcp"] = 2.90;
  targs["chcp"] = 5.27;

  // engy
  targs["ebcc"] = -10.0898;
  targs["efcc"] = -9.7707;
  targs["ehcp"] = -9.7952;

  targs["bcc2hcp"] = 0.295;
  targs["bcc2fcc"] = 0.320;

  // elastic
  targs["c11"] = 250.;
  targs["c12"] = 133.;
  targs["c44"] = 30.;

  // surf
  targs["suf100"] = 2.35;
  targs["suf110"] = 2.10;
  targs["suf111"] = 2.40;

  weigh["lat"] = 1e5;
  weigh["c11"] = 10;
  weigh["c12"] = 10;
  weigh["c44"] = 10;
  weigh["suf"] = 1.;
}