/*
 * @Author: chaomy
 * @Date:   2017-10-23 20:10:54
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-16 20:03:21
 */

#include "pfHome.h"
#include "pfLmpDrv.h"

// #define EPS 0.1
#define EPS 10
#define TEMPVAR 0.85
#define STEPVAR 2.0
// #define KMAX  500   // 1000
#define KMAX 100
#define INVSQRT2PI 0.39894228040143267794
#define GAUSS(a) (INVSQRT2PI * (exp(-((a) * (a)) / 2.0)))
#define NEPS 4
#define NSTEP 3  // 20
#define NTEMP (3 * 40)
#define RESFREQ 10

using std::cout;
using std::endl;
using std::string;
using std::vector;

void pfHome::randomize(vector<double>& vv, const int n,
                       const vector<double>& v) {
  const double width = fabs(randNormal());
  const double height = randNormal() * v[n];  // lower limit for step

  if (n > endps.back()) {  // change gradient
    // vv[n] += GAUSS(double(n) / width) * height;
  } else {
    int ff = 0;
    double w2 = 1.0 + 3 * width;                // original verison 4.0 * width
    while (ff < nfuncs) {                       //
      if (n >= startps[ff] && n < endps[ff]) {  // update function values
        chid = ff;
        for (int i = 0; i <= w2; i++) {
          int j = n + i;
          if (j < endps[ff]) vv[j] += GAUSS(double(i) / width) * height;
          j = n - i;
          if (j >= startps[ff]) vv[j] += GAUSS(double(i) / width) * height;
        }
        break;
      }
      ff++;
    }  // while
  }    //  else
}

void pfHome::simAnnealEAM() {
  int loopcnt = 0;
  int loopagain = 1;
  double T = dparams["temp"];
  double err = forceEAM(ini);
  double tmper = err, opter = err;

  vector<double> v(nvars, dparams["istep"]);
  vector<int> naccs(nvars, 0);
  vector<double> tmpvv = ini;
  vector<double> optvv = ini;
  vector<double> errold(NEPS, err);

  while (loopcnt < iparams["kmax"] && loopagain) {
    for (int m = 0; m < NTEMP; m++) {
      for (int j = 0; j < NSTEP; j++) {
        for (int h = 0; h < nvars; h++) {
          tmpvv = ini;
          randomize(tmpvv, h, v);
          tmper = forceEAM(tmpvv);

          if (tmper <= err) {
            ini = tmpvv;
            err = tmper;
            naccs[h]++;
            if (tmper < opter) {
              optvv = tmpvv;
              opter = tmper;
              writePot(optvv);
            }
          } else if (randUniform() < (exp((err - tmper) / T))) {
            ini = tmpvv;
            err = tmper;
            naccs[h]++;
          }
        }  // h
      }    // steps

      for (int n = 0; n < nvars; n++) { /* step adjustment */
        if (naccs[n] > 0.6 * NSTEP)
          v[n] *= (1 + STEPVAR * ((double)naccs[n] / NSTEP - 0.6) / 0.4);
        else if (naccs[n] < 0.4 * NSTEP)
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccs[n] / NSTEP) / 0.4);
        naccs[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\t%f\n", loopcnt, T, m + 1, error["frc"],
             error["punish"], opter);
      fflush(stdout);

      if (iparams["resfreq"] > 0 && (m + 1) % iparams["resfreq"] == 0) {
        updaterho(ini);
        writeLMPS(optvv);
        // if (rescaleRHO(ini) == 1) {
        //   printf("before rescale = %f\n", err);
        //   err = forceEAM(ini);
        //   printf("min = %f ; max = %f ; ave = %f \n", ominrho, omaxrho,
        //          oaverho);
        //   printf("after rescale = %f\n", err);
        // }
        //   if (rescaleEMF(ini) == 1) {
        //     printf("before rescale = %f\n", err);
        //     err = forceMEAM(ini);
        //     printf("min = %f ; max = %f ; ave = %f \n", ominrho, omaxrho,
        //            oaverho);
        //     printf("after rescale = %f\n", err);
        //   }
      } else if ((m + 1) % 20 == 0)
        updaterho(ini);
    }  // temp

    T *= TEMPVAR;
    loopcnt++;

    for (int i = 0; i < NEPS - 1; i++) errold[i] = errold[i + 1];

    errold[NEPS - 1] = err;
    loopagain = 0;

    for (int n = 0; n < NEPS - 1; n++) {
      if (fabs(err - errold[n]) > (EPS * err * 0.01)) {
        loopagain = 1;
        break;
      }
    }

    if (!loopagain && ((err - opter) > (EPS * err * 0.01))) {
      ini = optvv;
      err = opter;
      loopagain = 1;
    }
  }  // while

  ini = optvv;
  // recorderr.push_back(opter);
  // recordStage(scnt++);

  v.clear();
  naccs.clear();
  tmpvv.clear();
  optvv.clear();
  errold.clear();
}

void pfHome::simAnneal() {
  int loopcnt = 0;
  int loopagain = 1;
  double T = dparams["temp"];
  double err = forceMEAM(ini);
  double tmper = err, opter = err;

  vector<double> v(nvars, dparams["istep"]);
  vector<int> naccs(nvars, 0);
  vector<double> tmpvv = ini;
  vector<double> optvv = ini;
  vector<double> errold(NEPS, err);

  while (loopcnt < iparams["kmax"] && loopagain) {
    for (int m = 0; m < NTEMP; m++) {
      for (int j = 0; j < NSTEP; j++) {
        for (int h = 0; h < nvars; h++) {
          tmpvv = ini;
          randomize(tmpvv, h, v);
          tmper = forceMEAM(tmpvv);

          if (tmper <= err) {
            ini = tmpvv;
            err = tmper;
            naccs[h]++;
            if (tmper < opter) {
              optvv = tmpvv;
              opter = tmper;
              writePot(optvv);
            }
          } else if (randUniform() < (exp((err - tmper) / T))) {
            ini = tmpvv;
            err = tmper;
            naccs[h]++;
          }
        }  // h
      }    // steps

      for (int n = 0; n < nvars; n++) { /* step adjustment */
        if (naccs[n] > 0.6 * NSTEP)
          v[n] *= (1 + STEPVAR * ((double)naccs[n] / NSTEP - 0.6) / 0.4);
        else if (naccs[n] < 0.4 * NSTEP)
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccs[n] / NSTEP) / 0.4);
        naccs[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\t%f\t%f\n", loopcnt, T, m + 1, error["frc"],
             error["lat"], error["ela"], opter);
      // error["pv"], opter);
      printf("%f %f %f %f\n", exprs["lat"], exprs["c11"], exprs["c12"],
             exprs["c44"]);
      // -mpcf["pv"][5].strs[0],
      // -mpcf["pv"][10].strs[0]);  // 60 GPa
      fflush(stdout);

      if (iparams["resfreq"] > 0 && (m + 1) % iparams["resfreq"] == 0) {
        updaterhoMEAM(ini);
        if (rescaleEMF(ini) == 1) {
          printf("before rescale = %f\n", err);
          err = forceMEAM(ini);
          printf("min = %f ; max = %f ; ave = %f \n", ominrho, omaxrho,
                 oaverho);
          printf("after rescale = %f\n", err);
        }
      } else if ((m + 1) % 20 == 0)
        updaterhoMEAM(ini);

      // if (iparams["lmpfreq"] > 0 && (m + 1) % iparams["lmpfreq"] == 0){
      //   writeLMPS(optvv);
      //   lmpdrv->calPhy();
      // }
    }  // temp

    T *= TEMPVAR;
    loopcnt++;

    for (int i = 0; i < NEPS - 1; i++) errold[i] = errold[i + 1];

    errold[NEPS - 1] = err;
    loopagain = 0;

    for (int n = 0; n < NEPS - 1; n++) {
      if (fabs(err - errold[n]) > (EPS * err * 0.01)) {
        loopagain = 1;
        break;
      }
    }

    if (!loopagain && ((err - opter) > (EPS * err * 0.01))) {
      ini = optvv;
      err = opter;
      loopagain = 1;
    }
  }  // while

  ini = optvv;
  // recorderr.push_back(opter);
  // recordStage(scnt++);

  v.clear();
  naccs.clear();
  tmpvv.clear();
  optvv.clear();
  errold.clear();
}