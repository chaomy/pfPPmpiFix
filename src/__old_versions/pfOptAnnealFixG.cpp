/*
 * @Author: chaomy
 * @Date:   2017-10-23 20:10:54
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-10-30 21:58:52
 */

#include "pfHome.h"

#define EPS 0.1
#define TEMPVAR 0.85
#define STEPVAR 2.0
#define KMAX 1000
#define INVSQRT2PI 0.39894228040143267794
#define GAUSS(a) (INVSQRT2PI * (exp(-((a) * (a)) / 2.0)))
#define NEPS 4
#define NSTEP 3 // 20
#define NTEMP (3 * 40)
#define RESFREQ 10 

using std::cout;
using std::vector;
using std::string;

void pfHome::randomize(vector<Func>& ff, const int n, const vector<double>& v) {
  const double width = fabs(randNormal());
  const double height = randNormal() * v[n];

  for (int i = 0; i <= 4.0 * width; i++) {
    int j = n + i;
    double delta = GAUSS(double(i) / width) * height;
    for (int k = 0; k < nfuncs; k++) {
      if (j >= 0 && j < ff[k].npts) ff[k].yy[j] += delta;
      else j -= ff[k].npts;
    }
  }
}

void pfHome::simAnneal() {
  int loopcnt = 0;
  int loopagain = 1;
  double T = 0.01; 

  double err = forceEAM(funcs, 0);
  double tmp = err, opt = err;

  vector<double> v(nvars, 0.1);
  vector<int> naccs(nvars, 0);
  vector<Func> tmpff = funcs;
  vector<Func> optff = funcs;

  double errold[NEPS];

  for (int n = 0; n < NEPS; n++) errold[n] = err;

  while (loopcnt < KMAX && loopagain) {
    for (int m = 0; m < NTEMP; m++) {
      for (int j = 0; j < NSTEP; j++) {
        for (int h = 0; h < nvars; h++) {
          tmpff = funcs; 
          randomize(tmpff, h, v);
          tmp = forceEAM(tmpff, 0);
          if (tmp <= err) {
            funcs = tmpff;
            err = tmp;
            naccs[h]++;

            if (tmp < opt) {
              optff = tmpff;
              opt = tmp; 
              writePot();
            }

          } else if (randUniform() < (exp((err - tmp) / T))) {
            funcs = tmpff;
            err = tmp;
            naccs[h]++;
          }
        }  // h
      }    // steps

      /* step adjustment */
      for (int n = 0; n < nvars; n++) {
        if (naccs[n] > 0.6 * NSTEP)
          v[n] *= (1 + STEPVAR * ((double)naccs[n] / NSTEP - 0.6) / 0.4);
        else if (naccs[n] < 0.4 * NSTEP)
          v[n] /= (1 + STEPVAR * (0.4 - (double)naccs[n] / NSTEP) / 0.4);
        naccs[n] = 0;
      }

      printf("%3d\t%f\t%3d\t%f\t%f\n", loopcnt, T, m + 1, err, opt);
      fflush(stdout);

      // rescale
      if ((m + 1) % RESFREQ == 0){
        if (rescale(optff) != 0.0){
          printf("before rescale = %f\n", err);
          err = forceEAM(tmpff, 0);
          printf("after rescale = %f\n", err); 
        }  
      } // rescale 
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

    if (!loopagain && ((err - opt) > (EPS * err * 0.01))){
      funcs = optff; 
      err = opt;
      loopagain = 1;
    } 
  } // while 

  /* wake up others */
  forceEAM(funcs, 1);
}