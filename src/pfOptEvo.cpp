/*
 * @Author: chaomy
 * @Date:   2017-12-30 14:13:52
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-13 12:21:11
 */

#include "pfHome.h"

#define EVOTH 1e-3
#define FL 0.1
#define FU 0.9
#define TF 0.1
#define TC 0.1

using std::cout;
using std::endl;
using std::vector;

void pfHome::diffEvo() {
  int a = 0, b = 0, c = 0, cn = 0;
  double mi = 10e10, mx = 0.0, sm = 0.0;
  double cr = 1e3;

  vector<double> hbd({3.0, 1.0, -0.0, 10.0});
  vector<double> lbd({-0.5, 0.0, -10.0, -10.0});

  int np = 15 * nvars;

  vector<vector<double>> p1(np, ini);
  vector<vector<double>> p2(np, ini);

  vector<double> fv(np, 0.);  // F
  vector<double> cv(np, 0.);  // CR

  vector<double> cs(np, 0.);     // cost vector
  vector<double> tt(nvars, 0.);  // trial
  vector<double> bs(ini.begin(), ini.end());

  for (int i = 1; i < np; i++) { /* init pop */
    int kk = 0;
    for (int j = 0; j < nfuncs; j++) {
      for (int l = 0; l < funcs[j].npts; l++)
        p1[i][kk++] = lbd[j] + randUniform() * (hbd[j] - lbd[j]);
    }
    while (kk < nvars)
      p1[i][kk++] = lbd.back() + randUniform() * (hbd.back() - lbd.back());

    for (double& e : p1[i]) cout << e << endl;
    fv[i] = FL + randUniform() * FU;
    cv[i] = randUniform();
  }

  for (int i = 0; i < np; i++) {
    // cs[i] = forceMEAM(p1[i]);
    cs[i] = forceEAM(p1[i]);
    bs = cs[i] < mi ? p1[i] : bs;
    mi = fmin(mi, cs[i]);
    mx = fmax(mx, cs[i]);
    sm += cs[i];
  }

  cr = mx - mi;
  printf("Loops\t\tOptimum\t\tAverage error sum\t\tMax-Min\n");
  printf("%5d\t\t%15f\t%20f\t\t%.2e\n", cn, mi, sm / (np), cr);

  while (cr >= EVOTH && mi >= EVOTH) {
    mx = 0.0;
    for (int i = 0; i < np; i++) { /* new populations */
      do
        a = (int)floor(randUniform() * np);
      while (a == i);

      do
        b = (int)floor(randUniform() * np);
      while (b == a || b == i);

      do
        c = (int)floor(randUniform() * np);
      while (c == a || c == b || c == i);

      /* adjust parameters */
      double tf;
      double tc;
      if (randUniform() < TF)
        tf = FL + randUniform() * FU;
      else
        tf = fv[i];

      if (randUniform() < TC)
        tc = randUniform();
      else
        tc = cv[i];

      int j = (int)floor(randUniform() * nvars);

      for (int k = 1; k <= nvars; k++) {
        if (randUniform() < tc || k == j) {
          /* DE/rand/1/exp */
          tt[j] = p1[c][j] + tf * (p1[a][j] - p1[b][j]);
          /* DE/best/1/exp */
          // tt[j] = bs[j] + tf * (p1[a][j] - p1[b][j]);
        } else
          tt[j] = p1[i][j];
        j = (j + 1) % nvars;
      }

      // double cc = forceMEAM(tt);
      double cc = forceEAM(tt);

      if (cc < mi) {
        bs = tt;
        writePot(bs);
        mi = cc;
      }

      if (cc <= cs[i]) {
        p2[i] = tt;
        cs[i] = cc;
        fv[i] = tf;
        cv[i] = tc;
      } else
        p2[i] = p1[i];

      if (cs[i] > mx) mx = cs[i];
    }

    sm = 0.0;
    for (double ee : cs) sm += ee;
    cr = mx - mi;
    printf("%5d\t\t%15f\t%20f\t\t%.2e\n", ++cn, mi, sm / np, cr);
    for (int i = 0; i < np; i++) p1[i] = p2[i];
  }  // while
}