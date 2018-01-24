/*
 * @Author: chaomy
 * @Date:   2017-11-10 14:44:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-16 21:21:31
 */

#include "pfLmpDrv.h"

#define SQRT2 1.4142135623730951
#define SMALLSTS 0.05  // Stress threshold 50 Mpa (0.05 Gpa) sqrt(2) * 0.05
#define MAXSTEP 500

using std::cout;
using std::endl;

inline void pmtx(const vector<vector<double>>& m) {
  for (unsigned int i = 0; i < m.size(); i++) {
    for (unsigned int j = 0; j < m.front().size(); j++) cout << m[i][j] << " ";
    cout << endl;
  }
}

inline void pvec(const vector<double>& v) {
  for (unsigned int i = 0; i < v.size(); i++) cout << v[i] << " ";
  cout << endl;
}

inline void matcopy(vector<vector<double>>& m, const double a[3][3]) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) m[i][j] = a[i][j];
}

inline vector<vector<double>> mat33mlt(const vector<vector<double>>& a,
                                       const vector<vector<double>>& b) {
  vector<vector<double>> c(a.size(), vector<double>(b.back().size(), 0));
  for (int i = 0; i < a.size(); i++)
    for (int j = 0; j < b.back().size(); j++) {
      c[i][j] = 0.0;
      for (int k = 0; k < a.front().size(); k++) c[i][j] += a[i][k] * b[k][j];
    }
  return c;
}

void pfLMPdrv::calIten() {
  vector<string> ms({"tpath", "opath"});
  int npts = 10;
  for (string tag : ms) {
    for (int i = 0; i < npts; i++) calItenOptLin(0.02 * i, tag);

    for (int i = npts - 1; i >= 0; i--) {
      iegyv[tag][i] -= iegyv[tag].front();
      iegyv[tag + "dft"][i] -= iegyv[(tag + "dft")].front();
    }

    for (int i = 0; i < npts; i++)
      cout << iegyv[tag][i] << " " << iegyv[tag + "dft"][i] << endl;
  }
}

void pfLMPdrv::calItenOptLin(const double& dlt, const string& tag) {
  vector<vector<double>> pbis(3, vector<double>(3, 0));
  vector<vector<double>> strm(3, vector<double>(3, 0));
  vector<double> coeff(3, 1.);
  vector<double> stsv(6);

  const ItenT& idft = pfhm->mele.itdftm[dlt][tag];
  double s12 = pfhm->mele.cijm["Nb"].sij[1], egy = 0.0;
  int cc = 0;

  /* initialize the perfect cell  */
  if (!tag.compare("tpath")) {
    for (int it = 0; it < 3; it++) pbis[it][it] = 1.0;
  } else if (!tag.compare("opath")) {
    pbis[0][0] = 1.0;
    pbis[1][1] = SQRT2;
    pbis[2][2] = SQRT2;
  }

  /** initialize strain  **/
  for (int it = 0; it < 3; it++) strm[it][it] = idft.strm[it][it];

  bsm = mat33mlt(strm, pbis); /** add strain on bsm **/
  calItenRun(tag, egy, stsv);

  while (true) {
    if (fmax(fabs(stsv[1]), fabs(stsv[2])) < SMALLSTS) {
      iegyv[tag].push_back(egy);
      isxxv[tag].push_back(stsv[0]);
      iegyv[(tag + "dft")].push_back(idft.egy);
      isxxv[(tag + "dft")].push_back(idft.stsv[0]);
      break;
    } else if (cc > MAXSTEP) {
      printf(" run too many times, tune your parameters \n");
      break;
    } else {
      for (int it = 1; it <= 2; it++) { /* ideal tensile only yy, zz maters */
        if (fabs(stsv[it]) > SMALLSTS)
          strm[it][it] += coeff[it] * s12 * stsv[it]; /* modify strain */
        else
          strm[it][it] += 0.00001 * s12 * stsv[it]; /* satisfies (small step) */
      }
      bsm = mat33mlt(strm, pbis); /** add strain on bsm **/
      calItenRun(tag, egy, stsv);
      cc++;
    } /*** err < EPSILON ***/
  }   /*** while true ***/
}