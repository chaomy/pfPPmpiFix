/*
 * @Author: chaomy
 * @Date:   2017-12-13 09:53:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-23 16:32:14
 */

#include "pfHome.h"

using std::cout;
using std::endl;
using std::to_string;
using std::vector;

void pfHome::loopBwth() {
  for (int i = 1; i < 20; i++) {
    double tm = 0.2 + 0.2 * i;
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        for (int it : {0, 1, 2})
          atm.fweigh[it] = exp(-tm * atm.frc[it] * atm.frc[it]);

    vector<double> mtm(3, 0);
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        for (int ii : {0, 1, 2}) mtm[ii] += fabs(atm.frc[ii]) * atm.fweigh[ii];

    for (int ii : {0, 1, 2}) mtm[ii] /= ftn;

    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        for (int ii : {0, 1, 2}) atm.fweigh[ii] *= (mfrc[ii] / mtm[ii]);

    string fm("fc" + to_string(i) + ".txt");
    FILE* fid = fopen(fm.c_str(), "w");
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms)
        fprintf(fid, "%f %f %f %f %f %f\n", atm.frc[X], atm.frc[Y], atm.frc[Z],
                atm.frc[X] * atm.fweigh[X], atm.frc[Y] * atm.fweigh[Y],
                atm.frc[Z] * atm.fweigh[Z]);
    fclose(fid);
  }  // i = 1 to 9
}

void pfHome::forceDis() {
  vector<double> mtm(3, 0);
  for (Config& tmpc : configs) {
    for (pfAtom& atm : tmpc.atoms)
      for (int ii : {0, 1, 2}) mtm[ii] += fabs(atm.frc[ii]) * atm.fweigh[ii];
  }
  // for (int ii : {0, 1, 2}) mtm[ii] /= ftn;
  // for (Config& tmpc : configs) {
  //   for (pfAtom& atm : tmpc.atoms)
  //     for (int ii : {0, 1, 2}) atm.fweigh[ii] *= (mfrc[ii] / mtm[ii]);
  // }
  fsm = square11(mfrc[X]) + square11(mfrc[Y]) + square11(mfrc[Z]);
}

void pfHome::deleteAtoms() { /* delete some atoms */
  for (Config& tmpc : configs) {
    vector<pfAtom> ntm;
    for (pfAtom& atm : tmpc.atoms)
      if (atm.nneighs == N3 && atm.nneighsFull == N3) ntm.push_back(atm);
    tmpc.atoms = ntm;
    tmpc.natoms = tmpc.atoms.size();
  }
}

void pfHome::cutoffNeighs() { /* for chossing a cutoff*/
  FILE* fid = fopen("cutoff.txt", "w");
  rocut = 4.60;

  /* loop over cutoff */
  while (rocut < 6.40) {
    for (Config& tmpc : configs)
      for (pfAtom& atm : tmpc.atoms) atm.neighsFull.clear();

    initNeighsFull();
    fprintf(fid, "%0.3f ", rocut);
    for (Config& tmpc : configs) {
      int aven = 0;
      for (pfAtom& atm : tmpc.atoms) aven += atm.nneighsFull;
      fprintf(fid, "%03.1f ", float(aven) / tmpc.natoms);
    }

    fprintf(fid, "\n");
    rocut += 0.01;
  }

  fclose(fid);
}