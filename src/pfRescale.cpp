/*
 * @Author: chaomy
 * @Date:   2017-10-30 21:34:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2018-01-24 00:04:33
 */

#include "pfHome.h"
#define NNEIGH 26

// unfinished
void pfHome::shiftRHO(vector<double>& vv) {
  double shift = (0 - ominrho) / NNEIGH;
  for (unsigned int i = 0; i < funcs[RHO].yy.size(); i++)
    vv[startps[RHO] + i] += shift;
}

void pfHome::shiftEMF(double shift) {
  for (unsigned int i = 0; i < funcs[EMF].yy.size(); i++)
    funcs[EMF].xx[i] += shift;
}

int pfHome::rescaleEMF(arma::mat& in) {
  Func& ff = funcs[EMF];
  int npts = ff.npts;
  double delt = ff.step;
  if (fabs(ominrho - ff.xx.front()) > 3. * delt ||
      fabs(omaxrho - ff.xx.back()) > 3. * delt) {
    arma::mat iterate = decodev(in);  // -> [a, b]
    double ndelt = (omaxrho - ominrho + 2. * delt) / (npts - 1);
    for (int i = 0; i < npts; i++)
      iterate[startps[EMF] + i] = ff.s(ominrho - delt + i * ndelt);
    for (int i = 0; i < npts; i++) ff.xx[i] = ominrho - delt + i * ndelt;
    ff.step = ndelt;
    in = encodev(iterate);
    return 1;
  } else
    return 0;
}

int pfHome::rescaleEMF(vector<double>& vv) {  // check if rescale
  Func& ff = funcs[EMF];
  int npts = ff.npts;
  double delt = ff.step;

  if (fabs(ominrho - ff.xx.front()) > delt ||
      fabs(omaxrho - ff.xx.back()) > delt) {
    for (int i = 0; i < npts; i++) ff.yy[i] = vv[startps[EMF] + i];

    ff.g1.front() = vv[nvars - 4];
    ff.g1.back() = vv[nvars - 2];

    splineNe(ff, gradRight[EMF]);  // update splines
    double ndelt = (omaxrho - ominrho) / (npts - 1);

    /* update y values (to vv) */
    for (int i = 0; i < npts; i++)
      splint(ff, ominrho + i * ndelt, vv[startps[EMF] + i]);

    /* update xx values */
    for (int i = 0; i < npts; i++) ff.xx[i] = ominrho + i * ndelt;
    ff.step = ndelt;
    return (1);
  } else
    return (0);
}

int pfHome::rescaleRHO(vector<double>& vv) {
  printf("max %f hilim %f min %f lolim %f oaverho %f\n", omaxrho, hirho,
         ominrho, lorho, oaverho);
  double aa = (omaxrho - ominrho) / (hirho - lorho);
  if (aa > 1.2 || aa < 1.00) {  // rescale + shift
    double coeff = 1.01 / (aa);
    // rescale
    for (unsigned int i = 0; i < funcs[RHO].yy.size() + 2; i++)
      vv[startps[RHO] + i] *= coeff;
    // shift
    double shift = -coeff * (omaxrho + ominrho - hirho - lorho) / (2 * NNEIGH);
    for (unsigned int i = 0; i < funcs[RHO].yy.size(); i++)
      vv[startps[RHO] + i] += shift;
    printf("a = %f ; rescale %f ; shift %f \n", aa, coeff, shift);

  } else if ((omaxrho - hirho) * (ominrho - lorho) > 0) {
    double shift = -(omaxrho + ominrho - hirho - lorho) / (2 * NNEIGH);
    for (unsigned int i = 0; i < funcs[RHO].yy.size(); i++)
      vv[startps[RHO] + i] += shift;
    printf("a = %f ; shift %f \n", aa, shift);
  } else
    return (0);
  return (1);
}