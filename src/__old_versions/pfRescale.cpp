/*
 * @Author: chaomy
 * @Date:   2017-10-30 21:34:42
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-18 15:19:23
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

int pfHome::rescaleEMF(vector<double>& vv) {
  // check if rescale
  double delt = funcs[EMF].xx[1] - funcs[EMF].xx[0];
  int npts = funcs[EMF].npts;

  if (fabs(ominrho - funcs[EMF].xx.front()) > delt ||
      fabs(omaxrho - funcs[EMF].xx.back()) > delt) {
    for (int i = 0; i < npts; i++) funcs[EMF].yy[i] = vv[startps[EMF] + i];
    funcs[EMF].g1.front() = vv[startps[EMF] + npts];
    funcs[EMF].g1.back() = vv[startps[EMF] + npts + 1];

    spline(funcs[EMF], gradRight[EMF]);  // update splines

    double ndelt = (omaxrho - ominrho) / (npts - 1);

    /* update y values (to vv) */
    for (int i = 1; i < npts - 1; i++)
      splint(funcs[EMF], ominrho + i * ndelt, vv[startps[EMF] + i]);
    splint(funcs[EMF], ominrho, vv[startps[EMF] + 0], vv[startps[EMF] + npts]);
    splint(funcs[EMF], omaxrho, vv[startps[EMF] + npts - 1],
           vv[startps[EMF] + npts + 1]);

    /* update the xx values */
    for (int i = 0; i < npts; i++) funcs[EMF].xx[i] = ominrho + i * ndelt;
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