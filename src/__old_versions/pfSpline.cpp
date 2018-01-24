/*
* @Author: chaomy
* @Date:   2017-11-12 02:13:44
* @Last Modified by:   chaomy
* @Last Modified time: 2017-11-12 02:13:57
*/

#include "pfHome.h"

void pfHome::splineEd(Func& func, int flag){
  double qn = 0.0;
  double un = 0.0;
  int n = func.xx.size();
  double xstep = func.step;

  const vector<double>& y = func.yy;
  vector<double>& y2 = func.g2;
  vector<double> u(n - 1, 0.0);

  const double yp1 = func.g1.front();
  if (flag == 0) func.g1.back() = 0.0;
  double ypn = func.g1.back();

  if (yp1 > 0.99e30) {
    y2[0] = 0.0;
    u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0 / (xstep)) * ((y[1] - y[0]) / (xstep)-yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    double p = 0.5 * y2[i - 1] + 2.0;
    y2[i] = (-0.5) / p;
    u[i] = (y[i + 1] - y[i]) / xstep - (y[i] - y[i - 1]) / (xstep);
    u[i] = (6.0 * u[i] / (2 * xstep) - 0.5 * u[i - 1]) / p;
  }

  if (ypn > 0.99e30) {
    qn = 0.0;
    un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 / (xstep)) * (ypn - (y[n - 1] - y[n - 2]) / (xstep));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);
  for (int k = n - 2; k >= 0; --k) y2[k] = y2[k] * y2[k + 1] + u[k];
} 

void pfHome::spline(Func& func, int flag) {
  double qn = 0.0;
  double un = 0.0;
  int n = func.xx.size();  // size

  const vector<double>& x = func.xx;
  const vector<double>& y = func.yy;
  vector<double>& y2 = func.g2;
  vector<double> u(n - 1, 0.0);

  const double yp1 = func.g1.front();
  if (flag == 0) func.g1.back() = 0.0;
  double ypn = func.g1.back();

  if (yp1 > 0.99e30) {
    y2[0] = 0.0;
    u[0] = 0.0;
  } else {
    y2[0] = -0.5;
    u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }

  for (int i = 1; i < n - 1; i++) {
    double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    double p = sig * y2[i - 1] + 2.0;
    y2[i] = (sig - 1.0) / p;
    u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
           (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
  }

  if (ypn > 0.99e30) {
    qn = 0.0;
    un = 0.0;
  } else {
    qn = 0.5;
    un = (3.0 / (x[n - 1] - x[n - 2])) *
         (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
  }

  y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

  for (int k = n - 2; k >= 0; k--) y2[k] = y2[k] * y2[k + 1] + u[k];
}

/****************************************************************
  Equal distance
 ****************************************************************/

void pfHome::splintEd(const Func& func, double r, double& val) {
  double rr = r - func.xx.front();
  /* indices into potential table */
  int k = (int)(rr * func.invstep);
  double b = (rr - k * func.step) * func.invstep;
  double a = 1.0 - b;
  double p1 = func.yy[k];
  double d21 = func.g2[k++];
  double p2 = func.yy[k];
  double d22 = func.g2[k];
  val = a * p1 + b * p2 +
        ((a * a * a - a) * d21 + (b * b * b - b) * d22) /
            (6.0 * func.invstep * func.invstep);
}

void pfHome::splintEd(const Func& func, double r, double& val, double& grad) {
  double rr = r - func.xx.front();
  int k = (int)(rr * func.invstep);
  double b = (rr - k * func.step) * func.invstep;
  double a = 1.0 - b;
  double p1 = func.yy[k];
  double d21 = func.g2[k++];
  double p2 = func.yy[k];
  double d22 = func.g2[k];
  grad = (p2 - p1) * func.invstep +
         ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) /
             (6.0 * func.invstep);
  val = a * p1 + b * p2 +
        ((a * a * a - a) * d21 + (b * b * b - b) * d22) /
            (6.0 * func.invstep * func.invstep);
}

/****************************************************************
  Non-equal distance
 ****************************************************************/
// given x value
void pfHome::splint(const Func& func, double r, double& val) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  while (hi - lo > 1) {
    int k = (hi + lo) >> 1;
    if (func.xx[k] > r)
      hi = k;
    else
      lo = k;
  }

  double h = func.xx[hi] - func.xx[lo];
  double b = (r - func.xx[lo]) / h;
  double a = (1.0 - b);

  val = a * func.yy[lo] + b * func.yy[hi] +
        ((a * a * a - a) * func.g2[lo] + (b * b * b - b) * func.g2[hi]) *
            (h * h) / 6.0;
}

// given x value
void pfHome::splint(const Func& func, double r, double& val, double& grad) {
  int lo = 0;
  int hi = func.xx.size() - 1;
  while (hi - lo > 1) {
    int k = (hi + lo) >> 1;
    if (func.xx[k] > r)
      hi = k;
    else
      lo = k;
  }

  double h = func.xx[hi] - func.xx[lo];
  double b = (r - func.xx[lo]) / h;
  double a = (1.0 - b);
  grad = (func.yy[hi] - func.yy[lo]) / h +
         ((3 * (b * b) - 1) * func.g2[hi] - (3 * (a * a) - 1) * func.g2[lo]) *
             h / 6.0;

  val = a * func.yy[lo] + b * func.yy[hi] +
        ((a * a * a - a) * func.g2[lo] + (b * b * b - b) * func.g2[hi]) *
            (h * h) / 6.0;
}

// given col and shift
void pfHome::splint(const Func& func, int k, double b, double step, double& val,
                    double& grad) {
  double a = 1.0 - b;
  double p1 = func.yy[k];
  double d21 = func.g2[k++];
  double p2 = func.yy[k];
  double d22 = func.g2[k];
  grad = (p2 - p1) / step +
         ((3 * (b * b) - 1) * d22 - (3 * (a * a) - 1) * d21) * step / 6.0;
  val = a * p1 + b * p2 +
        ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
}

// given col and shift
void pfHome::splint(const Func& func, int k, double b, double step, double& val) {
  double a = 1.0 - b;
  double p1 = func.yy[k];
  double d21 = func.g2[k++];
  double p2 = func.yy[k];
  double d22 = func.g2[k];
  val = a * p1 + b * p2 +
        ((a * a * a - a) * d21 + (b * b * b - b) * d22) * (step * step) / 6.0;
}