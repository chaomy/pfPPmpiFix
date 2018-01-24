#ifndef pfDefines
#define pfDefines

#include <boost/mpi/collectives/broadcast.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/datatype_fwd.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpl/and.hpp>
#include <boost/serialization/export.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/shared_ptr.hpp>
#include <functional>
#include <numeric>
#include <vector>
#include "spline.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

#define EXT 0
#define DIM 3
#define MAXLEN 1024
#define PWEIGHT 100
#define EWEIGH 500
#define FRCEPS 0.1
#define PI 3.14159265359
#define INVPI 0.31830988618
#define PFROOT 0
#define PFERROR 1
#define PFSUCCESS 0
#define LMPPNTS 10000
#define EVA3_GPA 160.21766208

// nearest neighbors
// 1st 8 ; 2nd 6; 3rd 12; 4th 24; 5th 8
enum { PHI = 0, RHO = 1, EMF = 2, ADPU = 3, ADPW = 4, ADPC = 5 };
enum { N1 = 8, N2 = 14, N3 = 26, N4 = 50, N5 = 58 };
enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };
enum { MEAMF = 3, MEAMG = 4 };

inline void printVec(double v[3]) {
  for (int i = 0; i < 3; i++) cout << v[i] << " ";
  cout << endl;
}
inline double square11(double x) { return x * x; };
inline double square33(const vector<double> &v) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
};
inline double innDot33(const vector<double> &a, const vector<double> b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};
inline double relerr(double a, double b) { return fabs(a - b) / b; }
inline void scaleVec(vector<double> &v, double s) {
  for (auto &ee : v) ee *= s;
};
inline void crossProd33(const vector<double> &a, const vector<double> &b,
                        vector<double> &res) {
  res[0] = a[1] * b[2] - a[2] * b[1];
  res[1] = a[2] * b[0] - a[0] * b[2];
  res[2] = a[0] * b[1] - a[1] * b[0];
}
inline double vecInnProd33(const vector<double> &a, const vector<double> &b) {
  double res = 0;
  for (int i = 0; i < 3; i++) res += a[i] * b[i];
  return res;
}

class Angle {
 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &gcos &gval &ggrad &step &shift &slot;
  }
  double gcos, gval, ggrad, step, shift;
  int slot;
  friend class pfAtom;
  friend class pfHome;
};

class Neigh {
 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &aid &r &invr &psum &phi &phig &rho &rhog &uval &ugrad &wval &wgrad &fval
        &fgrad &dist &dist2r &dis2mat &slots &steps &shifts;
  }

  int aid;
  double r, invr, psum;
  double phi, phig, rho, rhog, uval, ugrad, wval, wgrad, fval, fgrad;
  vector<double> dist, dist2r, dis2mat;  // xx yy zz xy yz zx
  vector<int> slots;
  vector<double> steps, shifts;

 public:
  Neigh() : dist(3), dist2r(3), dis2mat(6) {}
  friend class pfAtom;
  friend class pfHome;
};

class pfAtom {
 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &id &tp &nneighs &nneighsFull &rho &crho &prho &gradF &absfrc &pst &prl
        &frc &fitfrc &fweigh &phifrc &rhofrc &trifrc &mu &nu &lambda &neighs
            &neighsFull &angMat;
  }

  int id, tp, nneighs, nneighsFull;
  double rho, crho, prho, gradF, absfrc;

  vector<double> pst, prl;
  vector<double> frc, fitfrc, fweigh, phifrc, rhofrc, trifrc;

  /* ADP */
  vector<double> mu;
  double nu;
  vector<double> lambda;
  vector<Neigh> neighs;

  /* MEAM */
  vector<Neigh> neighsFull;
  vector<vector<Angle>> angMat;

 public:
  pfAtom()
      : id(0),
        pst(3),
        prl(3),
        frc(3),
        fitfrc(3),
        fweigh(3),
        phifrc(3),
        rhofrc(3),
        trifrc(3),
        mu(3),
        lambda(6){};
  pfAtom(int n)
      : id(n),
        pst(3),
        prl(3),
        frc(3),
        fitfrc(3),
        fweigh(3),
        phifrc(3),
        rhofrc(3),
        trifrc(3),
        mu(3),
        lambda(6){};
  ~pfAtom(){};
  friend class Config;
  friend class pfHome;
};

class Config {
 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &cfgid &natoms &weigh &engy &fitengy &phiengy &emfengy &vol &rhomx &rhomi
        &scale &bvx &tvx &bvy &tvy &bvz &tvz &strs &atoms &natomsv &nelemsv;
  }

  int cfgid, natoms;
  double weigh, engy;
  double fitengy, phiengy, emfengy;
  double vol, rhomx, rhomi;
  vector<int> scale;
  vector<double> bvx, tvx, bvy, tvy, bvz, tvz;
  vector<double> strs;
  vector<pfAtom> atoms;
  vector<int> natomsv;
  vector<string> nelemsv;

 public:
  Config()
      : cfgid(0),
        natoms(0),
        scale(3),
        bvx(3),
        tvx(3),
        bvy(3),
        tvy(3),
        bvz(3),
        tvz(3),
        strs(6) {
    weigh = 0.0, engy = 0.0;
  };
  Config(int n)
      : cfgid(n),
        natoms(0),
        scale(3),
        bvx(3),
        tvx(3),
        bvy(3),
        tvy(3),
        bvz(3),
        tvz(3),
        strs(6) {
    weigh = 0.0, engy = 0.0;
  };
  ~Config(){};
  friend class pfHome;
};

class Func {
 private:
  int npts, bnd;
  double rng;      // range ymax - ymin
  double step;     // equidist; has't initilaized yet
  double invstep;  // equidist
  vector<double> xx, yy, g1, g2;
  tk::spline s;
  friend class pfHome;
};

// namespace boost {
// namespace mpi {
// template <>
// struct is_mpi_datatype<Angle>
//     : public mpl::and_<is_mpi_datatype<int>, is_mpi_datatype<float>> {};
// }  // namespace mpi
// }  // namespace boost

#endif