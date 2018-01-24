#ifndef pfDefines
#define pfDefines

#include <vector>
using std::string;
using std::vector;

#define DIM 3
#define MAXLEN 1024
#define PWEIGHT 100
#define EWEIGH 500
#define FRCEPS 0.1
#define PFROOT 0
#define PFERROR 1
#define PFSUCCESS 0

enum {
  PHI = 0,
  RHO = 1,
  EMF = 2,
  ADPU = 3,
  ADPW = 4,
  ADPC = 5,
  LMPPNTS = 10000
};
enum { XX = 0, YY = 1, ZZ = 2, XY = 3, YZ = 4, ZX = 5 };
enum { X = 0, Y = 1, Z = 2 };
enum { MEAMF = 3, MEAMG = 4 };
typedef double vec3[3];
typedef double position[3];
typedef double force[3];

inline double square11(double x) { return x * x; };
inline double square33(double v[3]) {
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
};

inline double innDot33(double a[3], double b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

inline double relerr(double a, double b) { return fabs(a - b) / b; }

inline void scaleVec(double v[3], double s) {
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
};

class Angle {
 private:
  double gcos;
  int slot;
  double step;
  double shift;
  double gval;
  double ggrad;
  friend class pfAtom;
  friend class pfHome;
};

class Neigh {
 private:
  int aid;
  double rhog;
  double invr;
  double r;
  double dist[3];
  double dist2r[3];
  double dis2mat[6];  // xx yy zz xy yz zx
  double uval, ugrad;
  double wval, wgrad;
  double fval, fgrad;
  vector<int> slots;
  vector<double> steps;
  vector<double> shifts;
  friend class pfAtom;
  friend class pfHome;
};

typedef struct {
  int id;
  position pst;
  force frc;
  double absfrc;
} AtomMPI;

class pfAtom {
 private:
  int id;
  int nneighs;
  int nneighsFull;
  position pst;
  force frc;
  force fitfrc;
  force fweigh;

  /* ADP */
  double mu[3];
  double nu;
  double lambda[6];

  vector<Neigh> neighs;

  /* MEAM */
  vector<Neigh> neighsFull;
  vector<vector<Angle>> angMat;

  double absfrc;
  double rho;
  double gradF;

 public:
  pfAtom() : id(0) {
    fitfrc[0] = 0.0;
    fitfrc[1] = 0.0;
    fitfrc[2] = 0.0;
  };
  pfAtom(int n) : id(n) {
    fitfrc[0] = 0.0;
    fitfrc[1] = 0.0;
    fitfrc[2] = 0.0;
  };
  ~pfAtom(){};
  friend class Config;
  friend class pfHome;
};

class Config {
 private:
  int natoms;
  int cfgid;
  double weigh;
  double engy;
  double fitengy;
  double vol;
  int scale[3];
  vec3 bvx, tvx;
  vec3 bvy, tvy;
  vec3 bvz, tvz;

  vector<pfAtom> atoms;

 public:
  Config() { cfgid = 0, natoms = 0, weigh = 0.0, engy = 0.0; };
  Config(int n) : cfgid(n) { natoms = 0, weigh = 0.0, engy = 0.0; };
  ~Config(){};
  friend class pfHome;
};

class Func {
 private:
  int npts;
  int bnd;
  double step;     // equidist; has't initilaized yet
  double invstep;  // equidist
  vector<double> xx;
  vector<double> yy;
  vector<double> g1;
  vector<double> g2;
  friend class pfHome;
};

#endif