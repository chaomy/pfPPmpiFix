#ifndef pfHome_H_
#define pfHome_H_

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <unordered_map>
#include <sys/stat.h>
#include <sys/types.h>
#include "lmpMatrix.h"
#include "nlopt.hpp"
#include "pfDefines.h"

using std::string;
using std::unordered_map;
using std::vector;
using namespace MMatrix;

class pfLMPdrv;
class pfOptimizer;

class pfHome {
 private:
  int gcnt;  // count the times of calling force calculation
  int scnt;  // count the times of calling simulated annealing 
  int nvars;
  int nfuncs;
  int nconfs;
  int locstt;
  int locend;
  double ricut;
  double rocut;
  double lorho;
  double hirho;
  double punish;
  double physic;

  double ominrho;
  double omaxrho;
  double oaverho;

  /* parameters */
  unordered_map<string, double> dparams;
  unordered_map<string, int> iparams;
  unordered_map<string, string> sparams;

  vector<Config> configs;
  vector<Func> funcs;
  vector<int> startps;
  vector<int> gradRight;
  pfLMPdrv* lmpdrv; 
  pfOptimizer* optdrv;

  vector<double> hil; // 5 + 1
  vector<double> lol; // 5 + 1
  vector<double> recorderr;
  vector<double> ini;
  vector<double> lob;
  vector<double> hib;

 public:
  pfHome(int argc, char* argv[]);
  ~pfHome();

  unordered_map<string, string> gsparams() const { return sparams; };
  unordered_map<string, double> gdparams() const { return dparams; };
  unordered_map<string, int> giparams() const { return iparams; };

  // initialization
  void pfInit();
  void parseArgs(int argc, char* argv[]);
  void readConfig();
  void readPot();
  void readParam();
  void initParam();
  void writePot();
  void writePot(const vector<double>& vv);
  void writeLMPS();
  void writeLMPS(const vector<double>& vv);
  void writeMEAM();
  void initBox();
  void initNeighs();
  void initNeighsFull(); /* meam */
  void initAngles();
  void setNeighslot(Neigh& n, Func f, double r);
  void updateNeighslot(Neigh& n, Func f, double r, int id);
  void setAngleslot(Angle& a, Func f, double r);

  // spline interpolation
  void spline(Func& func, int flag);
  void splineEd(Func& func, int flag);
  void splintEd(const Func& func, double r, double& val);
  void splintEd(const Func& func, double r, double& val, double& grad);
  void splint(const Func& func, double r, double& val, double& grad);
  void splint(const Func& func, double r, double& val);
  void splint(const Func& func, int k, double b, double step, double& val);
  void splint(const Func& func, int k, double b, double step, double& val,
              double& grad);

  // force calculation
  void doShift();  // make the rho / emf good
  void run();
  void run(int argc, char* argv[]);
  void calErr();
  void updaterho(vector<double>& vv);
  void updaterhoMEAM(vector<double>& vv);
  double errFunct(const vector<double>& x);
  double errFunctGrad(const vector<double>& x, vector<double>& g);
  double forceEAM(vector<Func>& ffs, int tag);
  double forceEAM(const vector<double>& vv, int tag);
  double forceADP(const vector<double>& vv, int tag);
  double forceMEAM(const vector<double>& vv, int tag);

  // optimizatio  
  void simAnneal();
  void randomize(vector<double>& vv, const int n, const vector<double>& v);
  int rescaleEMF(vector<double>& vv);
  int rescaleRHO(vector<double>& vv);
  void shiftRHO(vector<double>& vv);
  void shiftEMF(double shift);
  void nloptGlobal();

  // random
  double randNormal();
  double randUniform();
  double randUniform(const double min, const double max);

  // increase nodes 
  void upgrade(int id);
  void increAnneal(); 
  void recordStage(int cnt);

  // MPI utilitis
  void pfMPIinit(int argc, char* argv[]);
  void pfMPIfinalize();
  int createMPIdataTypes();

  // utils 
  void outMkdir(string mdir);

  // debug
  void testSpline();
  friend class pfLMPdrv;
  friend class pfOptimizer;
};

class pfUtil {
 public:
  void split(const string& s, const char* delim, vector<string>& v);
  friend class pfHome;
};

#endif  // pfHome_H_