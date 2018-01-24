#ifndef _M_ELEMENTS_H
#define _M_ELEMENTS_H

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using std::string;
using std::unordered_map;
using std::vector;

class pfLMPdrv;
class pfHome;

typedef struct {
  double cij[3];
  double sij[3];
} ElasT;

typedef struct {
  double egy;
  double stsv[6];
  double strm[3][3];
} ItenT;

class Melem {
 public:
  Melem();
  void initLat();
  void initElastic();
  void initMass();
  void initDFTiten();
  void initPV();
  void initGSF();

 private:
  unordered_map<double, unordered_map<string, ItenT>> itdftm;
  unordered_map<string, double> massm;
  unordered_map<string, ElasT> cijm;
  unordered_map<string, vector<double>> pvm;
  unordered_map<string, vector<double>> gsf;
  unordered_map<string, vector<double>> lat; 
  friend class pfLMPdrv;
  friend class pfHome;
};

#endif /** _M_ELEMENTS_H **/
