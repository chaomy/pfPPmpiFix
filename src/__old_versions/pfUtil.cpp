/*
 * @Author: chaomy
 * @Date:   2017-11-04 14:53:19
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-17 08:29:04
 */

#include "pfHome.h"

using std::vector;

void pfUtil::split(const string& s, const char* delim, vector<string>& v) {
  // first duplicate the original string and return a char pointer then free the
  // memory
  char* dup = strdup(s.c_str());
  char* token = strtok(dup, delim);
  while (token != NULL) {
    v.push_back(string(token));
    // the call is treated as a subsequent calls to strtok:
    // the function continues from where it left in previous invocation
    token = strtok(NULL, delim);
  }
  free(dup);
}

void pfHome::outMkdir(string mdir){
    struct stat buf;
    if ((stat(mdir.c_str(), &buf) == 0) == 1)
        printf("kmcOutput already exists.\n");
    else
        mkdir(mdir.c_str(), S_IRWXU); 
}

void pfHome::testSpline() {
  vector<Func>& ffs = funcs;
  for (int i = 0; i < nfuncs; i++) spline(ffs[i], gradRight[i]);  // update splines

  int nums = 1000;
  double begin = 0.0, delta = 1.0 / nums;
  vector<double> xxs(nums);
  vector<double> res(nums);
  vector<double> grs(nums);

  for (int i = 0; i < 1000; i++) {
    xxs[i] = begin + i * delta;
    splint(ffs[EMF], xxs[i], res[i], grs[i]);
  }

  // output
  string fname("spl.txt.");
  FILE* fid = fopen(fname.c_str(), "w");
  for (int i = 0; i < nums; i++)
    fprintf(fid, "%f %f %f\n", xxs[i], res[i], grs[i]);
  fclose(fid);
}