/*
 * @Author: chaomy
 * @Date:   2017-11-23 08:17:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-25 10:58:48
 */

#include "pfEle.h"

using std::vector;

void Melem::initGSF() {
  gsf["Nb111z110"] =
      vector<double>({0.0, 0.0055, 0.0183, 0.0292, 0.0389, 0.0442, 0.0389,
                      0.0292, 0.0183, 0.0055, 0.0});
  gsf["Nb111z211"] =
      vector<double>({0.0, 0.0065, 0.0221, 0.0368, 0.0474, 0.0515, 0.0448,
                      0.0334, 0.0214, 0.0066, 0.0000});
}