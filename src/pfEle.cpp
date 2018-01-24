/*
 * @Author: chaomy
 * @Date:   2017-11-22 14:08:35
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-12-16 21:10:56
 */

#include "pfEle.h"

/* -----------------------------------------------------------------------------
 * Atomic weight data taken from:
 * Pure Appl. Chem., Vol. 83, No. 2, pp. 359–396, 2011.
 * Atomic weights of the elements 2009 (IUPAC Technical Report)
 * ---------------------------------------------------------------------------*/
const int numMax = 112;
const double weight[] = {
    0.,                                                      // Null
    1.00797,   4.0026,     6.939,      9.012182,    10.811,  // H - B
    12.01115,  14.0067,    15.9994,    18.9984032,  20.17976,
    22.989769,                                                  // C - Na
    24.30506,  26.9815386, 28.086,     30.973762,   32.064,     // Mg - S
    35.453,    39.948,     39.0983,    40.078,      44.955912,  // Cl - Sc
    47.867,    50.9415,    51.9961,    54.938045,   55.845,     // Ti - Fe
    58.933195, 58.6934,    63.546,     65.38,       69.723,     // Co - Ga
    72.63,     74.9216,    78.96,      79.904,      83.798,     // Ge - Kr
    85.4678,   87.62,      88.90585,   91.224,      92.90638,   // Rb - Nb
    95.96,     98.9062,    101.07,     102.9055,    106.42,     // Mo - Pd
    107.8682,  112.411,    114.818,    118.71,      121.76,     // Ag - Sb
    127.6,     126.90447,  131.293,    132.9054519, 137.327,    // Te - Ba
    138.90547, 140.116,    140.90765,  144.242,     147.,       // La - Pm
    150.36,    151.964,    157.25,     158.92535,   162.5,      // Sm - Dy
    164.93032, 167.259,    168.93421,  173.054,     174.9668,   // Ho - Lu
    178.49,    180.94788,  183.84,     186.207,     190.23,     // Hf - Os
    192.217,   195.084,    196.966569, 200.59,      204.383,    // Ir - Tl
    207.2,     208.9804,   209.,       210.,        222.,       // Pb - Rn
    223.,      226.025,    227.028,    232.03806,   231.03588,  // Fr - Pa
    238.02891, 237.048,    244.,       243.,        247.,       // U - Cm
    247.,      251.,       252.,       257.,        258.,       // Bk - Md
    259.,      260.,       261.11,     262.11,      263.12,     // No - Sg
    262.12,    265.,       266.,       269.,        272.,
    285.  // Bh - Cn
};

const char symbol[][3] = {
    "X",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na",
    "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",
    "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
    "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag",
    "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi",
    "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am",
    "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh",
    "Hs", "Mt", "Ds", "Rg", "Cn"};

Melem::Melem() {
  initLat();
  initElastic();
  initMass();
  initDFTiten();
  initPV();
  initGSF();
}

void Melem::initMass() {
  for (int i = 0; i < numMax; i++) massm[string(symbol[i])] = weight[i];
}

void Melem::initElastic() {
  ElasT elm;
  elm.cij[0] = 124.000000;
  elm.cij[1] = 93.400000;
  elm.cij[2] = 46.100000;
  elm.sij[0] = 0.022900;
  elm.sij[1] = -0.009800;
  elm.sij[2] = 0.021700;
  cijm["Ag"] = elm;

  elm.cij[0] = 107.300000;
  elm.cij[1] = 60.900000;
  elm.cij[2] = 28.300000;
  elm.sij[0] = 0.015800;
  elm.sij[1] = -0.005700;
  elm.sij[2] = 0.035300;
  cijm["Al"] = elm;

  elm.cij[0] = 192.900000;
  elm.cij[1] = 163.800000;
  elm.cij[2] = 41.500000;
  elm.sij[0] = 0.023500;
  elm.sij[1] = -0.010800;
  elm.sij[2] = 0.024100;
  cijm["Au"] = elm;

  elm.cij[0] = 168.400000;
  elm.cij[1] = 121.400000;
  elm.cij[2] = 75.400000;
  elm.sij[0] = 0.015000;
  elm.sij[1] = -0.006300;
  elm.sij[2] = 0.013300;
  cijm["Cu"] = elm;

  elm.cij[0] = 580.000000;
  elm.cij[1] = 242.000000;
  elm.cij[2] = 256.000000;
  elm.sij[0] = 0.002300;
  elm.sij[1] = -0.000700;
  elm.sij[2] = 0.003900;
  cijm["Ir"] = elm;

  elm.cij[0] = 246.500000;
  elm.cij[1] = 147.300000;
  elm.cij[2] = 127.400000;
  elm.sij[0] = 0.007300;
  elm.sij[1] = -0.002700;
  elm.sij[2] = 0.007800;
  cijm["Ni"] = elm;

  elm.cij[0] = 49.500000;
  elm.cij[1] = 42.300000;
  elm.cij[2] = 14.900000;
  elm.sij[0] = 0.095100;
  elm.sij[1] = -0.043800;
  elm.sij[2] = 0.067100;
  cijm["Pb"] = elm;

  elm.cij[0] = 227.100000;
  elm.cij[1] = 176.000000;
  elm.cij[2] = 71.700000;
  elm.sij[0] = 0.013600;
  elm.sij[1] = -0.005900;
  elm.sij[2] = 0.013900;
  cijm["Pd"] = elm;

  elm.cij[0] = 346.700000;
  elm.cij[1] = 250.700000;
  elm.cij[2] = 76.500000;
  elm.sij[0] = 0.007300;
  elm.sij[1] = -0.003100;
  elm.sij[2] = 0.013100;
  cijm["Pt"] = elm;

  elm.cij[0] = 339.800000;
  elm.cij[1] = 58.600000;
  elm.cij[2] = 99.000000;
  elm.sij[0] = 0.003100;
  elm.sij[1] = -0.000500;
  elm.sij[2] = 0.010100;
  cijm["Cr"] = elm;

  elm.cij[0] = 231.400000;
  elm.cij[1] = 134.700000;
  elm.cij[2] = 116.400000;
  elm.sij[0] = 0.007600;
  elm.sij[1] = -0.002800;
  elm.sij[2] = 0.008600;
  cijm["Fe"] = elm;

  elm.cij[0] = 4.140000;
  elm.cij[1] = 2.210000;
  elm.cij[2] = 2.630000;
  elm.sij[0] = 0.384400;
  elm.sij[1] = -0.133800;
  elm.sij[2] = 0.380200;
  cijm["K"] = elm;

  elm.cij[0] = 13.500000;
  elm.cij[1] = 11.440000;
  elm.cij[2] = 8.780000;
  elm.sij[0] = 0.332800;
  elm.sij[1] = -0.152600;
  elm.sij[2] = 0.113900;
  cijm["Li"] = elm;

  elm.cij[0] = 440.800000;
  elm.cij[1] = 172.400000;
  elm.cij[2] = 121.700000;
  elm.sij[0] = 0.002900;
  elm.sij[1] = -0.000800;
  elm.sij[2] = 0.008200;
  cijm["Mo"] = elm;

  elm.cij[0] = 6.150000;
  elm.cij[1] = 4.960000;
  elm.cij[2] = 5.920000;
  elm.sij[0] = 0.581000;
  elm.sij[1] = -0.259400;
  elm.sij[2] = 0.168900;
  cijm["Na"] = elm;

  elm.cij[0] = 240.200000;
  elm.cij[1] = 125.600000;
  elm.cij[2] = 28.200000;
  elm.sij[0] = 0.006500;
  elm.sij[1] = -0.002200;
  elm.sij[2] = 0.035500;
  cijm["Nb"] = elm;

  elm.cij[0] = 260.200000;
  elm.cij[1] = 154.500000;
  elm.cij[2] = 82.600000;
  elm.sij[0] = 0.006900;
  elm.sij[1] = -0.002600;
  elm.sij[2] = 0.012100;
  cijm["Ta"] = elm;

  elm.cij[0] = 228.000000;
  elm.cij[1] = 118.700000;
  elm.cij[2] = 42.600000;
  elm.sij[0] = 0.006800;
  elm.sij[1] = -0.002300;
  elm.sij[2] = 0.023500;
  cijm["V"] = elm;

  elm.cij[0] = 522.400000;
  elm.cij[1] = 204.400000;
  elm.cij[2] = 160.800000;
  elm.sij[0] = 0.002500;
  elm.sij[1] = -0.000700;
  elm.sij[2] = 0.006200;
  cijm["W"] = elm;

  elm.cij[0] = 949.000000;
  elm.cij[1] = 151.000000;
  elm.cij[2] = 521.000000;
  elm.sij[0] = 0.001100;
  elm.sij[1] = -0.000200;
  elm.sij[2] = 0.001900;
  cijm["C"] = elm;

  elm.cij[0] = 128.400000;
  elm.cij[1] = 48.200000;
  elm.cij[2] = 66.700000;
  elm.sij[0] = 0.009800;
  elm.sij[1] = -0.002700;
  elm.sij[2] = 0.015000;
  cijm["Ge"] = elm;

  elm.cij[0] = 166.200000;
  elm.cij[1] = 64.400000;
  elm.cij[2] = 79.800000;
  elm.sij[0] = 0.007700;
  elm.sij[1] = -0.002100;
  elm.sij[2] = 0.012500;
  cijm["Si"] = elm;

  elm.cij[0] = 118.800000;
  elm.cij[1] = 53.700000;
  elm.cij[2] = 59.400000;
  elm.sij[0] = 0.011700;
  elm.sij[1] = -0.003600;
  elm.sij[2] = 0.016800;
  cijm["GaAs"] = elm;

  elm.cij[0] = 141.200000;
  elm.cij[1] = 62.500000;
  elm.cij[2] = 70.500000;
  elm.sij[0] = 0.009700;
  elm.sij[1] = -0.003000;
  elm.sij[2] = 0.014200;
  cijm["GaP"] = elm;

  elm.cij[0] = 102.200000;
  elm.cij[1] = 57.600000;
  elm.cij[2] = 46.000000;
  elm.sij[0] = 0.016500;
  elm.sij[1] = -0.005900;
  elm.sij[2] = 0.021700;
  cijm["InP"] = elm;

  elm.cij[0] = 39.500000;
  elm.cij[1] = 4.900000;
  elm.cij[2] = 6.300000;
  elm.sij[0] = 0.026000;
  elm.sij[1] = -0.002900;
  elm.sij[2] = 0.158700;
  cijm["KCl"] = elm;

  elm.cij[0] = 114.000000;
  elm.cij[1] = 47.700000;
  elm.cij[2] = 63.600000;
  elm.sij[0] = 0.011600;
  elm.sij[1] = -0.003400;
  elm.sij[2] = 0.015700;
  cijm["LiF"] = elm;

  elm.cij[0] = 287.600000;
  elm.cij[1] = 87.400000;
  elm.cij[2] = 151.400000;
  elm.sij[0] = 0.004100;
  elm.sij[1] = -0.000900;
  elm.sij[2] = 0.006600;
  cijm["MgO"] = elm;

  elm.cij[0] = 49.600000;
  elm.cij[1] = 12.400000;
  elm.cij[2] = 12.900000;
  elm.sij[0] = 0.022400;
  elm.sij[1] = -0.004500;
  elm.sij[2] = 0.077500;
  cijm["NaCl"] = elm;

  elm.cij[0] = 500.000000;
  elm.cij[1] = 113.000000;
  elm.cij[2] = 175.000000;
  elm.sij[0] = 0.002200;
  elm.sij[1] = -0.000400;
  elm.sij[2] = 0.005700;
  cijm["TiC"] = elm;
}
