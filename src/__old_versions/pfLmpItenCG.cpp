/*
 * @Author: chaomy
 * @Date:   2017-11-10 14:44:56
 * @Last Modified by:   chaomy
 * @Last Modified time: 2017-11-22 11:21:27
 */

#include "lmpStruct.h"
#include "pfLMPdrv.h"

#define SQRT2 1.4142135623730951
#define MAXLINE_SEARCH 200
#define MAXBFGS_STEP 1000
#define PERTERB 0.2
#define STRESS_SMALL 0.05
#define STRAIN_WEIGHT_COE 2e3
// 0.070710678118654766 [>**  Stress threshold 50 Mpa (0.05 Gpa) sqrt(2) * 0.05
// **<]

#define EPSILON (2 * STRESS_SMALL * STRESS_SMALL)
#define VARN 2
#define SMALL_STEP 0.00005 /*** for calculate derivative ***/
#define INV_SMALL_STEP (1. / SMALL_STEP)
#define vect2(a, b) (sqrt(a * a + b * b))
#define dim2Dot(vec1, vec2) ((vec1[0]) * (vec2[0]) + (vec1[1]) * (vec2[1]))
#define EPSILON0 1e-9

/*** calculate dim2 vector out product, the result is mtx[2][2] ***/
inline void dim2outprod(double va[2], double vb[2], double mtx[2][2]) {
  mtx[0][0] = va[0] * vb[0];
  mtx[0][1] = va[0] * vb[1];
  mtx[1][0] = va[1] * vb[0];
  mtx[1][1] = va[1] * vb[1];
  return;
}

void pfLMPdrv::loopcalIten() {
  char tagT[] = "tpath";
  char tagO[] = "opath";

  errItenTE = errItenTS = 0.0;
  errItenOE = errItenOS = 0.0;

  if (mycalOpt.itentag == 1) {
    for (int i = 0; i < mycalOpt.itenn; i++) {
      /*** calculate tpath ***/
      calIdealtenLinear(itenstrain[i], tagT, i);

      itentpstss[i] = stressVec[0];
      tpstrainyy[i] = gstrainMtx[1][1];
      tpstrainzz[i] = gstrainMtx[2][2];

      errItenTE += itenEngyErrList[i];
      errItenTS += itenStssErrList[i];
      errItenTStr += itenStrErrList[i];

      /*** calculate opath ***/
      calIdealtenLinear(itenstrain[i], tagO, i);

      itenopstss[i] = stressVec[0];
      opstrainyy[i] = gstrainMtx[1][1];
      opstrainzz[i] = gstrainMtx[2][2];

      errItenOE += itenEngyErrList[i + TPOINTS];
      errItenOS += itenStssErrList[i + TPOINTS];
      errItenOStr += itenStrErrList[i + TPOINTS];
    }
    outputStressInfoSim("tensile.log", "w", "ideal");
  }
  if (mycalOpt.itentag == 2) {
    for (int i = 0; i < mycalOpt.itenn; i++) {
      /*** calculate tpath ***/
      calIdealtenLinear(itenstrain[i], tagT, i);
      itentpstss[i] = stressVec[0];
      tpstrainyy[i] = gstrainMtx[1][1];
      tpstrainzz[i] = gstrainMtx[2][2];

      /*** calculate opath ***/
      calIdealtenLinear(itenstrain[i], tagO, i);
      itenopstss[i] = stressVec[0];
      opstrainyy[i] = gstrainMtx[1][1];
      opstrainzz[i] = gstrainMtx[2][2];
    }
    outputStressInfoSim("tensile.log", "w", "ideal");
  }
}

/****************************************************************
 * in LAMMPS, use self-defined a1 a2 a3 in lattice command
 * to add the lattice strain (including ideal tensile and shear)
 ****************************************************************/
void lmpDriver::calIdealten(double delta, char tag[]) {
  /*** declare and define local variables  ***/
  double correctStrain[DIM][DIM];

  // double elasticSij[STRESS_N][STRESS_N];
  double tempMtx[DIM][DIM];
  double perfBase[DIM][DIM];

  Melem myEle;
  myEle.initElastic();

  Mmatrix<double>::initMtx33(correctStrain);
  Mmatrix<double>::initMtx33(strainMtx);
  Mmatrix<double>::initMtx33(baseVect);
  Mmatrix<double>::initMtx33(tempMtx);

  /** initialize the base vector  **/
  if (strcasecmp(tag, "tpath") == 0) {
    baseVect[0][0] = 1.0;
    baseVect[1][1] = 1.0;
    baseVect[2][2] = 1.0;

    /*** store the pefect Base_vector ***/
    Mmatrix<double>::mtx33copymtx33(perfBase, baseVect);

  } else if (strcasecmp(tag, "opath") == 0) {
    baseVect[0][0] = 1.0;
    baseVect[1][1] = SQRT2;
    baseVect[2][2] = SQRT2;

    // [>** store the pefect Base_vector **<]
    Mmatrix<double>::mtx33copymtx33(perfBase, baseVect);
  }

  /*** loop of ideal tensile ***/
  int count = 0;
  FILE* ptFile;

  /***************************************************************
   * defines parameters for optimization with BFGS method
   * **************************************************************/
  double x[VARN]; /** store inital strain[1][1] strain[2][2] **/
  double s[VARN]; /*** for update H (for secant equation)  ***/
  double y[VARN]; /*** for update H (for secant equation)  ***/
  double rhob;    /*** rho_k = 1./dot(yk,sk) ***/

  /*** line search param ***/
  double alpha[VARN];      /** step length **/
  double p[VARN];          /** trial direction **/
  double g[VARN];          /** graidiant at xk **/
  double gn[VARN];         /** graidiant at xk+1 **/
  int cnt[VARN];           /** cnt the times of finding step length **/
  double last_alpha[VARN]; /** store the previous stepLength **/

  /*** alhpa0 = I in Newton and quasi-Newton methods ***/
  last_alpha[0] = last_alpha[1] = 1;
  cnt[0] = cnt[1] = 0;

  int id = 0;
  double rho = 0.6;
  double c = 0.4;  // c should < 0.5
  double curval;
  double temp = 0;

  double H[VARN][VARN];  // H is B^-1  (inverse of Hessian approximations B)
  double Tmp[VARN][VARN];

  /***************************************************************
   *  use BFGS method to find the optimial strain
   * **************************************************************/

  /** initialize strain  **/
  strainMtx[0][0] = 1 + delta;
  strainMtx[1][1] = 1 + PERTERB * delta;  // some prterbation
  strainMtx[2][2] = 1 - PERTERB * delta;  // some perterbation

  /** add strain on baseVect **/
  Mmatrix<double>::mtx33copymtx33(tempMtx, baseVect);
  Mmatrix<double>::mtx33Multmtx33(strainMtx, tempMtx, baseVect);

  /*** initialize the inverse Hessian approximation H0 ***/
  // try unit matrix with coeff
  double coeff = 1. / 1000.;  // tuned parameter

  H[0][0] = H[1][1] = 1.0 * coeff;
  H[0][1] = H[1][0] = 0.0;

  // initialize  x
  x[0] = strainMtx[1][1];
  x[1] = strainMtx[2][2];

  /*** update strain ***/
  strainMtx[1][1] = x[0];
  strainMtx[2][2] = x[1];

  /** add the strain on baseVect **/
  Mmatrix<double>::mtx33Multmtx33(strainMtx, perfBase, baseVect);

  /*** run lammps ***/
  idealLmpRun(tag);

  /*** currently curval value ***/
  curval = vect2(stressVec[1], stressVec[2]);

  /*** print the info ***/
  if ((ptFile = fopen("relax.log", "w")) == NULL) {
    fprintf(stderr, "error when open relax.log. %s:%d", __FILE__, __LINE__);
    exit(1);
  };

  fprintf(ptFile, "target %f\n", curval);
  fprintf(ptFile, "[%9.8f %9.8f %9.8f] ", strainMtx[0][0], strainMtx[1][1],
          strainMtx[2][2]);
  fprintf(ptFile, "[%9.8f %9.8f %9.8f]\n", stressVec[0], stressVec[1],
          stressVec[2]);

  fflush(ptFile);
  fclose(ptFile);

  H[0][0] = -coeff * stressVec[1] /
            fabs(stressVec[1]);  // negative stress -> decrease the strain
  H[1][1] = -coeff * stressVec[2] /
            fabs(stressVec[2]);  // positive stress -> increase the strain

  while (true) {
    if (curval < EPSILON)
      break; /** Syy and Szz are small ***/
    else if (count > MAXBFGS_STEP) {
      printf(" run too much times, tune your parameters \n");
      break;
    } else {
      /* backtracking line search method */
      temp = curval;

      /*** first calculate derivative g[i] at this point ***/
      for (int i = 0; i < VARN; i++) {
        id = i + 1;

        // update the strain
        strainMtx[id][id] += SMALL_STEP;

        /** add the strain on baseVect **/
        Mmatrix<double>::mtx33Multmtx33(strainMtx, perfBase, baseVect);

        /** run **/
        idealLmpRun(tag);

        /** cal graidiant **/
        g[i] = INV_SMALL_STEP * (vect2(stressVec[1], stressVec[2]) - temp);

        /** recover strainMtx **/
        strainMtx[id][id] -= SMALL_STEP;

      } /*** loop VARN ***/

      /*** update the direction vector p ***/

      /*** pk = H * gv ***/
      p[0] = -(H[0][0] * g[0] + H[0][1] * g[1]);
      p[1] = -(H[1][0] * g[0] + H[1][1] * g[1]);

      double gdotp = (g[0] * p[0] + g[1] * p[1]);  // graident dot direction

      /****************************************************************
       * Backtracking Line Search to find step length
       ****************************************************************/
      alpha[0] = 1.;
      alpha[1] = 1.;

      for (int i = 0; i < VARN; i++) {  // i = 0, 1;
        cnt[i] = 0;

        while (true) {
          /*** i = 0; change matx[1][1] ; i = 1; change mtx[2][2];
           * ***/
          strainMtx[i + 1][i + 1] = x[i] + alpha[i] * p[i];

          /** add the strain on baseVect **/
          Mmatrix<double>::mtx33Multmtx33(strainMtx, perfBase, baseVect);

          if ((fabs(baseVect[1][1]) > 3.0) || (fabs(baseVect[2][2]) > 3.0)) {
            printf("b[1][1] = %f  b[2][2] = %f", baseVect[1][1],
                   baseVect[2][2]);
            exit(LMP_ERR);
          }

          if (myRank == LmpRoot)
            Mmatrix<double>::printMtx3x3("strain.log", "base", strainMtx);
          if (myRank == LmpRoot)
            Mmatrix<double>::printMtx3x3("stress.log", "base", baseVect);

          /** run lammps **/
          idealLmpRun(tag);

          /*** check the wolfe conditions ***/
          double left = vect2(stressVec[1], stressVec[2]);  // f(xk + alpha pk)
          double right =
              curval + c * alpha[i] * gdotp;  // f(xk) + c alpha fk' pk

          /*** f(xk + alpha pk) <= f(xk) + c alpha fk' pk ***/
          if (left <= right) break;

          alpha[i] *= rho;
          cnt[i] += 1;

          if (cnt[i] > MAXLINE_SEARCH) {
            printf("too many steps in line search %s:%d", __FILE__, __LINE__);
            break;
          } /***  end if ***/
        }   /*** end while ***/
        /** recoever the strainMtx after finding alpha **/
        strainMtx[i + 1][i + 1] = x[i];

      } /*** loop VARN ***/

      /*** store last_alpha ***/
      last_alpha[0] = alpha[0];
      last_alpha[1] = alpha[1];

      printf("The cnt[0] = %d cnt[1] = %d \n", cnt[0], cnt[1]);
      printf("The alpha[0] = %f; alpha[1] = %f\n", alpha[0], alpha[1]);

      /*** update the x ***/
      s[0] = alpha[0] * p[0];
      s[1] = alpha[1] * p[1];

      x[0] = x[0] + s[0];
      x[1] = x[1] + s[1];

      /*** calculate derivative gn[i] at new x ***/

      /*** update the strain ***/
      strainMtx[1][1] = x[0];
      strainMtx[2][2] = x[1];

      /** add the strain on baseVect **/
      Mmatrix<double>::mtx33Multmtx33(strainMtx, perfBase, baseVect);

      /** run **/
      idealLmpRun(tag);

      /*** print the info ***/
      if (myRank == LmpRoot) {
        /* Mmatrix<double>::printMtx3x3("strain.log",
         "strain",strainMtx);
         Mmatrix<double>::printMtx3x3("stress.log", "base", baseVect);
         Mmatrix<double>::printVec6(  "stress.log",
         "stress",stressVec);
        */

        if ((ptFile = fopen("relax.log", "a")) == NULL) {
          fprintf(stderr, "error when open relax.log. %s:%d", __FILE__,
                  __LINE__);
          exit(1);
        };

        fprintf(ptFile, "%d -> %f\n", count, curval);
        fprintf(ptFile, "alpha[0]=%9.8f, alpha[1]=%9.8f ", alpha[0], alpha[1]);
        fprintf(ptFile, "p[0]= %9.8f, p[1]=%9.8f\n", p[0], p[1]);
        fprintf(ptFile, "[%9.8f %9.8f %9.8f]  ", strainMtx[0][0],
                strainMtx[1][1], strainMtx[2][2]);
        fprintf(ptFile, "[%9.8f %9.8f %9.8f]\n", stressVec[0], stressVec[1],
                stressVec[2]);

        fprintf(ptFile, "g[0] = %f, g[1] = %f\n", g[0], g[1]);
        fprintf(ptFile, "gn[0]= %f, gn[1]= %f\n", gn[0], gn[1]);

        fflush(ptFile);
        fclose(ptFile);
      }

      curval = temp =
          vect2(stressVec[1], stressVec[2]);  // function value at xk+1;

      /** judge **/
      if (curval < EPSILON) {
        break;  // jump out of while
      }

      for (int i = 0; i < VARN; i++) {
        id = i + 1;

        /*** update the strain ***/
        strainMtx[id][id] += SMALL_STEP;

        /** add the strain on baseVect **/
        Mmatrix<double>::mtx33Multmtx33(strainMtx, perfBase, baseVect);

        if ((fabs(baseVect[1][1]) > 3.0) || (fabs(baseVect[2][2]) > 3.0)) {
          printf("b[1][1] = %f  b[2][2] = %f", baseVect[1][1], baseVect[2][2]);
          exit(LMP_ERR);
        }

        /** run **/
        idealLmpRun(tag);

        /** cal graidiant **/
        gn[i] = INV_SMALL_STEP * (vect2(stressVec[1], stressVec[2]) - temp);

        /** recover strainMtx **/
        strainMtx[id][id] -= SMALL_STEP;

        /*** update y (yk = f'k+1 - f'k) ***/
        y[i] = gn[i] - g[i];

      } /*** loop VARN ***/

      /*******************************************************************
       * update H
       *******************************************************************/

      /*** update rhob = 1./(yk * s[k]) ***/
      rhob = 1. / (dim2Dot(y, s));

      /*** ( I - rho s y ) H ( I - rho y s ) + rho s s ***/
      double HL[2][2];
      double HR[2][2];
      double LF[2][2];

      dim2outprod(s, y, HL);
      dim2outprod(y, s, HR);
      dim2outprod(s, s, LF);

      /*** HL = I - rho * HL; IR = I - rho * HR ***/
      HL[0][0] = 1 - rhob * HL[0][0];
      HL[0][1] = -rhob * HL[0][1];
      HL[1][0] = -rhob * HL[1][0];
      HL[1][1] = 1 - rhob * HL[1][1];

      HR[0][0] = 1 - rhob * HR[0][0];
      HR[0][1] = -rhob * HR[0][1];
      HR[1][0] = -rhob * HR[1][0];
      HR[1][1] = 1 - rhob * HR[1][1];

      /*** LF *= rho **/
      LF[0][0] *= rhob;
      LF[0][1] *= rhob;
      LF[1][0] *= rhob;
      LF[1][1] *= rhob;

      Mmatrix<double>::mtx22Multmtx22(HL, H, Tmp);
      Mmatrix<double>::mtx22Multmtx22(Tmp, HR, H);

      H[0][0] += LF[0][0];
      H[0][1] += LF[0][1];
      H[1][0] += LF[1][0];
      H[1][1] += LF[1][1];
    }
    count += 1;
  }
}

void lmpDriver::idealLmpRun(const char tag[]) {
  char lmpcmds[100][MAXSTRL];
  int i = 0;
  if (myRank == LmpRoot) {
    //  --------------------- INITIALIZAITION ---------------------
    sprintf(lmpcmds[i++], "clear");
    sprintf(lmpcmds[i++], "units  metal");
    sprintf(lmpcmds[i++], "dimension  3");
    sprintf(lmpcmds[i++], "boundary p p p");
    // sprintf(lmpcmds[i++],"variable  a equal  %.7f", lattice);
    // --------------------- ATOM DEFINITION ---------------------
    /*** strain is added by change the basis vector ***/
    if (((strcasecmp(tag, "tpath")) == 0) &&
        ((strcasecmp(structname, "bcc")) == 0)) {
      sprintf(lmpcmds[i++],
              "lattice custom %f "
              "a1 %8.6f %8.6f %8.6f a2 %8.6f %8.6f %8.6f a3 %8.6f %8.6f "
              "%8.6f "
              "basis 0 0 0 basis 0.5 0.5 0.5",
              lattice, baseVect[0][0], baseVect[0][1], baseVect[0][2],
              baseVect[1][0], baseVect[1][1], baseVect[1][2], baseVect[2][0],
              baseVect[2][1], baseVect[2][2]);
    }

    else if (((strcasecmp(tag, "opath")) == 0) &&
             ((strcasecmp(structname, "bcc")) == 0)) {
      sprintf(lmpcmds[i++],
              "lattice custom %f "
              "a1 %8.6f %8.6f %8.6f a2 %8.6f %8.6f %8.6f a3 %8.6f %8.6f "
              "%8.6f "
              "basis 0.0 0.0 0.0 basis 0.5 0.5 0.0 "
              "basis 0.5 0.0 0.5 basis 0.0 0.5 0.5 ",
              lattice, baseVect[0][0], baseVect[0][1], baseVect[0][2],
              baseVect[1][0], baseVect[1][1], baseVect[1][2], baseVect[2][0],
              baseVect[2][1], baseVect[2][2]);
    }

    /*** use units lattice ! ***/
    sprintf(lmpcmds[i++],
            "region    box   block  0  1  0   1   0   1  units lattice");

    sprintf(lmpcmds[i++], "create_box    1   box");
    sprintf(lmpcmds[i++], "create_atoms  1   region   box");

    // --------------------- FORCE FIELDS ---------------------
    sprintf(lmpcmds[i++], "pair_style  %s", pairStyle);
    sprintf(lmpcmds[i++], "pair_coeff  * * %s %s", potName, eleName);
    sprintf(lmpcmds[i++], "neighbor   1.0  bin");
    sprintf(lmpcmds[i++], "neigh_modify  every 1  delay  0 check yes");

    sprintf(lmpcmds[i++], "compute  mystress  all stress/atom  NULL");

    sprintf(lmpcmds[i++], "compute p1 all reduce sum c_mystress[1]");
    sprintf(lmpcmds[i++], "compute p2 all reduce sum c_mystress[2]");
    sprintf(lmpcmds[i++], "compute p3 all reduce sum c_mystress[3]");
    sprintf(lmpcmds[i++], "compute p4 all reduce sum c_mystress[4]");
    sprintf(lmpcmds[i++], "compute p5 all reduce sum c_mystress[5]");
    sprintf(lmpcmds[i++], "compute p6 all reduce sum c_mystress[6]");

    /*** stress unts is changed from bars to Mpa x 0.1 to Gpa x 1e-4 ***/
    sprintf(lmpcmds[i++], "variable  Sxx equal 1e-4*(c_p1)/vol");  // xx
    sprintf(lmpcmds[i++], "variable  Syy equal 1e-4*(c_p2)/vol");  // yy
    sprintf(lmpcmds[i++], "variable  Szz equal 1e-4*(c_p3)/vol");  // zz
    sprintf(lmpcmds[i++], "variable  Sxy equal 1e-4*(c_p4)/vol");  // zz
    sprintf(lmpcmds[i++], "variable  Sxz equal 1e-4*(c_p5)/vol");  // zz
    sprintf(lmpcmds[i++], "variable  Syz equal 1e-4*(c_p6)/vol");  // zz

    // ---------------------- THERMO  --------------------------
    sprintf(lmpcmds[i++], "thermo  10000");
    sprintf(lmpcmds[i++],
            "thermo_style custom  step  pe  etotal  lx   ly   lz   v_Sxx  "
            "v_Syy  v_Szz  v_Sxy  v_Sxz  v_Syz");

    // sprintf(lmpcmds[i++],"dump  1   all  xyz  5000  xyz/bcc.xyz");
    // sprintf(lmpcmds[i++], "dump    1   all     cfg     5000    cfg/W*.cfg
    // mass type xs ys zs");

    // ---------------------- RELAX   --------------------------
    sprintf(lmpcmds[i++], "min_style  cg");
    sprintf(lmpcmds[i++], "minimize  1e-10  1e-10  100000  100000");
  }

  for (int iter = 0; iter < i; iter++) lammps_command(lmp, lmpcmds[iter]);

  /*** extract energy ***/
  lmpIdealTensileEngy = lammps_get_thermo(lmp, ETOTAL);

  /*******************************************************************
   * extract stress
   * *****************************************************************/
  stressVec[0] = *(double*)lammps_extract_variable(lmp, "Sxx", NULL);
  stressVec[1] = *(double*)lammps_extract_variable(lmp, "Syy", NULL);
  stressVec[2] = *(double*)lammps_extract_variable(lmp, "Szz", NULL);
  stressVec[3] = *(double*)lammps_extract_variable(lmp, "Sxy", NULL);
  stressVec[4] = *(double*)lammps_extract_variable(lmp, "Sxz", NULL);
  stressVec[5] = *(double*)lammps_extract_variable(lmp, "Syz", NULL);
}

/***************************************************************
 * update strain based on calculated stress
 * strainVect = -SijMtx6x6 * StressVect * coeff
 * say 20% strain stress 5e4 [Mpa] -> 50 Gpa
 * 50 Gpa tune -> 0.005  in strain
 * so coeff = 0.005 / 50 = 0.0001
 * **************************************************************/

void lmpDriver::calIdealtenLinear(double delta, char tag[], int id) {
  /*** preparation: declare and define local variables  ***/

  bool screentag = false;
  double tempMtx[DIM][DIM];
  double perfBase[DIM][DIM];
  double s11, s12, s44;
  double curval;
  double strainVect[STRESS_N] = {0., 0., 0., 0., 0., 0.};
  double coeff[DIM] = {2., 2., 2.};

  Melem myEle;
  myEle.initElastic();

  for (int i = 0; i < int(myEle.cijlist.size()); i++) {
    if ((strcasecmp("W", myEle.cijlist[i].name)) == 0) {
      s11 = myEle.cijlist[i].cubicSij[0];
      s12 = myEle.cijlist[i].cubicSij[1];
      s44 = myEle.cijlist[i].cubicSij[2];
    }
  }

  if (screentag) printf("c11 %f  c12 %f  c44 %f\n", s11, s12, s44);

  /** initialize the base vector  **/
  if (strcasecmp(tag, "tpath") == 0) {
    baseVect[0][0] = 1.0;
    baseVect[1][1] = 1.0;
    baseVect[2][2] = 1.0;

    /*** store the pefect Base_vector ***/
    Mmatrix<double>::mtx33copymtx33(perfBase, baseVect);

  } else if (strcasecmp(tag, "opath") == 0) {
    baseVect[0][0] = 1.0;
    baseVect[1][1] = SQRT2;
    baseVect[2][2] = SQRT2;

    /*** store the pefect Base_vector ***/
    Mmatrix<double>::mtx33copymtx33(perfBase, baseVect);
  }

  /*** loop of ideal tensile ***/
  int count = 0;
  FILE* ptFile;

  /** initialize strain  **/
  strainMtx[0][0] = 1 + delta;
  strainMtx[1][1] = 1.0 + PERTERB * delta;  // some prterbation
  strainMtx[2][2] = 1.0 - PERTERB * delta;  // some perterbation

  /** add strain on baseVect **/
  Mmatrix<double>::mtx33copymtx33(tempMtx, baseVect);
  Mmatrix<double>::mtx33Multmtx33(strainMtx, tempMtx, baseVect);

  /*** run lammps ***/
  idealLmpRun(tag);

  /*** currently curval value ***/
  curval = MAX(fabs(stressVec[1]), fabs(stressVec[2]));

  /*** print the info ***/
  if (screentag) {
    if (myRank == LmpRoot) {
      if ((ptFile = fopen("relax.log", "w")) == NULL) {
        fprintf(stderr, "error when open relax.log. %s:%d", __FILE__, __LINE__);
        exit(1);
      };

      fprintf(ptFile, "target %f\n", curval);
      fprintf(ptFile, "[%9.8f %9.8f %9.8f] ", strainMtx[0][0], strainMtx[1][1],
              strainMtx[2][2]);
      fprintf(ptFile, "[%9.8f %9.8f %9.8f]\n", stressVec[0], stressVec[1],
              stressVec[2]);

      fflush(ptFile);
      fclose(ptFile);
    }
  }

  while (true) {
    if (curval < STRESS_SMALL) {
      /** Syy and Szz are both smaller than 0.05 ***/
      Mmatrix<double>::mtx33copymtx33(gstrainMtx, strainMtx);

      /*** check whether this is equilibrium calculation ***/
      if (fabs(baseVect[0][0] - 1.0) > EPSILON0) {
        /** stored the original energy difference for observation **/
        if (strcasecmp(tag, "tpath") == 0) {
          /*** calculate the relative energy to equilibrium ***/
          double diff1 = lmpIdealTensileEngy - equiDFTData.tpathEnergy;

          engyListLMP[id] = engyDiffLmp = diff1;

          if (mycalOpt.itentag == 1) {
            /*** use quadratic function ***/

            /*** calculate the relative difference for idealEngyErr
             * ***/
            diff1 -= dftDataList[id].tpathEngyIncre;
            idealEngyErr = diff1 * diff1;

            /*** calculate the stresErr ***/
            double temp;
            idealStssErr = idealStrnErr = 0;

            for (int k = 0; k < 3; k++) {
              temp = stressVec[k] - dftDataList[id].tpathStressV[k];
              idealStssErr += (temp * temp);
            }
            for (int k = 1; k < 3; k++) {
              temp = strainMtx[k][k] - dftDataList[id].tpathBaseMtx[k][k];
              idealStrnErr += (temp * temp);
            }

          } else if (mycalOpt.itentag == 2) {
            /*** use self tuned function ***/
            idealEngyErr = targetFunct(engyDiffLmp, tp_engy[id]);

            /*** currently only compare the sigma_xx ***/
            idealStssErr = targetFunct(stressVec[0], tp_stress[id]);
          }

          /*** store energy Err ***/
          itenEngyErrList[id] = weighTEngy[id] * idealEngyErr;

          /*** store stress Err ***/
          itenStssErrList[id] = weighTStss[id] * idealStssErr;

          /*** store strain Err ***/
          itenStrErrList[id] =
              0.2 * STRAIN_WEIGHT_COE * weighTStss[id] * idealStrnErr;

        } else if (strcasecmp(tag, "opath") == 0) {
          double diff1 = lmpIdealTensileEngy - equiDFTData.opathEnergy;

          /** stored the original energy difference for observation
           * **/
          engyListLMP[id + TPOINTS] = engyDiffLmp = diff1;

          if (mycalOpt.itentag == 1) {
            /*** use quadratic function ***/

            /*** calculate the relative difference for idealEngyErr
             * ***/
            diff1 -= dftDataList[id].opathEngyIncre;
            idealEngyErr = diff1 * diff1;

            /*** calculate the stresErr ***/
            double temp;
            idealStrnErr = idealStssErr = 0.0;

            for (int k = 0; k < 3; k++) {
              temp = stressVec[k] - dftDataList[id].opathStressV[k];
              idealStssErr += (temp * temp);
            }

            for (int k = 1; k < 3; k++) {
              temp =
                  strainMtx[k][k] - dftDataList[id].opathBaseMtx[k][k] / SQRT2;
              idealStrnErr += (temp * temp);
            }

          } else if (mycalOpt.itentag == 2) {
            /*** use self tuned function ***/
            idealEngyErr = targetFunct(engyDiffLmp, tp_engy[id]);

            /*** currently only compare the sigma_xx ***/
            idealStssErr = targetFunct(stressVec[0], tp_stress[id]);
          }
          /*** store energy Err ***/
          itenEngyErrList[id + TPOINTS] = weighOEngy[id] * idealEngyErr;

          /*** store stress Err ***/
          itenStssErrList[id + TPOINTS] = weighOStss[id] * idealStssErr;

          /*** store strain Err ***/
          itenStrErrList[id + TPOINTS] =
              STRAIN_WEIGHT_COE * weighOStss[id] * idealStrnErr;
        }

      } else {
        /** equilibrium case we store energy **/
        if (strcasecmp(tag, "tpath") == 0) {
          equiDFTData.tpathEnergy = lmpIdealTensileEngy;
        } else if (strcasecmp(tag, "opath") == 0) {
          equiDFTData.opathEnergy = lmpIdealTensileEngy;
        }
      }
      /* exit */
      break;

    } else if (count > MAXBFGS_STEP) {
      printf(" run too many times, tune your parameters \n");
      break;
    } else {
      /*** update the strain ***/
      // strainVect[0] = coeff1 * (s11 * 0 + s12 * stressVec[1] + s12 *
      // stressVec[2]); strainVect[1] = coeff1 * (s12 * 0 + s11 *
      // stressVec[1] + s12 * stressVec[2]); strainVect[2] = coeff2 * (s12
      // * 0 + s12 * stressVec[1] + s11 * stressVec[2]);

      // strainVect[3] = coeff1 * (s44 * stressVec[3]);
      // strainVect[4] = coeff1 * (s44 * stressVec[4]);
      // strainVect[5] = coeff1 * (s44 * stressVec[5]);

      /*** non coupleing ***/
      for (int it = 1; it <= 2; it++) {
        if (fabs(stressVec[it]) > STRESS_SMALL) {
          /*** modify strain ***/
          strainVect[it] = coeff[it] * s12 * stressVec[it];
        } else {
          /*** already satisfies (smaller step )***/
          strainVect[it] = 0.00001 * s12 * stressVec[it];
        }
      }

      strainMtx[1][1] += strainVect[1];
      strainMtx[2][2] += strainVect[2];

      /** add the strain on baseVect **/
      Mmatrix<double>::mtx33Multmtx33(strainMtx, perfBase, baseVect);

      /*** ideal tensile we know only yy, zz maters ***/
      if (screentag) {
        printf(
            "stressV    %.6f %.6f\n"
            "strainVect %.6f %.6f\n"
            "baseVect[1][1], baseVect[2][2]  %.6f  %.6f\n",
            stressVec[1], stressVec[2], strainVect[1], strainVect[2],
            baseVect[1][1], baseVect[2][2]);
      }

      /** run **/
      idealLmpRun(tag);
      count += 1;

      /*** print the info ***/
      if (screentag) {
        if (myRank == LmpRoot) {
          if ((ptFile = fopen("relax.log", "a")) == NULL) {
            fprintf(stderr, "error when open relax.log. %s:%d", __FILE__,
                    __LINE__);
            exit(1);
          };

          fprintf(ptFile, "%d -> %f\n", count, curval);
          fprintf(ptFile, "[%9.8f %9.8f %9.8f]  ", strainMtx[0][0],
                  strainMtx[1][1], strainMtx[2][2]);
          fprintf(ptFile, "[%9.8f %9.8f %9.8f]\n", stressVec[0], stressVec[1],
                  stressVec[2]);

          fflush(ptFile);
          fclose(ptFile);
        }
      }
      curval = temp = MAX(fabs(stressVec[1]),
                          fabs(stressVec[2]));  // function value at xk+1;

    } /*** curval < EPSILON ***/
  }   /*** while true ***/

  return;
}
