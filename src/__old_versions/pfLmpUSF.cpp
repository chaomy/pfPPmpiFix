/*
* @Author: chaomy
* @Date:   2017-11-10 14:40:55
* @Last Modified by:   chaomy
* @Last Modified time: 2017-11-22 00:31:09
*/

#include "lmpDriver.h"

/*******************************************************************
* calculate the generalized stacking fault energy
* gsf curve  [111](110); [111](211); [111](123);
* *****************************************************************/

void lmpDriver::calGSFenergy(){
    char lmpcmds[100][MAXSTRL];
    double disp  = 0;
    const int    samplePoints = GSF_POINTS;
    const double delta = 1./(samplePoints - 1);
    double gsf111_sur110[samplePoints];
    double gsf111_sur211[samplePoints];
    double equiEnergy = 0;
    double lx, ly, area;
    double dx = 0.0;
    int i;
    errGsf = 0.0;
    /*******************************************************************
    * gsf curve  [111](110);   20 points
    * x [-1  1  1]
    * y [1, -1  2]
    * z [1,  1, 0]
    * *****************************************************************/
    for (int iter = 0; iter < samplePoints; iter ++){
        i = 0;
        disp = double(iter * delta);

        if (myRank == 0){
            //  --------------------- INITIALIZAITION ---------------------
            sprintf(lmpcmds[i++],"clear");
            sprintf(lmpcmds[i++],"units  metal");
            sprintf(lmpcmds[i++],"dimension  3");
            sprintf(lmpcmds[i++],"boundary p p p");
            sprintf(lmpcmds[i++],"atom_style atomic");
            sprintf(lmpcmds[i++],"variable  a equal  %.7f", lattice);
            sprintf(lmpcmds[i++],"variable  widthx  equal ${a}*sqrt(3)/2");
            sprintf(lmpcmds[i++],"variable  widthy  equal ${a}*sqrt(6)");
            sprintf(lmpcmds[i++],"variable  zlow   equal  ${a}*sqrt(2)*10");
            sprintf(lmpcmds[i++],"variable  zhigh  equal  ${a}*sqrt(2)*20");
            sprintf(lmpcmds[i++],"variable  ztotal  equal ${a}*sqrt(2)*23");
            sprintf(lmpcmds[i++],"variable  x_displace equal  %f*${widthx}", disp);

            // --------------------- ATOM DEFINITION ---------------------
            sprintf(lmpcmds[i++],"lattice bcc  ${a}");
            sprintf(lmpcmds[i++],"region	whole  block  0  ${widthx}  0  ${widthy}  0  ${ztotal}  units box");
            sprintf(lmpcmds[i++],"create_box  1   whole");
            sprintf(lmpcmds[i++],"lattice   bcc  ${a} orient x -1 1 1 orient y 1 -1  2 orient  z 1 1 0");
            sprintf(lmpcmds[i++],"region	mybody  block  0  ${widthx}  0  ${widthy}  0  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"create_atoms   1  region   mybody");
            sprintf(lmpcmds[i++],"region 	 top   block   0  ${widthx}  0  ${widthy}  ${zlow}  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"group    gtop  region   top");

            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++],"pair_style  %s", pairStyle);
            sprintf(lmpcmds[i++],"pair_coeff  *  *  %s %s", potName, eleName);
            sprintf(lmpcmds[i++],"neighbor 2.0 bin");
            sprintf(lmpcmds[i++],"neigh_modify  every 1  delay  5  check yes");

            // ---------------------- THERMO  --------------------------
            sprintf(lmpcmds[i++],"thermo 5000");
            sprintf(lmpcmds[i++],"thermo_style  custom  lx  ly  lz step  etotal");
            // sprintf(lmpcmds[i++], "dump  1  all   cfg   500000  cfg/Nb_*_%d.cfg  mass  type xs ys zs", iter);
            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++], "displace_atoms  gtop    move   ${x_displace}  0  0  units box");

            // ---------------------- RELAX   --------------------------
            sprintf(lmpcmds[i++], "fix  1  all  setforce  0  0  0");
            sprintf(lmpcmds[i++],"min_style  cg");
            sprintf(lmpcmds[i++],"minimize  1e-12  1e-12  100000  100000");
            sprintf(lmpcmds[i++], "unfix  1");
        }

        MPI_Bcast(&i, 1, MPI_INT, LmpRoot, lmp_comm);

        for (int iter = 0; iter < i; iter ++){
            MPI_Bcast(lmpcmds[iter], MAXSTRL, MPI_CHAR, LmpRoot, lmp_comm);
            lammps_command(lmp, lmpcmds[iter]);
        }

        if (myRank == LmpRoot){
            /*******************************************************************
            * extract the lattice constant
            * *****************************************************************/
            gsf111_sur110[iter] = lammps_get_thermo(lmp, "etotal");
            if (iter == 0){
                lx = lammps_get_thermo(lmp, "lx");
                ly = lammps_get_thermo(lmp, "ly");
                area = lx * ly;
                equiEnergy = gsf111_sur110[0];
            }
            /* 
            printf("###################################\n"
                    "the energy total  is  %f\n"
                    "###################################\n", gsf111_sur110[iter]);
            */
            gsf111_sur110[iter] = (gsf111_sur110[iter] - equiEnergy) / area;
            dx = (gsf111_sur110[iter] - dft_gsf_sur110[iter]);
            errGsf += (gsf_weight * dx * dx);
        }
        lammps_command(lmp, "clear");

    MPI_Barrier(lmp_comm);
    }

    /*******************************************************************
    * gsf curve  [111](211);   20 points
    * x [-1  1  1]
    * y [ 0 -1  1
    * z [ 2  1  1]
    * *****************************************************************/
    for (int iter = 0; iter < samplePoints; iter ++){
        i = 0;
        disp = double(iter * delta);

        if (myRank == LmpRoot){
            //  --------------------- INITIALIZAITION ---------------------
            sprintf(lmpcmds[i++],"clear");
            sprintf(lmpcmds[i++],"units  metal");
            sprintf(lmpcmds[i++],"dimension  3");
            sprintf(lmpcmds[i++],"boundary p p p");
            sprintf(lmpcmds[i++],"atom_style atomic");
            sprintf(lmpcmds[i++],"variable  a equal  %.7f", lattice);
            sprintf(lmpcmds[i++],"variable  widthx  equal  ${a}*sqrt(3)/2");
            sprintf(lmpcmds[i++],"variable  widthy  equal  ${a}*sqrt(2)");
            sprintf(lmpcmds[i++],"variable  zlow    equal  ${a}*sqrt(6)*7");
            sprintf(lmpcmds[i++],"variable  zhigh   equal  ${a}*sqrt(6)*14");
            sprintf(lmpcmds[i++],"variable  ztotal  equal  ${a}*sqrt(6)*18");
            sprintf(lmpcmds[i++],"variable  x_displace equal %f*${widthx}", disp);

            // --------------------- ATOM DEFINITION ---------------------
            sprintf(lmpcmds[i++],"lattice bcc  ${a}");
            sprintf(lmpcmds[i++],"region	whole  block  0  ${widthx}  0  ${widthy}  0  ${ztotal}  units box");
            sprintf(lmpcmds[i++],"create_box  1   whole");
            sprintf(lmpcmds[i++],"lattice   bcc  ${a} orient x -1 1 1  orient y 0 -1 1  orient  z 2 1 1");
            sprintf(lmpcmds[i++],"region	mybody  block  0  ${widthx}  0  ${widthy}  0  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"create_atoms   1  region   mybody");
            sprintf(lmpcmds[i++],"region 	 top   block   0  ${widthx}  0  ${widthy}  ${zlow}  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"group    gtop  region   top");

            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++],"pair_style  %s", pairStyle);
            sprintf(lmpcmds[i++],"pair_coeff  *  *  %s %s", potName, eleName);
            sprintf(lmpcmds[i++],"neighbor 2.0 bin");
            sprintf(lmpcmds[i++],"neigh_modify  every 1  delay  0  check yes");

            // ---------------------- THERMO  --------------------------
            sprintf(lmpcmds[i++],"thermo 5000");
            sprintf(lmpcmds[i++],"thermo_style  custom  lx  ly  lz step  etotal");
            // sprintf(lmpcmds[i++], "dump  1  all   cfg   500000  cfg/Nb_*_%d.cfg  mass  type xs ys zs", iter);
            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++], "displace_atoms  gtop    move   ${x_displace}  0  0  units box");

            // ---------------------- RELAX   --------------------------
            sprintf(lmpcmds[i++], "fix  1  all  setforce  0  0  0");
            sprintf(lmpcmds[i++],"min_style  cg");
            sprintf(lmpcmds[i++],"minimize  1e-12  1e-12  100000  100000");
            sprintf(lmpcmds[i++], "unfix  1");
        }

        MPI_Bcast(&i, 1, MPI_INT, LmpRoot, lmp_comm);

        for (int iter = 0; iter < i; iter ++){
            MPI_Bcast(lmpcmds[iter], MAXSTRL, MPI_CHAR, LmpRoot, lmp_comm);
            lammps_command(lmp, lmpcmds[iter]);
        }

        if (myRank == LmpRoot){
            /*******************************************************************
            * extract the gsf energy and area
            * *****************************************************************/
            gsf111_sur211[iter] = lammps_get_thermo(lmp, "etotal");

            if (iter == 0){
                lx = lammps_get_thermo(lmp, "lx");
                ly = lammps_get_thermo(lmp, "ly");
                /*******************************************************************
                 * get area and equilibrium energy
                * *****************************************************************/
                area = lx * ly;
                equiEnergy = gsf111_sur211[0];
            }

            gsf111_sur211[iter] = (gsf111_sur211[iter] - equiEnergy) / area;

            dx = (gsf111_sur211[iter] - dft_gsf_sur211[iter]);
            errGsf += (gsf_weight * dx * dx);
        }
        lammps_command(lmp, "clear");

    MPI_Barrier(lmp_comm);
    }

    if (myRank == LmpRoot){
        for (int j = 0; j < samplePoints; j ++){
            printf("gsf[%.4f]  gsf[%.4f]\n",
                    gsf111_sur110[j], gsf111_sur211[j]);
        }
    }

    errPhy += errGsf;

    MPI_Barrier(lmp_comm);
    return;
}

void lmpDriver::calUstackingFault(){
    /*******************************************************************
    * calculate the unstable stacking fault
    * gsf curve  [111](110); [111](211); [111](123);
    * *****************************************************************/
    char lmpcmds[100][MAXSTRL];
    double disp  = 0;

    /*** initialization of global variable ***/
    ustacking_fault_110 = 0.0;
    ustacking_fault_211 = 0.0;
    errUssf = 0.0;

    double lx, ly, area;
    double dx = 0.0;
    double delta = 0.5;
    double equiEnergy = 0.0;
    int    i;

    /*******************************************************************
    * ustf  [111](110); 2 point
    * x [-1  1  1]
    * y [1, -1  2]
    * z [1,  1, 0]
    * *****************************************************************/
    for (int iter = 0; iter < 2; iter++){
        i = 0;
        disp = double(iter * delta);

        if (myRank == 0){
            //  --------------------- INITIALIZAITION ---------------------
            sprintf(lmpcmds[i++],"clear");
            sprintf(lmpcmds[i++],"units  metal");
            sprintf(lmpcmds[i++],"dimension  3");
            sprintf(lmpcmds[i++],"boundary p p p");
            sprintf(lmpcmds[i++],"atom_style atomic");
            sprintf(lmpcmds[i++],"variable  a equal  %.7f", lattice);
            sprintf(lmpcmds[i++],"variable  widthx  equal ${a}*sqrt(3)/2");
            sprintf(lmpcmds[i++],"variable  widthy  equal ${a}*sqrt(6)");
            sprintf(lmpcmds[i++],"variable  zlow   equal  ${a}*sqrt(2)*10");
            sprintf(lmpcmds[i++],"variable  zhigh  equal  ${a}*sqrt(2)*20");
            sprintf(lmpcmds[i++],"variable  ztotal  equal ${a}*sqrt(2)*23");
            sprintf(lmpcmds[i++],"variable  x_displace equal  %f*${widthx}", disp);

            // --------------------- ATOM DEFINITION ---------------------
            sprintf(lmpcmds[i++],"lattice bcc  ${a}");
            sprintf(lmpcmds[i++],"region	whole  block  0  ${widthx}  0  ${widthy}  0  ${ztotal}  units box");
            sprintf(lmpcmds[i++],"create_box  1   whole");
            sprintf(lmpcmds[i++],"lattice   bcc  ${a} orient x -1 1 1 orient y 1 -1  2 orient  z 1 1 0");
            sprintf(lmpcmds[i++],"region	mybody  block  0  ${widthx}  0  ${widthy}  0  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"create_atoms   1  region   mybody");
            sprintf(lmpcmds[i++],"region 	 top   block   0  ${widthx}  0  ${widthy}  ${zlow}  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"group    gtop  region   top");

            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++],"pair_style  %s", pairStyle);
            sprintf(lmpcmds[i++],"pair_coeff  *  *  %s %s", potName, eleName);
            sprintf(lmpcmds[i++],"neighbor 2.0 bin");
            sprintf(lmpcmds[i++],"neigh_modify  every 1  delay  0  check yes");

            // ---------------------- THERMO  --------------------------
            sprintf(lmpcmds[i++],"thermo 5000");
            sprintf(lmpcmds[i++],"thermo_style  custom  lx  ly  lz step  etotal");
            // sprintf(lmpcmds[i++], "dump  1  all   cfg   500000  cfg/Nb_*_%d.cfg  mass  type xs ys zs", iter);
            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++], "displace_atoms  gtop    move   ${x_displace}  0  0  units box");

            // ---------------------- RELAX   --------------------------
            sprintf(lmpcmds[i++],"fix    1   all  setforce  0  0  0");
            sprintf(lmpcmds[i++],"min_style  cg");
            sprintf(lmpcmds[i++],"minimize   1e-12  1e-12  100000  100000");
            sprintf(lmpcmds[i++],"unfix      1");
        }

        MPI_Bcast(&i, 1, MPI_INT, LmpRoot, lmp_comm);

        for (int iter = 0; iter < i; iter ++){
            MPI_Bcast(lmpcmds[iter], MAXSTRL, MPI_CHAR, LmpRoot, lmp_comm);
            lammps_command(lmp, lmpcmds[iter]);
        }
        if (myRank == LmpRoot){
            /*******************************************************************
             * extract gsf energy
            * *****************************************************************/
            ustacking_fault_110 = lammps_get_thermo(lmp, "etotal");

            if (iter == 0){
                /*******************************************************************
                 * get area and equilibrium energy
                * *****************************************************************/

                lx = lammps_get_thermo(lmp, "lx");
                ly = lammps_get_thermo(lmp, "ly");
                area = lx * ly;
                equiEnergy = ustacking_fault_110;
            }else{
                ustacking_fault_110 = (ustacking_fault_110 - equiEnergy) / area;

                if (mycalOpt.usstag == 1){
                    /*** use quadratic function ***/
                    dx = (ustacking_fault_110 - dft_ustacking_fault_110);
                    ussfErr[0] = (ussf_weight * dx * dx);
                    errUssf += ussfErr[0];

                }else if (mycalOpt.usstag == 2){
                    /*** use tuned function ***/
                    ussfErr[0] = (ussf_weight * targetFunct(ustacking_fault_110, ussf110));
                    errUssf += ussfErr[0];
                }
            }
        }
        lammps_command(lmp, "clear");

    MPI_Barrier(lmp_comm);
    }

    /*******************************************************************
    * ussf  [111](211);   2 points
    * x [-1  1  1]
    * y [ 0 -1  1
    * z [ 2  1  1]
    * *****************************************************************/
    for (int iter = 0; iter < 2; iter ++){
        i = 0;
        disp = double(iter * delta);

        if (myRank == LmpRoot){
            //  --------------------- INITIALIZAITION ---------------------
            sprintf(lmpcmds[i++],"clear");
            sprintf(lmpcmds[i++],"units  metal");
            sprintf(lmpcmds[i++],"dimension  3");
            sprintf(lmpcmds[i++],"boundary p p p");
            sprintf(lmpcmds[i++],"atom_style atomic");
            sprintf(lmpcmds[i++],"variable  a equal  %.7f", lattice);
            sprintf(lmpcmds[i++],"variable  widthx  equal  ${a}*sqrt(3)/2");
            sprintf(lmpcmds[i++],"variable  widthy  equal  ${a}*sqrt(2)");
            sprintf(lmpcmds[i++],"variable  zlow    equal  ${a}*sqrt(6)*7");
            sprintf(lmpcmds[i++],"variable  zhigh   equal  ${a}*sqrt(6)*14");
            sprintf(lmpcmds[i++],"variable  ztotal  equal  ${a}*sqrt(6)*18");
            sprintf(lmpcmds[i++],"variable  x_displace equal %f*${widthx}", disp);

            // --------------------- ATOM DEFINITION ---------------------
            sprintf(lmpcmds[i++],"lattice bcc  ${a}");
            sprintf(lmpcmds[i++],"region	whole  block  0  ${widthx}  0  ${widthy}  0  ${ztotal}  units box");
            sprintf(lmpcmds[i++],"create_box  1   whole");
            sprintf(lmpcmds[i++],"lattice   bcc  ${a} orient x -1 1 1 orient y 0 -1 1  orient  z 2 1 1");
            sprintf(lmpcmds[i++],"region	mybody  block  0  ${widthx}  0  ${widthy}  0  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"create_atoms   1  region   mybody");
            sprintf(lmpcmds[i++],"region 	 top   block   0  ${widthx}  0  ${widthy}  ${zlow}  ${zhigh}  units box");
            sprintf(lmpcmds[i++],"group    gtop  region   top");

            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++],"pair_style  %s", pairStyle);
            sprintf(lmpcmds[i++],"pair_coeff  *  *  %s %s", potName, eleName);
            sprintf(lmpcmds[i++],"neighbor 2.0 bin");
            sprintf(lmpcmds[i++],"neigh_modify  every 1  delay  0  check yes");

            // ---------------------- THERMO  --------------------------
            sprintf(lmpcmds[i++],"thermo 5000");
            sprintf(lmpcmds[i++],"thermo_style  custom  lx  ly  lz step  etotal");
            // sprintf(lmpcmds[i++], "dump  1  all   cfg   500000  cfg/Nb_*_%d.cfg  mass  type xs ys zs", iter);
            // --------------------- FORCE FIELDS ---------------------
            sprintf(lmpcmds[i++], "displace_atoms  gtop    move   ${x_displace}  0  0  units box");

            // ---------------------- RELAX   --------------------------
            sprintf(lmpcmds[i++], "fix  1  all  setforce  0  0  0");
            sprintf(lmpcmds[i++],"min_style  cg");
            sprintf(lmpcmds[i++],"minimize  1e-12   1e-12  100000  100000");
            sprintf(lmpcmds[i++], "unfix  1");
        }
        MPI_Bcast(&i, 1, MPI_INT, LmpRoot, lmp_comm);

        for (int iter = 0; iter < i; iter ++){
            MPI_Bcast(lmpcmds[iter], MAXSTRL, MPI_CHAR, LmpRoot, lmp_comm);
            lammps_command(lmp, lmpcmds[iter]);
        }

        if (myRank == LmpRoot){
            /*******************************************************************
            * extract the gsf energy and area
            * *****************************************************************/
            ustacking_fault_211 = lammps_get_thermo(lmp, "etotal");

            if (iter == 0){
                lx = lammps_get_thermo(lmp, "lx");
                ly = lammps_get_thermo(lmp, "ly");
                /*******************************************************************
                 * get area and equilibrium energy
                * *****************************************************************/
                area = lx * ly;
                equiEnergy = ustacking_fault_211;
            }else{
                ustacking_fault_211 = (ustacking_fault_211 - equiEnergy) / area;
                if (mycalOpt.usstag == 1){
                    dx = (ustacking_fault_211 - dft_ustacking_fault_211);
                    ussfErr[1] = (ussf_weight * dx * dx);
                    errUssf += ussfErr[1];

                }else if (mycalOpt.usstag == 2){
                    // printf("ulmp , udft = %f, %f \n",ustacking_fault_211, ussf211.targ_x);
                    ussfErr[1] = (ussf_weight * targetFunct(ustacking_fault_211, ussf211));
                    errUssf += ussfErr[1];
                }
            }
        }
        lammps_command(lmp, "clear");

    errPhy += errUssf;
    MPI_Barrier(lmp_comm);
    }
}
