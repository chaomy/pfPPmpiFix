/*
* @Author: yangchaoming
* @Date:   2017-10-23 15:52:29
* @Last Modified by:   chaomy
* @Last Modified time: 2017-11-14 16:41:30
*/

#include "pfHome.h"

double pfHome::forceEAM(vector<Func>& ffs, int tag){
	while(true){
		double gerr = 0.0;	
		for (int i = 0; i <nfuncs; i++){
			MPI_Bcast(&ffs[i].yy.front(), ffs[i].yy.size(), MPI_DOUBLE, PFROOT, MPI_COMM_WORLD);
			MPI_Bcast(&ffs[i].g1.front(), ffs[i].g1.size(), MPI_DOUBLE, PFROOT, MPI_COMM_WORLD);
		}

		MPI_Bcast(&tag, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (tag == 1) break; 
		
		// update splines
		for (int i = 0; i < nfuncs; i++) spline(ffs[i]);	
		double error = 0.0;
		double maxrho = -1e3;
		double minrho = 1e3; 
		double gminrho = 0.0; 
		double gmaxrho = 0.0; 

		for (int cc = locstt; cc < locend; cc++){
			Config& cnf = configs[cc]; 
			double fitengy = 0.0; 

			/*-------------- first loop over atoms to reset values  --------------*/
			for (int ii = 0; ii < cnf.natoms; ii++){
				pfAtom& atm = cnf.atoms[ii];	
				atm.rho = 0.0; 
				atm.fitfrc[0] = -atm.frc[0]; 
				atm.fitfrc[1] = -atm.frc[1]; 
				atm.fitfrc[2] = -atm.frc[2]; 
			} // ii

			/*-------------- second loop over atoms pairs, and densities  --------------*/
			for (int ii = 0; ii < cnf.natoms; ii++){
				pfAtom& atm = cnf.atoms[ii];	

				for (int nn = 0; nn < atm.nneighs; nn++){
					Neigh& ngb = atm.neighs[nn];   

					if (ngb.r < ffs[PHI].xx.back()){
						double phi = 0, phigrad = 0;	
						splint(ffs[PHI], ngb.slots[0], ngb.shifts[0], ngb.steps[0], phi, phigrad);

						fitengy += phi;
						double tmp[3]; 
						tmp[0] = ngb.dist2r[0] * phigrad;
						tmp[1] = ngb.dist2r[1] * phigrad;
						tmp[2] = ngb.dist2r[2] * phigrad;

						atm.fitfrc[0] += tmp[0];	
						atm.fitfrc[1] += tmp[1];
						atm.fitfrc[2] += tmp[2];

						cnf.atoms[ngb.aid].fitfrc[0] -= tmp[0];	
						cnf.atoms[ngb.aid].fitfrc[0] -= tmp[1];	
						cnf.atoms[ngb.aid].fitfrc[0] -= tmp[2];	
					}

					if (ngb.r < ffs[RHO].xx.back()){
						double rho; 	
						splint(ffs[RHO], ngb.slots[1], ngb.shifts[1], ngb.steps[1], rho, ngb.rhog);
						atm.rho += rho;
						cnf.atoms[ngb.aid].rho += rho; 
					}	
				} // nn

				double embE = 0.0; 
				double extra = 0.0; 
				if (atm.rho > ffs[EMF].xx.back()){
					extra = atm.rho - ffs[EMF].xx.back();	
					splint(ffs[EMF], ffs[EMF].xx.back(), embE, atm.gradF);
					embE += atm.gradF * extra;  	
				}

				if (atm.rho < ffs[EMF].xx.front()){
				    extra =	atm.rho - ffs[EMF].xx.front(); 
					splint(ffs[EMF], ffs[EMF].xx[0] - 0.5 * (ffs[EMF].xx[1] - ffs[EMF].xx[0]), 
						embE, atm.gradF);
					embE += atm.gradF * extra; 
				}else splint(ffs[EMF], atm.rho, embE, atm.gradF);

				maxrho = fmax(maxrho, atm.rho);
				minrho = fmin(minrho, atm.rho);	
				fitengy += embE;	
			} // ii

			/*-------------- third loop over atoms eambedding forces  --------------*/
			for (int ii = 0; ii < cnf.natoms; ii++){
				pfAtom& atm = cnf.atoms[ii];	

				for (int nn = 0; nn < atm.nneighs; nn++){
					Neigh& ngb = atm.neighs[nn];
					
					if (ngb.r < ffs[RHO].xx.back()){
						double rhoj = ngb.rhog;
						double emf = ngb.rhog * atm.gradF + rhoj * cnf.atoms[ngb.aid].gradF; 
						double tmp[3]; 

						tmp[0] = ngb.dist2r[0] * emf;
						tmp[1] = ngb.dist2r[1] * emf;
						tmp[2] = ngb.dist2r[2] * emf;

						atm.fitfrc[0] += tmp[0];
						atm.fitfrc[1] += tmp[1];
						atm.fitfrc[2] += tmp[2];
					}
				} // nn 
				scaleVec(atm.fitfrc, 1./(FRCEPS + atm.absfrc));
				error += square11(atm.fitfrc[0]) * atm.fweigh[0]; 		
				error += square11(atm.fitfrc[1]) * atm.fweigh[1];
				error += square11(atm.fitfrc[2]) * atm.fweigh[2];
			} // ii 	
			error += square11(fitengy - cnf.engy); 
		} // cc	

		MPI_Reduce(&error, &gerr, 1, MPI_DOUBLE, MPI_SUM, PFROOT, MPI_COMM_WORLD);	
		MPI_Reduce(&minrho, &gminrho, 1, MPI_DOUBLE, MPI_MIN, PFROOT, MPI_COMM_WORLD);
		MPI_Reduce(&maxrho, &gmaxrho, 1, MPI_DOUBLE, MPI_MAX, PFROOT, MPI_COMM_WORLD);
		if (mrank == PFROOT){
			gerr += LIMWEIGH*square11(ffs[EMF].xx.front() - gminrho);	
			gerr += LIMWEIGH*square11(ffs[EMF].xx.back() - gmaxrho); 
			return gerr;
		} 
	} // infinity loop 
	return -1;	
}