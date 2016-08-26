#include "FFT3d.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
int calc_dFdP_volAvg(ten4th dFdP_volAvg)
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
/* INFO: Inverts dP/dF_volAvg (stiffness) to get dF/dP_volAvg (compliance) 
	 * METHOD:
	 * 1) Convert 4th tensor to 9x9 matrix 
	 * 2) Reduce 9x9 to NxN matrix, where N is the number of "true" elements in P_BC_Flag
	 * 3) Invert reduced matrix
	 * 4) Reintroduce NxN matrix back to 9x9, other slots are filled with zeros. 
	 * 5) Convert 9x9 back to 4th tensor */ 
{
	int status;
	int i,j,k,m,n;
	double aux99[9][9] = {{0.0}};
	vec9 aux9;
	
	M33ToV9(P_BC_mask, aux9);
	M3333ToM99(dPdF_volAvg, aux99);
	
	// Reduce dPdF-9x9 based on P_BC
	k = -1;
	for(n=0;n<9;n++){
		if((int)aux9[n]){
			k++;
			j = -1;
			for(m=0;m<9;m++){
				if((int)aux9[m]){
					j++;
					dPdF_reduced[k][j] = aux99[n][m];
				}
			}
		}
	}
	
	// Invert dPdF --> dPdF
	status = matrixInverse(dPdF_reduced, dFdP_reduced, size_reduced);
	
	if(status){
		printf(BOLDRED "\nError: \"matrixInverse\" failed to invert dPdF_reduced" );
		printf("\ndPdF_reduced matrix is of dimensions (%dx%d) with following values:" RESET, size_reduced, size_reduced);
		for(i=0;i<size_reduced;i++){
			printf("\n");
			for(j=0;j<size_reduced;j++){
				printf("%.4f  ", dPdF_reduced[i][j]);
			}
		}
		printf("\n");
		return(status);
	}
	
	// Calc dFdP (zeros except where P_BC specified)
	for(i=0;i<9;i++){
		for(j=0;j<9;j++){
			aux99[i][j] = 0.0;
		}
	}
	
	k = -1; 
	for(n=0;n<9;n++){
		if(aux9[n]){
			k++;
			j = -1; 
			for(m=0;m<9;m++){
				if(aux9[m]){
					j++;
					aux99[n][m] = dFdP_reduced[k][j];
				}
			}
		}
	}
	
	// 9x9 --> 3x3x3x3 compliance
	M99ToM3333(aux99, dFdP_volAvg);
	
	return(0);
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void spectral_solver(void)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// Spectral solver for one load case only
{
	int i;
	int istep=0;
	int status=0;
	int jph;						// Phase index
	int status_writeLabel=1;			// 1: Write labels in output	  || 0: Don't write labels in output
	int status_reducedTime=0;		// 1: Reduced time increment    || 0: Full time increment
	int status_guessTrajectory=0;		// 1: Guess F on previous rates || 0: Homogenous forwarding of F
	int status_cutBack=0;			// 1: Reduce time increment     || 0: Keep current time increment
	int status_continueCalc=0;		// 1: Continues calculations    || 0: Stop with error
	int cutBack_level=0;			// The level of reduced time increment reduction
	int counter_convergedFull=0;		// Number of full time increments converged under the reduced time increment
	double timeInc_fraction;			// Fraction of a time step
	double time = 0.0;				// Simulation time
	
	// Pre-simulation calculations
	timeInc = timeTotal/N_steps;
	timeInc_old = timeInc;
	timeInc_fraction = 1.0;
	
	if(rank==0){
		printf(BOLDCYAN "\n------------------------\n");
		printf("Running Spectral Solver");
		printf(BOLDCYAN "\n------------------------\n" RESET);
	}
	
	// Simulation loop
	while(istep < N_steps){
		istep++;
		
		timeInc_fraction = 0.0;
		while(timeInc_fraction < 1.0){
			
			time += timeInc;
			timeInc_fraction += 1.0/pow(cutBack_factor, (double)cutBack_level); 
			
			if(rank ==0){
				if(timeInc_fraction == 1.0){
					printf(BOLDCYAN "\n\n\nStep %d: Time %e\n" RESET, istep, time);
				}
				else{
					printf(BOLDYELLOW "\n\n\nStep %3.3f: Time %e\t\tFull increments converged: %d\n", 
						 (double)istep+timeInc_fraction-1.0, time, counter_convergedFull);
					printf("Cutback factor: %f\t\tCutback level: %d\n" RESET, cutBack_factor, cutBack_level);
				}
			}
			
			// Spectral solver (the main course) ------------------------------------------
			switch(spectral_ID){
				case 1:
					spectral_forwardBASIC(istep, status_cutBack, status_guessTrajectory);
					status = spectral_solverBASIC(istep);
					break;
				case 2:
					if(rank == 0) printf(BOLDRED "Error: Augmented Lagrangian method not implemented yet.\n" RESET);
					return;
			}
			
			
			// If failed solution, set program to reduce time increment 
			// Otherwise, update hardening and do post processing for data.--------------------------------------
			status_cutBack = 0;
			if(status){
				if(cutBack_level < cutBack_lvlMax){
					if(rank==0) printf(BOLDYELLOW "Warning: Numerical instability detected. Re-evaluation with reduced time step.\n" RESET);
					status_cutBack = 1;
					counter_convergedFull = 0;
					status_reducedTime = 1;
					timeInc_fraction -= 1.0/pow(cutBack_factor, (double)cutBack_level);
					cutBack_level++;
					
					time -= timeInc;
					timeInc = timeInc/cutBack_factor;
				}
				else if(status_continueCalc){
					goto converged;
				}
				else{
					if(rank==0){
						printf(BOLDRED "\nError: Failed to achieve solution despite reduced time increments.");
						printf("\n       Spectral solver will now terminate. \n" RESET);
					}
					return;
				}
			}
			else{
				converged:
				// Record iteration variables
				status_guessTrajectory=1;
				timeInc_old = timeInc;
				
				// Increasing timeInc (if possible)
				if(status_reducedTime==1 && timeInc_fraction==1.0){
					counter_convergedFull++;
					
					if(counter_convergedFull == N_time2Increase){
						counter_convergedFull = 0;
						cutBack_level--;
						timeInc = timeInc*cutBack_factor;
						
						if(rank == 0) printf(BOLDGREEN "\nUpdate: Stable conditions detected. Attempting to increase time increment.");
						
						if(cutBack_level == 0){
							status_reducedTime = 0;
						}
					}
				}
				
				// Update hardening and 2nd PK stress-------------------------------------------
				MPI_Barrier(MPI_COMM_WORLD);
				local_loop{
					jph = findINT(phase_ID, node_phaseID[pIDX]);
					
					update_hardening(S[pIDX], Schmid0_xt[jph], crss_lastInc[pIDX], crss[pIDX], 
										gammaAccum_lastInc[pIDX], gammaAccum[pIDX], phase_NSYS[jph], jph);
					
					for(i=0;i<phase_NSYS[jph];i++){
						crss_lastInc[pIDX][i] = crss[pIDX][i];
						gammaAccum_lastInc[pIDX] = gammaAccum[pIDX];
					}
				}
				
				// Post processing for converged timestep---------------------------------------
				// 1. Field average information
				if (rank == 0){
					calc_vonMises();
					write_vonMises(istep, time, status_writeLabel);
					write_Favg(istep, time, status_writeLabel);
					write_Pavg(istep, time, status_writeLabel);
					status_writeLabel = 0;
				}
				
				// 2, Material point information
				if(istep % writeFreq == 0){
					writeP_MPI(istep);
					writeF_MPI(istep);
				}
			}	
		}
	}
	
	
	return;
	
}