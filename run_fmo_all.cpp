/* AWGL */
#include "fmr.h"
#include "run.h"
#include "atom.h"
#include "state.h"
#include "input.h"

#define MAX_LENGTH 1024
#define MAX_SIZE 12474

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Perform all the FMO calculations 
-----------------------------------------------------------------*/
void Run::do_qchem_calculations(int FORCE)
{
  // Check if we should branch to the env approximation
  if (EnvApprox) {
    do_qchem_calculations_env();
    return;
  }

  // Divide the list up across MPI ranks
  // Then, run the FMO calculations with call to serial Q-Chem in parallel

  int natoms     = fmr->atom->natoms;
  int nfragments = fmr->atom->nfragments;
  int nstates    = fmr->atom->nstates;
  int my_rank    = fmr->my_rank;
  int nf2        = nfragments*nfragments;

  int na	 = 2*fmr->atom->na + 1;
  int nb 	 = 2*fmr->atom->nb + 1;
  int nc 	 = 2*fmr->atom->nc + 1;

  int xa	 = fmr->atom->na;
  int xb	 = fmr->atom->nb;
  int xc 	 = fmr->atom->nc;

  int nafield	 = 2*fmr->atom->afield + 1;
  int nbfield	 = 2*fmr->atom->bfield + 1;
  int ncfield	 = 2*fmr->atom->cfield + 1; 

  int afield     = fmr->atom->afield;
  int bfield     = fmr->atom->bfield;
  int cfield     = fmr->atom->cfield;

  n_monomers = nstates * nfragments * na*nb*nc;
  // Assuming all states have equal number of dimers, for now
  n_dimers = nstates * (nf2 * (na*nb*nc-1) + (nfragments * (nfragments-1)) / 2);
  n_dimers_sq = nstates * nf2 *na*nb*nc; // inclues self

  if (fmr->master_rank) {
    printf("Preparing to run FMO calculations:\n");
    printf("Monomer FMO calculations: %d\n", n_monomers);
    printf("Dimer FMO calculations:   %d\n", n_dimers);
  }

  // ** If number of states changed, need to deallocate memory and reallocate below ** //
  if (fmr->state->flag_state_number_change) {
    delete [] fmo_energies;
    delete [] monomer_energies;
    delete [] dimer_energies;
    fmo_energies = monomer_energies = dimer_energies = NULL;
    if (FORCE) {
      delete [] fmo_gradients;
      delete [] monomer_gradients;
      delete [] dimer_gradients;
      fmo_gradients = monomer_gradients = dimer_gradients = NULL;
    }
  }

  // ** Allocate energies and zero ** //
  if (fmo_energies == NULL)     fmo_energies     = new double [nstates];
  if (monomer_energies == NULL) monomer_energies = new double [n_monomers];
  if (dimer_energies == NULL)   dimer_energies   = new double [n_dimers_sq];
  for (int i=0; i<nstates; ++i)     fmo_energies[i]     = 0.0;
  for (int i=0; i<n_monomers; ++i)  monomer_energies[i] = 0.0;
  for (int i=0; i<n_dimers_sq; ++i) dimer_energies[i]   = 0.0;

  if (FORCE) {
    // ** Allocate gradients and zero ** //
    if (fmo_gradients == NULL)     fmo_gradients     = new double[nstates*3*natoms];
    if (monomer_gradients == NULL) monomer_gradients = new double[n_monomers*3*natoms];
    if (dimer_gradients == NULL)   dimer_gradients   = new double[n_dimers_sq*3*natoms];
    for (int i=0; i<nstates*3*natoms; ++i)     fmo_gradients[i]     = 0.0;
    for (int i=0; i<n_monomers*3*natoms; ++i)  monomer_gradients[i] = 0.0;
    for (int i=0; i<n_dimers_sq*3*natoms; ++i) dimer_gradients[i]   = 0.0;
  }

  // ** Determine load balance ** //
  int div = n_monomers / fmr->world_size;
  int rem = n_monomers % fmr->world_size;
  int ifrom_mono, ito_mono;
  if (my_rank < rem) {
    ifrom_mono = my_rank*div + my_rank;
    ito_mono   = ifrom_mono + div + 1;
  } else {  
    ifrom_mono = my_rank*div + rem;
    ito_mono   = ifrom_mono + div;
  }
  div = n_dimers / fmr->world_size;
  rem = n_dimers % fmr->world_size;
  int ifrom_dim, ito_dim;
  if (my_rank < rem) {
    ifrom_dim = my_rank*div + my_rank;
    ito_dim   = ifrom_dim + div + 1;
  } else {  
    ifrom_dim = my_rank*div + rem;
    ito_dim   = ifrom_dim + div;
  }

  // Clock
  MPI_Barrier(fmr->world);
  double FMO_start = MPI_Wtime();
  if (fmr->master_rank) {
    printf("Running FMO calculations...\n");
  }

  // ** Make the system calls to run each FMO calculation now ** //
  char command[MAX_LENGTH];
  int ierr;
  int index_mono = 0;
  int index_dim  = 0;

  for (int istate=0; istate<nstates; ++istate) {

    char state_directory[256];
    char snum[16];
    sprintf(snum, "%02d", istate);
    sprintf(state_directory, "state_%02d", istate);

    // ********** Handle FMO monomers here ************ //
    for (int x=-xa; x<=xa; x++) {
      for (int y=-xb; y<=xb; y++) {
        for (int z=-xc; z<=xc; z++) {

          //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
          //if (abs(x)!=abs(y) || abs(y)!=abs(z) || abs(x)!=abs(z)) continue;
          //if (x==abs(xa) && y==abs(xb) && z==abs(xc)) continue;
          //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG

	  for (int ifrag=0; ifrag<nfragments; ++ifrag) {
            if (ifrom_mono <= index_mono && index_mono < ito_mono) {

              char filename[256];
              char inum[16];

	      sprintf(inum,"%03d",ifrag);

	      char cname[16];
	      sprintf(cname,"cell.%d.%d.%d", x+xa, y+xb, z+xc);

              sprintf(filename, "fmo_st%s_m%s_%s", snum, inum, cname);
              sprintf(command, "%s %s/%s.in %s/%s/ > %s/%s.out", 
                      exec,
                      state_directory,
                      filename,
                      scratch_dir,
                      filename,
                      state_directory,
                      filename
                     );
              //printf("Rank %d: %s\n", my_rank, command);

              // ** The system call ** //
              ierr = system(command);

              // ** Check for error ** //
              if (ierr) {
                printf("Q-Chem run error on rank %d:\n", fmr->my_rank);
                fmr->error(FLERR, command);
              }

              // ** Open output file and get the energy ** //
              char output_file[MAX_LENGTH];
              sprintf(output_file, "%s/%s.in.energy", state_directory, filename);
              FILE *fs = fopen(output_file, "r");
              if (fs == NULL) {
                char tmpstr[MAX_LENGTH];
                sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
                fmr->error(FLERR, tmpstr);
              }
              char line[MAX_LENGTH];
              double en;

	      
              while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
      	        if ( sscanf(line, "%lf", &en) == 1 ) {

		  //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
		  monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag] = en;
		  //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
	        }
              }
              fclose(fs);
	      //printf("Rank %d: energy\n", my_rank);

              if (FORCE) {
                // ** Get gradient from file ** // 
                sprintf(output_file, "%s/%s.in.gradient", state_directory, filename);
                fs = fopen(output_file, "r");
                if (fs == NULL) {
                  char tmpstr[MAX_LENGTH];
                  sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
                  fmr->error(FLERR, tmpstr);
                }
                char line[MAX_LENGTH];
                int iatom, atnum;
                atnum = 0; // index of QM atom for storing gradient
                double gx, gy, gz;
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                  // Advance atnum until it matches as a QM atom index for this monomer fragment
                  while ( !fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z) ) {
                    atnum++;
                  }
	          if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {
                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                    monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx;                    
                    monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                    monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
		    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
	          }
                  // Increment atnum for the next round
                  atnum++;
                }
                fclose(fs);
		//printf("Rank %d: grad\n", my_rank);


                // ** Get field from file ** // 
                sprintf(output_file, "%s/%s.in.field", state_directory, filename);
                fs = fopen(output_file, "r");
                if (fs == NULL) {
                  char tmpstr[MAX_LENGTH];
                  sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
                  fmr->error(FLERR, tmpstr);
                }
                atnum = 0; // index of non-QM atom for storing gradient
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                  // Advance atnum until it matches as a non-QM atom index for this monomer fragment
                  while ( fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z, afield, bfield, cfield) ) {
		  //while ( fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z) ) {
                    atnum++; 
                  }
		  
		  //printf("atnum:%d ifrag:%d istate:%d istate:%d x:%d y:%d z:%d iatom:%d gx:%d gy:%d gz:%d \n",atnum,ifrag,istate,x,y,z,iatom,gx,gy,gz);
	
  	          if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {

		    if (fmr->atom->AtomInCell(atnum,istate,0,0,0, afield, bfield, cfield)) {	
                      // gx,gy,gz = the electric field
                      // multiply by charge to get force (i.e. negative gradient) on atom
                      double mmq = fmr->atom->getCharge(atnum%natoms, istate);
                      gx *= -mmq;
                      gy *= -mmq;
                      gz *= -mmq;
	  	      //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                      monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx;
  		      monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
		      monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                      //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
		    }
		  }
                  // Increment atnum for the next round
                  atnum++;
                }
                fclose(fs);
		//printf("Rank %d: efield\n", my_rank);
	
		
              }

            }
            ++index_mono;
          }

	}
      }
    }

    // ********** Handle FMO dimers in this loop ************ //
    for (int x=-xa; x<=xa; x++) {
      for (int y=-xb; y<=xb; y++) {
        for (int z=-xc; z<=xc; z++) {

	  //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
          //if (x!=abs(xa) && y==abs(xb) && z==abs(xc)) continue;
	  //if (abs(x)!=abs(y) || abs(y)!=abs(z) || abs(x)!=abs(z)) continue;
	  //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG

          for (int ifrag=0; ifrag<nfragments; ++ifrag) {
            for (int jfrag=0; jfrag<nfragments; ++jfrag) {
	
	      if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;

              if (ifrom_dim <= index_dim && index_dim < ito_dim) {

  	        char filename[256];
                //char inum[16];
                //char jnum[16];

                sprintf(filename, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
   	        sprintf(command, "%s %s/%s.in %s/%s/ > %s/%s.out", 
		        exec,
                        state_directory,
   		        filename,
		        scratch_dir,
		        filename,
                        state_directory,
		        filename
		       );
	        //printf("Rank %d: %s\n", my_rank, command);

                // ** The system call ** //
	        ierr = system(command);

                // ** Check for error ** //
                if (ierr) {
                  printf("Q-Chem run error on rank %d:\n", fmr->my_rank);
                  fmr->error(FLERR, command);
                }

	        // ** Open output file and get the energy ** //
                char output_file[MAX_LENGTH];
                //sprintf(output_file, "%s/%s.out", state_directory, filename);
                sprintf(output_file, "%s/%s.in.energy", state_directory, filename);
	        FILE *fs = fopen(output_file, "r");
	        if (fs == NULL) {
	          char tmpstr[MAX_LENGTH];
	          sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
	          fmr->error(FLERR, tmpstr);
	        }
	        char line[MAX_LENGTH];
	        double en=0.0;
	        while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
	          if ( sscanf(line, "%lf", &en) == 1 ) {
                    // save symmetrized

                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
		    // save symmetrized for zeroth unit cell
		    if (x==0 && y==0 && z==0) {
                      dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag] = en;
		      dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag] = en;
		    } else {
		      dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag] = en;
		    }
                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG

	          }
	        }
	        fclose(fs);

                if (FORCE) {
                  // ** Get gradient from file ** // 
                  sprintf(output_file, "%s/%s.in.gradient", state_directory, filename);
                  fs = fopen(output_file, "r");
                  if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                  }
                  char line[MAX_LENGTH];
                  int iatom; // dummy index
                  int atnum = 0; // index of QM atom for storing gradient
                  double gx, gy, gz;
                  while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    // Advance atnum until it matches as a QM atom index for this dimer fragment
                    while ( !(fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0) || 
                              fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z)) ) {
                      atnum++; 
                    }

  	            if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {

		      if (fmr->atom->AtomInCell(atnum,istate,0,0,0)) {
		      //if (fmr->atom->AtomInFragment(atnum, jfrag, istate, 0, 0, 0) || fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0)) {		

                        // store symmetrically for zeroth unit cell
		        if (x==0 && y==0 && z==0) {
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz; 
	  	        } else {
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
			}
		      }
	            }
                    // Increment atnum for the next round
                    atnum++;
                  }
                  fclose(fs);

                  // ** Get field from file ** // 
                  sprintf(output_file, "%s/%s.in.field", state_directory, filename);
                  fs = fopen(output_file, "r");
                  if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read Q-Chem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                  }
                  atnum = 0; // index of non-QM atom for storing gradient
                  while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    // Advance atnum until it matches as a non-QM atom index for this dimer fragment
                    //while ( (fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0) ||
                    //         fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z)) ) {
                    while ( (fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0, afield, bfield, cfield) || 
                             fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z, afield, bfield, cfield)) ) {
                      atnum++; 
                    }

 	            if ( sscanf(line, "%d %lf %lf %lf", &iatom, &gx, &gy, &gz) == 4 ) {

		      if (fmr->atom->AtomInCell(atnum,istate,0,0,0,afield,bfield,cfield)) {
                        // gx,gy,gz = the electric field
                        // multiply by charge to get force (i.e. negative gradient) on atom
                        double mmq = fmr->atom->getCharge(atnum%natoms, istate);
                        gx *= -mmq;
                        gy *= -mmq;
                        gz *= -mmq;
                        // store symmetrically for zeroth unit cell
                        if (x==0 && y==0 && z==0) { 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy; 
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz; 
	  	        } else {
			  dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                          dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
		        }
                      }
		    }
                    // Increment atnum for the next round
                    atnum++;
                  }
                  fclose(fs);

                } 
	      }
              ++index_dim;
            }
          }

	}
      }
    }
  }

  // Clock
  double FMO_end = MPI_Wtime();
  MPI_Barrier(fmr->world);
  if (fmr->master_rank) {
    printf("Finished with FMO calculations.\n\n");
  }

  // *** Reduce the energies from the parallel Q-Chem calls *** //
  double *rbuffer;
  //if (FORCE) rbuffer = new double [MAX_SIZE];
  //else       rbuffer = new double [MAX_SIZE];
  if (FORCE) rbuffer = new double [n_dimers_sq*3*natoms]; 
  else       rbuffer = new double [n_dimers_sq];

  // Monomers
  for (int i=0; i<n_monomers; ++i) rbuffer[i] = 0.0;
  MPI_Allreduce(monomer_energies, rbuffer, n_monomers, MPI_DOUBLE, MPI_SUM, fmr->world);
  for (int i=0; i<n_monomers; ++i) monomer_energies[i] = rbuffer[i];
  // Dimers
  for (int i=0; i<n_dimers_sq; ++i) rbuffer[i] = 0.0;
  MPI_Allreduce(dimer_energies, rbuffer, n_dimers_sq, MPI_DOUBLE, MPI_SUM, fmr->world);
  for (int i=0; i<n_dimers_sq; ++i) dimer_energies[i] = rbuffer[i];
  
  if (FORCE) {
    // Monomers
    for (int i=0; i<n_monomers*3*natoms; ++i) rbuffer[i] = 0.0;
    MPI_Allreduce(monomer_gradients, rbuffer, n_monomers*3*natoms, MPI_DOUBLE, MPI_SUM, fmr->world);
    for (int i=0; i<n_monomers*3*natoms; ++i) monomer_gradients[i] = rbuffer[i];
    // Dimers
    for (int i=0; i<n_dimers_sq*3*natoms; ++i) rbuffer[i] = 0.0;
    MPI_Allreduce(dimer_gradients, rbuffer, n_dimers_sq*3*natoms, MPI_DOUBLE, MPI_SUM, fmr->world);
    for (int i=0; i<n_dimers_sq*3*natoms; ++i) dimer_gradients[i] = rbuffer[i];
  }

  delete [] rbuffer;

  // *** Compute the FMO energies/forces for each state *** //
  
  
  FILE *fs = fopen("fmr_calc.log", "a");
    
  if (fmr->master_rank) {
    for (int istate=0; istate<nstates; ++istate) {
      printf("----- State %4d -----\n", istate);
      if (fmr->print_level > 0) {
        fprintf(fs,"Monomer | EI\n");
      }
      double en_fmo1 = 0.0;
      for (int x=-xa; x<=xa; x++) {
        for (int y=-xb; y<=xb; y++) {
          for (int z=-xc; z<=xc; z++) {

	    if (x==0 && y==0 && z==0) {
              for (int ifrag=0; ifrag<nfragments; ++ifrag) {
                if (fmr->print_level > 0) {
 		  fprintf(fs,"FRAG I:%4d | %16.10f\n", ifrag, monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag]);
                }
                en_fmo1 += monomer_energies[nfragments*istate + ifrag];
	      }
            }

	  }
	}
      }
      printf("Total monomer energy: %16.10f\n\n", en_fmo1);
      if (fmr->print_level > 0) {
        fprintf(fs,"Dimer | EIJ EI EJ (EIJ - EI - EJ)\n");
      }
      double en_fmo2 = 0.0;

      for (int x=-xa; x<=xa; x++) {
        for (int y=-xb; y<=xb; y++) {
          for (int z=-xc; z<=xc; z++) {

            for (int ifrag=0; ifrag<nfragments; ++ifrag) {
              for (int jfrag=0; jfrag<nfragments; ++jfrag) {

		if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;
		
		//overcount for unit cell 
		double oc=0.5; if (x==0 && y==0 && z==0) oc=1.0;

                double etmp = dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag] - 
                              monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(xa) +   nc*nfragments*(xb) +   nfragments*(xc) +   ifrag] -
                              monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag];
                if (fmr->print_level > 0) {
                  fprintf(fs,"FRAG I:%02d--J:%02d  CELL x:%2d y:%2d z:%2d | %16.10f %16.10f %16.10f %16.10f\n", ifrag, jfrag,x,y,z,
			  dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag],                            
                          monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(xa) +   nc*nfragments*(xb) +   nfragments*(xc) +   ifrag],
                          monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag], etmp
                          //dimer_energies[nfragments*nfragments*istate + nfragments*ifrag + jfrag],
                          //monomer_energies[nfragments*istate + ifrag],
                          //monomer_energies[nfragments*istate + jfrag], etmp
                        );
                }
                en_fmo2 += etmp*oc;
              }
            }

	  }
        }
      }

      printf("Total dimer energy:   %16.10f\n", en_fmo2);
      fmo_energies[istate] = en_fmo1 + en_fmo2;
      printf("Total FMO energy:     %16.10f\n", en_fmo1 + en_fmo2);
      printf("----------------------\n");

      // ** Do forces for this state ** //
      if (FORCE) {
        double gx, gy, gz;
        gx = gy = gz = 0.0;
        // Monomer force
        for (int x=-xa; x<=xa; x++) {
          for (int y=-xb; y<=xb; y++) {
            for (int z=-xc; z<=xc; z++) {

	      if (x==0 && y==0 && z==0) {
                for (int ifrag=0; ifrag<nfragments; ++ifrag) {
                  for (int i=0; i<natoms; ++i) {
                    gx = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i];
                    gy = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i+1];
                    gz = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i+2];
                    fmo_gradients[3*natoms*istate + 3*i]   += gx;
                    fmo_gradients[3*natoms*istate + 3*i+1] += gy;
                    fmo_gradients[3*natoms*istate + 3*i+2] += gz;
                  }
                }
	      }

	    }
	  }
        }
        // Dimer force
        //
        for (int x=-xa; x<=xa; x++) { 
          for (int y=-xb; y<=xb; y++) { 
            for (int z=-xc; z<=xc; z++) {

              for (int ifrag=0; ifrag<nfragments; ++ifrag) {
                for (int jfrag=0; jfrag<nfragments; ++jfrag) {

		  if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;

  		  //overcount for unit cell 
  		  double oc=1.0; if (x==0 && y==0 && z==0) oc=1.0;

                  for (int i=0; i<natoms; ++i) {
                    double tmpx, tmpy, tmpz;

                    //double dx = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i];
                    //double dy = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+1];
                    //double dz = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+2];

                    gx = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i];  
                    gy = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+1];
                    gz = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+2];

		    //double m1x = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i];
		    //double m1y = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i+1];
                    //double m1z = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i+2];

                    //double m2x = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i];
                    //double m2y = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i+1];
                    //double m2z = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i+2];     

                    gx -= monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i] + 
                          monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i];
                    gy -= monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i+1] + 
                          monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i+1];
                    gz -= monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i+2] + 
                          monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i+2];

                    fmo_gradients[3*natoms*istate + 3*i]   += gx*oc;
                    fmo_gradients[3*natoms*istate + 3*i+1] += gy*oc;
                    fmo_gradients[3*natoms*istate + 3*i+2] += gz*oc;

		    //fprintf(fs,"%15.10f,%15.10f,%15.10f,%15.10f,\n",dx,m1x,m2x,gx);
		    //fprintf(fs,"%15.10f,%15.10f,%15.10f,%15.10f,\n",dy,m1y,m2y,gy);
		    //fprintf(fs,"%15.10f,%15.10f,%15.10f,%15.10f,\n",dz,m1z,m2z,gz);

		  }
		}
	      }  

            }
          }
        }
      }
    }
  }
  fprintf(fs,"\n");
  fclose(fs);
  // Broadcast the FMO energy/force to worker ranks
  MPI_Bcast(fmo_energies, nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);
  if (FORCE) {
    MPI_Bcast(fmo_gradients, nstates*3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
  }

//#ifdef FMR_DEBUG
  if (fmr->master_rank && FORCE) {
    for (int x=-xa; x<=xa; x++) {
      for (int y=-xb; y<=xb; y++) {
        for (int z=-xc; z<=xc; z++) {
          printf("Monomer gradients:\n");
          for (int istate=0; istate<nstates; ++istate) {
            for (int ifrag=0; ifrag<nfragments; ++ifrag) {
              printf("State %d fragment %d cell %d %d %d:\n", istate, ifrag, x, y, z);
              for (int i=0; i<natoms; ++i) {
                printf("%3d %12.8f %12.8f %12.8f\n", i, 
		       monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i],
		       monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i+1],
		       monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i+2]

                       //monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i],
                       //monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+1],
                       //monomer_gradients[(nfragments*istate + ifrag)*3*natoms + 3*i+2]
                      );
	      }
	    }
	  }
        }
      }
    }
    printf("Dimer gradients:\n");
    for (int x=-xa; x<=xa; x++) {
      for (int y=-xb; y<=xb; y++) {
        for (int z=-xc; z<=xc; z++) {

          for (int istate=0; istate<nstates; ++istate) {
            for (int ifrag=0; ifrag<nfragments; ++ifrag) {
              for (int jfrag=0; jfrag<nfragments; ++jfrag) {

		if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;

                printf("State %d fragments %d--%d cell %d %d %d:\n", istate, ifrag, jfrag, x, y, z);
                for (int i=0; i<natoms; ++i) {
                  printf("%3d %12.8f %12.8f %12.8f\n", i, 
			 dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i],
			 dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+1],
			 dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+2]

                         //dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i],  
                         //dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i+1],
                         //dimer_gradients[(nf2*istate + nfragments*ifrag + jfrag)*3*natoms + 3*i+2]
                        );
	        }
	      }
	    }
          }

        }
      }
    }
    printf("FMO gradients:\n");
    for (int istate=0; istate<nstates; ++istate) {
      printf("State %d:\n", istate);
      for (int i=0; i<natoms; ++i) {
        double gx, gy, gz;
        gx = fmo_gradients[3*natoms*istate + 3*i];   
        gy = fmo_gradients[3*natoms*istate + 3*i+1]; 
        gz = fmo_gradients[3*natoms*istate + 3*i+2];
        printf("%d : %f %f %f\n", i, gx, gy, gz);
      }
    }
  }
//#endif

  if (fmr->master_rank) {
    printf(" Time for FMO: %.4f seconds\n\n", FMO_end - FMO_start);
    //printf(" Proc %d: Time for FMO: %.4f seconds\n", my_rank, FMO_end - FMO_start);
  }
}

