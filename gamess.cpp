/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>
#include <unistd.h>

#define MAX_LENGTH 1024
#define MAX_SIZE 12474

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Write Gamess inputs for each state's FMO calculations 
-----------------------------------------------------------------*/
void State::write_gamess_inputs(int jobtype)
{
    // Writes a separate input file for all monomers and all dimers
    // Master rank does all the work here
    
    if (fmr->master_rank) {
        
        printf("Writing Gamess inputs.\n");
        printf("Read MOs: %d\n", flag_read_MOs);
        
        Run *run       = fmr->run;
        Atom *atom     = fmr->atom;
        int natoms     = atom->natoms;
        int nstates    = atom->nstates;
        int nfragments = atom->nfragments;
        
        int cellA      = fmr->atom->cellA;
        int cellB      = fmr->atom->cellB;
        int cellC      = fmr->atom->cellC;
        
        int xa         = fmr->atom->na;
        int xb         = fmr->atom->nb;
        int xc         = fmr->atom->nc;
        
        int afield     = fmr->atom->afield;
        int bfield     = fmr->atom->bfield;
        int cfield     = fmr->atom->cfield;
        
        // ***** Loop over states ***** //
        for (int istate=0; istate<nstates; ++istate) {
            
            // Determine the charged reactive fragment for this state
            int chgfrag = 0;
            for (int i=0; i<natoms; ++i) {
                if (atom->reactive[istate*natoms + i]) {
                    chgfrag = atom->fragment[istate*natoms + i];
                    break;
                }
            }
            
            // Put files in directory for organization
            char state_directory[256];
            char snum[16];
            char make_directory[256];
            
            sprintf(snum, "%02d", istate);
            sprintf(state_directory, "state_%02d", istate);
            // Make the directory...
            sprintf(make_directory, "mkdir -p %s", state_directory);
            int ierr = system(make_directory);
            
            char jobname[256];
            char filename[256];
            
            // Get name of job
            sprintf(jobname, "fmo_st%s", snum);
            
            sprintf(filename, "%s/%s.inp", state_directory, jobname);
            
            FILE *fs = fopen(filename, "w");
            if (fs == NULL) {
                char tmpstr[256];
                sprintf(tmpstr, "Failure to write Gamess input for file %s", filename);
                fmr->error(FLERR, tmpstr);
            }
            
            // Comment for labeling
            fprintf(fs, "! State %s\n", snum);
            
            // contrl section
            fprintf(fs, " $contrl scftyp=rhf runtyp=gradient $end\n");
	    int tmlm = nfragments * nfragments * 0.5 / 60;
	    if ( tmlm < 1 ) tmlm=1;
            fprintf(fs, " $system timlim=%d $end\n", tmlm);
            
            // scf section
            fprintf(fs, " $scf conv=1.0d-06 dirscf=.true. $end\n");
            
            // basis section
            fprintf(fs, " $basis  gbasis=sto ngauss=3 $end\n");
            
            // fmo section
            //fprintf(fs, " $fmo mplevl=2 nfrag=3 icharg(1)=1,0,0\n");
            fprintf(fs, " $fmo mplevl=2 nfrag=%d ", nfragments);
	    int count = 0;
            fprintf(fs, "icharg(1)=");
            for (int ifrag=0; ifrag<nfragments; ++ifrag) {
                if (ifrag == chgfrag) fprintf(fs, "1");
                else                  fprintf(fs, "0");
                
                if (ifrag == nfragments-1) { fprintf(fs, "\n"); break; }
                else                         fprintf(fs, ",");
		++count;
		if (count >= 10) { fprintf(fs, "\n                                       "); count = 0; } // keep frag charges from overunning

            }
            fprintf(fs, "      indat(1)=");
 	    count = 0;
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (iatom == natoms-1) { fprintf(fs, "%0d\n",atom->fragment[istate*natoms+iatom]+1); break; }
                else                     fprintf(fs, "%0d,",atom->fragment[istate*natoms+iatom]+1);
		++count;
		if (count >= 5) { fprintf(fs, "\n               "); count = 0; } // keep frag ids from overunning
            }
            //fprintf(fs, "      indat(1)=1,2,2,2,1,1,1,3,3,3\n");
            //fprintf(fs, "      indat(1)=");
            fprintf(fs, " $end\n");
            fprintf(fs, " $fmoprp nprint=0 $end\n");
            fprintf(fs, " $fmoxyz\n");
            
            // geometry section
            for (int iatom=0; iatom<natoms; ++iatom) {
                fprintf(fs, "%c %c %20.10lf %20.10lf %20.10lf\n",
                        atom->symbol[iatom],
                        atom->symbol[iatom],
                        atom->coord[3*iatom],
                        atom->coord[3*iatom+1],
                        atom->coord[3*iatom+2]
                        );
            }
            fprintf(fs, " $end\n");
            fprintf(fs, " $data\n");
            fprintf(fs, "Basis set input, with no atomic coordinates\n");
            fprintf(fs, "C1\n");
            fprintf(fs, "h-1 1\n");
            fprintf(fs, "o-1 8\n");
            fprintf(fs, " $end\n");
            
            fclose(fs);
            
        } // close loop over states
        
        printf("Done writing Gamess inputs.\n");
    }
    
    // Hold up
    MPI_Barrier(fmr->world);
}


/*-----------------------------------------------------------------
  Perform all the FMO calculations
-----------------------------------------------------------------*/
void Run::do_gamess_calculations(int FORCE)
{
    // Check if we should branch to the env approximation
    if (EnvApprox) {
        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
        //do_gamess_calculations_env();
        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
        return;
    }
    
    // Divide the list up across MPI ranks
    // Then, run the FMO calculations with call to serial Gamess in parallel
    
    
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
        printf("State FMO calculations: %d\n", nstates);
        
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
    int div = nstates / fmr->world_size;
    int rem = nstates % fmr->world_size;
    int ifrom_state, ito_state;
    if (my_rank < rem) {
        ifrom_state = my_rank*div + my_rank;
        ito_state   = ifrom_state + div + 1;
    } else {
        ifrom_state = my_rank*div + rem;
        ito_state   = ifrom_state + div;
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
    int nnodes=4;
    char verno[256];
    
    int index_state = 0;

    for (int istate=0; istate<nstates; ++istate) {
        if (fmr->master_rank) {
            
            char state_directory[256];
            char snum[16];
            sprintf(snum, "%02d", istate);
            sprintf(state_directory, "state_%02d", istate);
            
            char jobname[256];
            char filename[256];
            
            sprintf(jobname, "fmo_st%s", snum);
            
            // change directory
            char directory[512];
            sprintf(directory, "%s", state_directory);
            chdir(directory);
            
            //printf("Number of Gamess ncores: %d\n", gamess_ncores); 
            sprintf(command, "%s %s.inp %s %d > %s.log 2>&1",
                    exec,
                    jobname,
	            gamess_version,
                    gamess_ncores,
                    jobname
                    );
            printf("%s\n", command); 
 /*
            ierr = system("ls");
            ierr = system("/bin/ls");

	    ierr = system("echo dir=$PWD");
	    ierr = system("pwd");

            ierr = system("echo $PATH");
	    

	    ierr = system("cat rungms-xt");
 	    ierr = system("PATH=$PATH:. ; rungms-xt");

	    ierr = system("sh -c 'echo hi'");
            ierr = system("./sys");
            ierr = system("./gamess.asis.x");
*/

	    // ** The system call ** //
            ierr = system(command);
            
            // ** Check for error ** //
            if (ierr) {
                printf("Gamess run error on rank %d:\n", fmr->my_rank);
                fmr->error(FLERR, command);
            }
            
            // ** Open output file and get the energy ** //
            char output_file[MAX_LENGTH];
            sprintf(output_file, "%s.log",  jobname);
            FILE *fs = fopen(output_file, "r");
            if (fs == NULL) {
                char tmpstr[MAX_LENGTH];
                sprintf(tmpstr, "Failure to read Gamess output file: %s", output_file);
                fmr->error(FLERR, tmpstr);
            }
            char line[MAX_LENGTH];
            double en;
            double gx, gy, gz;
            
            char tmp0[16],tmp1[16],tmp2[16],tmp3[16];
            
            
            while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                //printf("%s",line);
                if ( strstr(line, "Two-body FMO properties") ) {
                    
                    while ( fgets(line, MAX_LENGTH, fs) != NULL) {
                        
                        if ( strstr(line, "ATOM# FRG#") ) {
                            
                            for (int iatom=0; iatom<natoms; ++iatom) {
                                
                                fgets(line, MAX_LENGTH, fs);
                                //printf("%s",line);
                                if ( sscanf(line, "%s %s %s %lf %lf %lf", tmp0, tmp1, tmp2, &gx, &gy, &gz) == 6 ) {
                                    
                                    fmo_gradients[3*natoms*istate + 3*iatom]   = gx;
                                    fmo_gradients[3*natoms*istate + 3*iatom+1] = gy;
                                    fmo_gradients[3*natoms*istate + 3*iatom+2] = gz;
                                    
                                }
                            }
                        }
                        else if ( strstr(line, "TOTAL ENERGY =") ) {
                            
                            //printf("%s",line);
                            if ( sscanf(line, "%s %s %s %lf", tmp0, tmp1, tmp2, &en) == 4 ) {
                                fmo_energies[istate] = en;
                            }
                        }
                    }
                    break;
                }
            }
        fclose(fs);
        

        // clean up directory
        sprintf(command, "rm -rf %s/%s/%s.dat",scratch_dir,state_directory,jobname);
	
        
        // ** The system call ** //
        ierr = system(command);
        
        // ** Check for error ** //
        if (ierr) {
            printf("Gamess run error on rank %d:\n", fmr->my_rank);
            fmr->error(FLERR, command);
        }
        
        chdir("../");
        }
    }
    
    // Clock
    double FMO_end = MPI_Wtime();
    MPI_Barrier(fmr->world);
    if (fmr->master_rank) {
        printf("Finished with FMO calculations.\n\n");
    }
    
    // *** Reduce the energies from the parallel Gamess calls *** //
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
    
    /*
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
    //fprintf(fs,"\n");
    //fclose(fs);
    */
    // Broadcast the FMO energy/force to worker ranks
    MPI_Bcast(fmo_energies, nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);
    if (FORCE) {
        MPI_Bcast(fmo_gradients, nstates*3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
    }
    
#ifdef FMR_DEBUG
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
#endif
    
    if (fmr->master_rank) {
        printf(" Time for FMO: %.4f seconds\n\n", FMO_end - FMO_start);
        //printf(" Proc %d: Time for FMO: %.4f seconds\n", my_rank, FMO_end - FMO_start);
    }
}
  
