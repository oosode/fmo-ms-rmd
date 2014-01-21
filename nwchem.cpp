/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Write NWChem inputs for each state's FMO calculations 
-----------------------------------------------------------------*/
void State::write_nwchem_inputs(int jobtype)
{
  // Writes a separate input file for all monomers and all dimers
  // Master rank does all the work here

  if (fmr->master_rank) {

    printf("Writing NWChem inputs.\n");
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
        
        // *** Monomers *** //
        for (int x=-xa; x<=xa; ++x) {
            for (int y=-xb; y<=xb; ++y) {
                for (int z=-xc; z<=xc; ++z) {
                    
                    for (int ifrag=0; ifrag<nfragments; ++ifrag){
                        // Get name of file to open
                        char filename[256];
                        sprintf(filename, "%s/fmo_st%s_m%03d_cell.%d.%d.%d.nw", state_directory, snum, ifrag, x+xa, y+xb, z+xc);
                        
                        FILE *fs = fopen(filename, "w");
                        if (fs == NULL) {
                            char tmpstr[256];
                            sprintf(tmpstr, "Failure to write NWChem input for file %s", filename);
                            fmr->error(FLERR, tmpstr);
                        }
                        
                        // Comment for labeling
                        fprintf(fs, "start energy_gradients\n");
                        fprintf(fs, "title 'State %d Monomer %d Cell %d %d %d'\n\n", istate, ifrag, x+xa, y+xb, z+xc);
                        
                        // geometry section
                        fprintf(fs, "geometry nocenter noautoz units angstrom\n");
                        fprintf(fs, "symmetry c1\n");
                        for (int iatom=0; iatom<natoms; ++iatom) {
                            if (atom->fragment[istate*natoms + iatom] == ifrag) {
                                fprintf(fs, "%c %3.1lf %20.10lf %20.10lf %20.10lf\n",
                                        atom->symbol[iatom],
                                        atom->getAtomicNumber(iatom),
                                        atom->charge[iatom],
                                        atom->coord[3*iatom] + x*cellA,
                                        atom->coord[3*iatom+1] + y*cellB,
                                        atom->coord[3*iatom+2] + z*cellC
                                        );
                            }
                        }
                        fprintf(fs, "end\n\n");
                        
                        
                        // $GUESS section
                        // Read previous step MO coeffs?
                        //if (flag_read_MOs) {
                        //    fprintf(fs, " $GUESS ");
                        //    fprintf(fs, "GUESS=MOREAD ");
                        //    fprintf(fs, "$END\n");
                        //}
                        
                        // bq section
                        fprintf(fs, "bq units angstrom\n");
                        for (int x0=-afield; x0<=afield; ++x0) {
                            for (int y0=-bfield; y0<=bfield; ++y0) {
                                for (int z0=-cfield; z0<=cfield; ++z0) {
                                    
                                    for (int iatom=0; iatom<natoms; ++iatom) {
                                        if (atom->fragment[istate*natoms + iatom] != ifrag || x0 != x || y0 != y || z0 != z) {
                                            double mmq = atom->getCharge(iatom, istate);
                                            fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                                                    atom->coord[3*iatom] + x0*cellA,
                                                    atom->coord[3*iatom+1] + y0*cellB,
                                                    atom->coord[3*iatom+2] + z0*cellC,
                                                    mmq
                                                    );
                                        }
                                    }
                                    
                                }
                            }
                        }
                        fprintf(fs, "end\n\n");
                        
                        // scf section
                        fprintf(fs, "scf\n");
                        fprintf(fs, "singlet\n");
                        fprintf(fs, "direct\n");
                        fprintf(fs, "thresh 1e-6\n");
                        fprintf(fs, "sym off\n");
                        fprintf(fs, "end\n\n");
                        
                        // correlation section
                        fprintf(fs, "mp2\n");
                        fprintf(fs, "freeze core atomic\n");
                        fprintf(fs, "end\n\n");
                        
                        // charge section
                        if (ifrag == chgfrag) {
                            fprintf(fs, "charge 1\n\n");
                        } else {
                            fprintf(fs, "charge 0\n\n");
                        }
                       
			// scratch section
			char fname[256];
			char scratch[256];
			sprintf(fname, "fmo_st%s_m%03d_cell.%d.%d.%d.nw", snum, ifrag, x+xa, y+xb, z+xc);
			sprintf(scratch, "%s/%s/",run->scratch_dir,fname);	
			fprintf(fs, "scratch_dir %s\n\n",scratch);
 
                        // basis set section
                        fprintf(fs, "basis\n");
                        fprintf(fs, "* library %s\n", run->basis);
                        fprintf(fs, "end\n\n");
                        
                        // task section
                        if (jobtype == RUN_ENERGY)
                            fprintf(fs, "task %s %s\n\n", run->correlation, "energy");
                        else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                            fprintf(fs, "task %s %s\n\n", run->correlation, "gradient");
                        
                        
                        /*
                         *
                         */
                        
                        // Comment for labeling
                        fprintf(fs, "\n\nstart field\n");
                
                        // geometry section
                        fprintf(fs, "geometry nocenter noautoz units angstrom\n");
                        fprintf(fs, "symmetry c1\n");
                        for (int iatom=0; iatom<natoms; ++iatom) {
                            if (atom->fragment[istate*natoms + iatom] == ifrag) {
                                fprintf(fs, "%c %3.1lf %20.10lf %20.10lf %20.10lf\n",
                                        atom->symbol[iatom],
                                        atom->getAtomicNumber(iatom),
                                        atom->charge[iatom],
                                        atom->coord[3*iatom] + x*cellA,
                                        atom->coord[3*iatom+1] + y*cellB,
                                        atom->coord[3*iatom+2] + z*cellC
                                        );
                            }
                        }
                        fprintf(fs, "end\n\n");
                        
                        
                        // $GUESS section
                        // Read previous step MO coeffs?
                        //if (flag_read_MOs) {
                        //    fprintf(fs, " $GUESS ");
                        //    fprintf(fs, "GUESS=MOREAD ");
                        //    fprintf(fs, "$END\n");
                        //}
                        
                        // bq section
                        fprintf(fs, "bq units angstrom\n");
                        fprintf(fs, "force fmo_st%s_m%03d_cell.%d.%d.%d.nw.field\n", snum, ifrag, x+xa, y+xb, z+xc);
                        for (int x0=-afield; x0<=afield; ++x0) {
                            for (int y0=-bfield; y0<=bfield; ++y0) {
                                for (int z0=-cfield; z0<=cfield; ++z0) {
                                    
                                    for (int iatom=0; iatom<natoms; ++iatom) {
                                        if (atom->fragment[istate*natoms + iatom] != ifrag || x0 != x || y0 != y || z0 != z) {
                                            double mmq = atom->getCharge(iatom, istate);
                                            fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                                                    atom->coord[3*iatom] + x0*cellA,
                                                    atom->coord[3*iatom+1] + y0*cellB,
                                                    atom->coord[3*iatom+2] + z0*cellC,
                                                    mmq
                                                    );
                                        }
                                    }
                                    
                                }
                            }
                        }
                        fprintf(fs, "end\n\n");
                        
                        // scf section
                        fprintf(fs, "scf\n");
                        fprintf(fs, "singlet\n");
                        fprintf(fs, "direct\n");
                        fprintf(fs, "thresh 1e-6\n");
                        fprintf(fs, "sym off\n");
                        fprintf(fs, "end\n\n");
                        
                        // charge section
                        if (ifrag == chgfrag) {
                            fprintf(fs, "charge 1\n\n");
                        } else {
                            fprintf(fs, "charge 0\n\n");
                        }
                        
			// scratch section
			fprintf(fs, "scratch_dir %s\n\n",scratch); 

                        // basis set section
                        fprintf(fs, "basis\n");
                        fprintf(fs, "* library %s\n", run->basis);
                        fprintf(fs, "end\n\n");
                        
                        // task section
                        if (jobtype == RUN_ENERGY)
                            fprintf(fs, "task scf %s\n\n", run->correlation, "energy");
                        else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                            fprintf(fs, "task scf %s\n\n", run->correlation, "gradient");

                        
                        fclose(fs);
                    } // close loop over fragments for monomers
                    
                }
            }
        }
        
        // *** Dimers *** //
        for (int x=-xa; x<=xa; ++x) {
            for (int y=-xb; y<=xb; ++y) {
                for (int z=-xc; z<=xc; ++z) {
                    
                    for (int ifrag=0; ifrag<nfragments; ++ifrag){
                        for (int jfrag=0; jfrag<nfragments; ++jfrag){
                            
                            if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;
                            
                            
                            // Get name of file to open
                            char filename[256];
                            char inum[16];
                            char jnum[16];
                            
                            sprintf(filename, "%s/fmo_st%s_d%03d-%03d_cell.%d.%d.%d.nw", state_directory, snum, ifrag, jfrag, x+xa, y+xb, z+xc);
                            
                            FILE *fs = fopen(filename, "w");
                            if (fs == NULL) {
                                char tmpstr[256];
                                sprintf(tmpstr, "Failure to write NWChem input for file %s", filename);
                                fmr->error(FLERR, tmpstr);
                            }
                            
                            // Comment for labeling
                            fprintf(fs, "start energy_gradients\n");
                            fprintf(fs, "title 'State %d Monomer %d Cell %d %d %d'\n\n", istate, ifrag, x+xa, y+xb, z+xc);
                            
                            // geometry section
                            fprintf(fs, "geometry nocenter noautoz units angstrom\n");
                            fprintf(fs, "symmetry c1\n");
        
                            int ra = 2*xa+1;
                            int rb = 2*xb+1;
                            int rc = 2*xc+1;
                            int totalatoms = ra*rb*rc*natoms;
                            int state_start = istate*totalatoms;
                            //int totalatoms = istate*ra*rb*rc*natoms;
                            
                            for (int iatom=state_start; iatom<state_start+totalatoms; ++iatom) {
                                if (atom->AtomInFragment(iatom,jfrag,istate,x,y,z)) {
                                    fprintf(fs, "  %c %4.2lf %20.10lf %20.10lf %20.10lf\n",
                                            atom->symbol[iatom%natoms],
                                            atom->getAtomicNumber(iatom%natoms),
                                            atom->coord[3*(iatom%natoms)]   + x*cellA,
                                            atom->coord[3*(iatom%natoms)+1] + y*cellB,
                                            atom->coord[3*(iatom%natoms)+2] + z*cellC
                                            );
                                    
                                }
                                else if (atom->AtomInFragment(iatom,ifrag,istate,0,0,0)) {
                                    fprintf(fs, "  %c %4.2lf %20.10lf %20.10lf %20.10lf\n",
                                            atom->symbol[iatom%natoms],
                                            atom->getAtomicNumber(iatom%natoms),
                                            atom->coord[3*(iatom%natoms)],
                                            atom->coord[3*(iatom%natoms)+1],
                                            atom->coord[3*(iatom%natoms)+2]
                                            );
                                }
                            }
                            fprintf(fs, "end\n\n");
                            
                            // $GUESS section
                            // Read previous step MO coeffs?
                            //if (flag_read_MOs) {
                            //    fprintf(fs, " $GUESS ");
                            //    fprintf(fs, "GUESS=MOREAD ");
                            //    fprintf(fs, "$END\n");
                            //}
                            
                            // bq section
                            fprintf(fs, "bq units angstrom\n");
                            for (int x0=-afield; x0<=afield; ++x0) {
                                for (int y0=-bfield; y0<=bfield; ++y0) {
                                    for (int z0=-cfield; z0<=cfield; ++z0) {
                                        
                                        for (int iatom=0; iatom<natoms; ++iatom) {
                                            if (atom->fragment[istate*natoms + iatom] != ifrag || x0!=0 || y0!=0 || z0!=0) {
                                                if (atom->fragment[istate*natoms + iatom] != jfrag || x0!=x || y0!=y || z0!=z) {
                                                    double mmq = atom->getCharge(iatom, istate);
                                                    fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                                                            atom->coord[3*iatom]   + x0*cellA,
                                                            atom->coord[3*iatom+1] + y0*cellB,
                                                            atom->coord[3*iatom+2] + z0*cellC,
                                                            mmq
                                                            );
                                                }
                                            }
                                        }
                                        
                                    }
                                }
                            }
                            fprintf(fs, "end\n\n");
                            
                            // scf section
                            fprintf(fs, "scf\n");
                            fprintf(fs, "singlet\n");
                            fprintf(fs, "direct\n");
                            fprintf(fs, "thresh 1e-6\n");
                            fprintf(fs, "sym off\n");
                            fprintf(fs, "end\n\n");
                            
                            // correlation section
                            fprintf(fs, "mp2\n");
                            fprintf(fs, "freeze core atomic\n");
                            fprintf(fs, "end\n\n");
                            
                            // charge section
                            if (ifrag == chgfrag && jfrag == chgfrag) {
                                fprintf(fs, "charge 2\n\n");
                            } else if (ifrag == chgfrag || jfrag == chgfrag) {
                                fprintf(fs, "charge 1\n\n");
                            } else {
                                fprintf(fs, "charge 0\n\n");
                            }
                        
			    // scratch section
			    char fname[256];
			    char scratch[256];
			    sprintf(fname, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d.nw", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
			    sprintf(scratch, "%s/%s/",run->scratch_dir,fname);
			    fprintf(fs, "scratch_dir %s\n\n",scratch);
    
                            // basis set section
                            fprintf(fs, "basis\n");
                            fprintf(fs, "* library %s\n", run->basis);
                            fprintf(fs, "end\n\n");
                            
                            // task section
                            if (jobtype == RUN_ENERGY)
                                fprintf(fs, "task %s energy\n\n", run->correlation);
                            else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                                fprintf(fs, "task %s gradient\n\n", run->correlation);
                            
                            /*
                             *
                             */
                            
                            
                            // Comment for labeling
                            fprintf(fs, "\n\nstart field\n");
                            fprintf(fs, "title 'State %d Monomer %d Cell %d %d %d'\n\n", istate, ifrag, x+xa, y+xb, z+xc);
                            
                            // geometry section
                            fprintf(fs, "geometry nocenter noautoz units angstrom\n");
                            fprintf(fs, "symmetry c1\n");
                            
                            int ra = 2*xa+1;
                            int rb = 2*xb+1;
                            int rc = 2*xc+1;
                            int totalatoms = ra*rb*rc*natoms;
                            int state_start = istate*totalatoms;
                            //int totalatoms = istate*ra*rb*rc*natoms;
                            
                            for (int iatom=state_start; iatom<state_start+totalatoms; ++iatom) {
                                if (atom->AtomInFragment(iatom,jfrag,istate,x,y,z)) {
                                    fprintf(fs, "  %c %4.2lf %20.10lf %20.10lf %20.10lf\n",
                                            atom->symbol[iatom%natoms],
                                            atom->getAtomicNumber(iatom%natoms),
                                            atom->coord[3*(iatom%natoms)]   + x*cellA,
                                            atom->coord[3*(iatom%natoms)+1] + y*cellB,
                                            atom->coord[3*(iatom%natoms)+2] + z*cellC
                                            );
                                    
                                }
                                else if (atom->AtomInFragment(iatom,ifrag,istate,0,0,0)) {
                                    fprintf(fs, "  %c %4.2lf %20.10lf %20.10lf %20.10lf\n",
                                            atom->symbol[iatom%natoms],
                                            atom->getAtomicNumber(iatom%natoms),
                                            atom->coord[3*(iatom%natoms)],
                                            atom->coord[3*(iatom%natoms)+1],
                                            atom->coord[3*(iatom%natoms)+2]
                                            );
                                }
                            }
                            fprintf(fs, "end\n\n");
                            
                            // $GUESS section
                            // Read previous step MO coeffs?
                            //if (flag_read_MOs) {
                            //    fprintf(fs, " $GUESS ");
                            //    fprintf(fs, "GUESS=MOREAD ");
                            //    fprintf(fs, "$END\n");
                            //}
                            
                            // bq section
                            fprintf(fs, "bq units angstrom\n");
			    fprintf(fs, "force fmo_st%s_d%03d-%03d_cell.%d.%d.%d.nw.field\n", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
                            for (int x0=-afield; x0<=afield; ++x0) {
                                for (int y0=-bfield; y0<=bfield; ++y0) {
                                    for (int z0=-cfield; z0<=cfield; ++z0) {
                                        
                                        for (int iatom=0; iatom<natoms; ++iatom) {
                                            if (atom->fragment[istate*natoms + iatom] != ifrag || x0!=0 || y0!=0 || z0!=0) {
                                                if (atom->fragment[istate*natoms + iatom] != jfrag || x0!=x || y0!=y || z0!=z) {
                                                    double mmq = atom->getCharge(iatom, istate);
                                                    fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                                                            atom->coord[3*iatom]   + x0*cellA,
                                                            atom->coord[3*iatom+1] + y0*cellB,
                                                            atom->coord[3*iatom+2] + z0*cellC,
                                                            mmq
                                                            );
                                                }
                                            }
                                        }
                                        
                                    }
                                }
                            }
                            fprintf(fs, "end\n\n");
                            
                            // scf section
                            fprintf(fs, "scf\n");
                            fprintf(fs, "singlet\n");
                            fprintf(fs, "direct\n");
                            fprintf(fs, "thresh 1e-6\n");
                            fprintf(fs, "sym off\n");
                            fprintf(fs, "end\n\n");
                            
                            // charge section
                            if (ifrag == chgfrag && jfrag == chgfrag) {
                                fprintf(fs, "charge 2\n\n");
                            } else if (ifrag == chgfrag || jfrag == chgfrag) {
                                fprintf(fs, "charge 1\n\n");
                            } else {
                                fprintf(fs, "charge 0\n\n");
                            }
                           
 			    // scratch section
 			    fprintf(fs, "scratch_dir %s\n\n",scratch);
 
                            // basis set section
                            fprintf(fs, "basis\n");
                            fprintf(fs, "* library %s\n", run->basis);
                            fprintf(fs, "end\n\n");
                            
                            // task section
                            if (jobtype == RUN_ENERGY)
                                fprintf(fs, "task scf energy\n\n", run->correlation);
                            else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                                fprintf(fs, "task scf gradient\n\n", run->correlation);

                            
                            fclose(fs);
                        } 
                    } // close loop over fragments for dimers
                    
                }
            }
        }
    } // close loop over states
      
      printf("Done writing NWChem inputs.\n");
  }
    
    // Hold up
    MPI_Barrier(fmr->world);
}


/*-----------------------------------------------------------------
  Perform all the FMO calculations 
-----------------------------------------------------------------*/
void Run::do_fmo_calculations(int FORCE)
{

}
