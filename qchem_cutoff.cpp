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
 Write NWChem inputs for each state's FMO calculations
 -----------------------------------------------------------------*/
void State::write_qchem_inputs_cutoff(int jobtype)
{
    // Writes a separate input file for all monomers and all dimers
    // Master rank does all the work here
    
    if (fmr->master_rank) { printf("Writing Q-Chem inputs.\n"); printf("Read MOs: %d\n", flag_read_MOs); }
    
    Run *run       = fmr->run;
    Atom *atom     = fmr->atom;
    int natoms     = atom->natoms;
    int nstates    = atom->nstates;
    int nfragments = atom->nfragments;
    int my_rank    = fmr->my_rank;
    int nf2        = nfragments*nfragments;
    
    int cellA      = fmr->atom->cellA;
    int cellB      = fmr->atom->cellB;
    int cellC      = fmr->atom->cellC;
    
    int na	 = 2*fmr->atom->na + 1;
    int nb 	 = 2*fmr->atom->nb + 1;
    int nc 	 = 2*fmr->atom->nc + 1;
    
    int xa         = fmr->atom->na;
    int xb         = fmr->atom->nb;
    int xc         = fmr->atom->nc;
    
    int afield     = fmr->atom->afield;
    int bfield     = fmr->atom->bfield;
    int cfield     = fmr->atom->cfield;
    
    bool python    = false;
    
    double comi[3],comj[3],comk[3];
    double massi,massj,massk,tmp,d2,d;
    
    int statedimers=0;
    int statemonomers=nfragments*nstates;
    int idx = 0;
    int pos,idxj,idxi;
    int skip;
 
    /*This is in the new branch
     Distributed
     */
    
    //run->n_monomers = nstates * nfragments * na*nb*nc;
    // Assuming all states have equal number of dimers and monomers, for now
    int nmonomers = nfragments * na*nb*nc;
    int ndimers = (nf2 * (na*nb*nc-1) + (nfragments * (nfragments-1)) / 2);

    run->n_monomers_tmp = nstates * nmonomers;
    //run->n_dimers_tmp = nstates * ndimers;
    run->n_dimers_sq = nstates * nf2 *na*nb*nc; // inclues self

    // initialize monomer list
    run->monomer_list = new int [run->n_monomers_tmp];
    for (int i=0; i<run->n_monomers_tmp; ++i) run->monomer_list[i] = -1;

    // include zeroth cell monomers in list
    //if (fmr->my_rank == 0) {
    for (int istate=0; istate<nstates; ++istate) {
        for (int ifrag=0; ifrag<nfragments; ++ifrag) {
            pos=atom->getMonomerIndex(istate,0,0,0,ifrag);
            //printf("POS:%d\n",pos);
            //run->monomer_list[istate*nmonomers + ifrag] = pos;
            run->monomer_list[istate*nfragments + ifrag] = pos;
        }
        //exit(0);
    }
    //}
    //exit(0);

    // initialize dimer list
    run->dimer_list = new int [run->n_dimers_sq];
    for (int i=0; i<run->n_dimers_sq; ++i) run->dimer_list[i] = -1;
    
    // initialize monomer queue list
    run->monomer_queue = new int [run->n_monomers_tmp];
    for (int i=0; i<run->n_monomers_tmp; ++i) run->monomer_queue[i] = 0;
    
    // include zeroth cell monomers in queue list
    for (int istate=0; istate<nstates; ++istate) {
        pos=atom->getMonomerIndex(istate,0,0,0,0);
        for (int ifrag=0; ifrag<nfragments; ++ifrag) run->monomer_queue[pos+ifrag] = 1;
    }
    
    // initialize dimer queue list
    run->dimer_queue = new int [run->n_dimers_sq];
    for (int i=0; i<run->n_dimers_sq; ++i) run->dimer_queue[i] = 0;


    // search for nearest dimers
    for (int x=-xa; x<=xa; ++x) {
        for (int y=-xb; y<=xb; ++y) {
            for (int z=-xc; z<=xc; ++z) {
                
                for (int ifrag=0; ifrag<nfragments; ++ifrag) {
                    
                    // center of mass of i fragment
                    comi[0] = comi[1] = comi[2] = massi = 0.0;
                    for (int i=0; i<natoms; ++i) {
                        if (atom->fragment[i] == ifrag) {
                            comi[0] += atom->mass[i] * atom->coord[3*i];
                            comi[1] += atom->mass[i] * atom->coord[3*i+1];
                            comi[2] += atom->mass[i] * atom->coord[3*i+2];
                            massi += atom->mass[i];
                        }
                    }
                    comi[0] /= massi;
                    comi[1] /= massi;
                    comi[2] /= massi;
                    
                    //printf("for dimers\n");
                    for (int jfrag=0; jfrag<nfragments; ++jfrag) {
                        
                        if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;
                        
                        // center of mass of j fragment
                        comj[0] = comj[1] = comj[2] = massj = 0.0;
                        for (int j=0; j<natoms; ++j) {
                            if (atom->fragment[j] == jfrag) {
                                comj[0] += atom->mass[j] * (atom->coord[3*j]   + x*cellA);
                                comj[1] += atom->mass[j] * (atom->coord[3*j+1] + y*cellB);
                                comj[2] += atom->mass[j] * (atom->coord[3*j+2] + z*cellC);
                                massj += atom->mass[j];
                            }
                        }
                        comj[0] /= massj;
                        comj[1] /= massj;
                        comj[2] /= massj;
                        
                        // distance between comi and comj
                        d2 = 0;
                        for (int l=0; l<3; ++l) {
                            tmp = comj[l]-comi[l];
                            d2 += tmp*tmp;
                        }
                        
                        d=sqrt(d2);
                        
                        if (d <= run->cut_dimer) {
                            
                            
                            //statedimers++;
                            //idx =atom->getDimerIndex(0,x,y,z,ifrag,jfrag);
                            for (int state=0; state<nstates; ++state) {
                                idx=atom->getDimerIndex(state,x,y,z,ifrag,jfrag);
                                run->dimer_list[statedimers] = idx;
                                run->dimer_queue[idx] = 1;
                                statedimers++;
                                
                            }
                            
                            // include j monomer in queue list
                            if (x != 0 || y != 0 || z != 0) {

				//printf("found....\n");
                                for (int state=0; state<nstates; ++state) {
                                    
                                    skip=0;
                                    idxj=atom->getMonomerIndex(state,x,y,z,jfrag);
                                    for (int l=0; l<statemonomers; ++l) {
                                        if (run->monomer_list[l] == idxj) {
                                            skip=1;
                                            break;
                                        } 
				    }
                                    if (skip == 0) {
                                      run->monomer_list[statemonomers] = idxj;
                                      run->monomer_queue[idxj] = 1;
                                      statemonomers++;
 				      //printf("added\n");
                                    } else { //printf("HELLO:%d\n",idxj); }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    MPI_Barrier(fmr->world);
    //if (fmr->my_rank == 0) {
    //printf("START\n");
    //for (int l=0; l<run->n_monomers_tmp; ++l) printf("%5d %d\n",l,run->monomer_list[l]);
    //printf("number of monomers:%d\n",statemonomers);
    //printf("END\n");
    //exit(0);
    //}
    
    // Assuming all states have equal number of dimers and monomers, for now
    run->n_monomers = statemonomers;//nstates * statemonomers;
    run->n_dimers   = statedimers;//nstates * statedimers;

    // ** Determine load balance ** //
    int div = run->n_monomers / fmr->world_size;
    int rem = run->n_monomers % fmr->world_size;
    //printf("n_monomers:%d\nworld_Size:%d\n",run->n_monomers,fmr->world_size);
    //printf("DIV:%d\nREM:%d\n",div,rem);

    //delete [] run->monomer_proc;
    //delete [] run->dimer_proc;
    run->monomer_proc = run->dimer_proc = NULL;
    int ifrom_mono, ito_mono;
    if (my_rank < rem) {
        ifrom_mono = my_rank*div + my_rank;
        ito_mono   = ifrom_mono + div + 1;
    } else {
        ifrom_mono = my_rank*div + rem;
        ito_mono   = ifrom_mono + div;
    }
    //printf("monomer initialize\n");
    //** Distributed Processors **//
    if (run->monomer_proc == NULL) { run->monomer_proc = new int [div+1]; }
    for (int i=0; i<div+1; ++i) run->monomer_proc[i] = -1;
    //printf("monomer_div\n");
    div = run->n_dimers / fmr->world_size;
    rem = run->n_dimers % fmr->world_size;
    int ifrom_dim, ito_dim;
    if (my_rank < rem) {
        ifrom_dim = my_rank*div + my_rank;
        ito_dim   = ifrom_dim + div + 1;
    } else {
        ifrom_dim = my_rank*div + rem;
        ito_dim   = ifrom_dim + div;
    }

    //** Distributed Processors **//
    if (run->dimer_proc == NULL) run->dimer_proc = new int [div+1];
    //for (int i=0; i<div+1; ++i) run->dimer_proc[i] = -1;
    //printf("dimer_div\n");
    int index_mono = 0;
    int index_dim  = 0;
   
    //printf("ifrom:%4d index_mono:%4d ito:%4d myrank:%d\n",ifrom_mono,index_mono,ito_mono,fmr->my_rank); 
    // ***** Loop through monomers ***** //
    for (int imon=0; imon<run->n_monomers_tmp; imon++) {
        //printf("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH...rank:%d\n",fmr->my_rank);
        //printf("imon:%d monomer:%d\n",imon,run->monomer_list[imon]);
        if (run->monomer_list[imon] == -1) continue;
        //if (run->monomer_queue[imon] == 0) continue;
        //printf("ifrom:%4d index_mono:%4d ito:%4d myrank:%d\n",ifrom_mono,index_mono,ito_mono,fmr->my_rank);
        if (ifrom_mono <= index_mono && index_mono < ito_mono) {

            run->monomer_proc[index_mono-ifrom_mono]=run->monomer_list[imon];

            int istate;
            int x,y,z;
            int ifrag;
            
            atom->getMonomerIndices(run->monomer_list[imon],istate,x,y,z,ifrag);
            //atom->getMonomerIndices(imon,istate,x,y,z,ifrag);
            //printf("testing:%5d rank:%2d\n", imon, fmr->my_rank);
            //printf("state:%d ifrag:%d x:%d y:%d z:%d rank:%d\n", istate, ifrag, x, y, z, my_rank);
            
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
            int midx = 0;
            int didx = 0;
            
            sprintf(snum, "%02d", istate);
            sprintf(state_directory, "state_%02d", istate);
            // Make the directory...
            //            
            sprintf(make_directory, "mkdir -p %s/%s", run->scratch_dir,state_directory);
            //printf("%s\n",make_directory);
            int ierr = system(make_directory);
            
            // Get name of file to open
            char jobname[256];
            char filename[256];

            // Get name of job
            sprintf(jobname, "fmo_st%s_m%03d_cell.%d.%d.%d", snum, ifrag, x+xa, y+xb, z+xc);
            sprintf(filename, "/%s/%s/fmo_st%s_m%03d_cell.%d.%d.%d.in", run->scratch_dir,state_directory, snum, ifrag, x+xa, y+xb, z+xc);
            //printf("%s\n",filename);

            FILE *fs = fopen(filename, "w");
            if (fs == NULL) {
                char tmpstr[256];
                sprintf(tmpstr, "Failure to write Q-Chem input for file %s", filename);
                fmr->error(FLERR, tmpstr);
            }
            
            // Comment for labeling
            fprintf(fs, "$comment\n");
            fprintf(fs, "State %d Monomer %d Cell %d %d %d\n", istate, ifrag, x+xa, y+xb, z+xc);
            fprintf(fs, "$end\n\n");
            
            // $rem section
            fprintf(fs, "$rem\n");
            if (jobtype == RUN_ENERGY)
                fprintf(fs, "jobtype sp\n");
            else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                fprintf(fs, "jobtype force\n");
            //fprintf(fs, "exchange pbe0\n");
            //fprintf(fs, "basis 6-31+G*\n");
            //fprintf(fs, "basis cc-pvdz\n");
            fprintf(fs, "basis %s\n", run->basis);
            //fprintf(fs, "aux_basis rimp2-cc-pvdz\n");
            //fprintf(fs, "exchange hf\n");
            fprintf(fs, "exchange %s\n", run->exchange);
            //fprintf(fs, "correlation mp2\n");
            fprintf(fs, "correlation %s\n", run->correlation);
            //fprintf(fs, "correlation rimp2\n");
            fprintf(fs, "scf_algorithm %s\n", run->algorithm);
            fprintf(fs, "scf_convergence 6\n");
            fprintf(fs, "max_scf_cycles 500\n");
            fprintf(fs, "qm_mm true\n");
            fprintf(fs, "print_input true\n");
            fprintf(fs, "sym_ignore true\n");
            fprintf(fs, "no_reorient true\n");
            fprintf(fs, "skip_charge_self_interact 1\n");
	    fprintf(fs, "mem_static 1000\n");
            //fprintf(fs, "gaussian_blur true\n");
            // Read previous step MO coeffs?
            //if (flag_read_MOs) fprintf(fs, "scf_guess read\n");
            fprintf(fs, "$end\n\n");
            
            // $molecule section
            fprintf(fs, "$molecule\n");
            if (ifrag == chgfrag) {
                fprintf(fs, "1 1\n");
            } else {
                fprintf(fs, "0 1\n");
            }
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (atom->fragment[istate*natoms + iatom] == ifrag) {
                    fprintf(fs, "%c %20.10lf %20.10lf %20.10lf\n",
                            atom->symbol[iatom],
                            atom->coord[3*iatom] + x*cellA,
                            atom->coord[3*iatom+1] + y*cellB,
                            atom->coord[3*iatom+2] + z*cellC
                            );
                }
            }
            fprintf(fs, "$end\n\n");
            
            // $external_charges section
            fprintf(fs, "$external_charges\n");
            for (int x0=-afield; x0<=afield; ++x0) {
                for (int y0=-bfield; y0<=bfield; ++y0) {
                    for (int z0=-cfield; z0<=cfield; ++z0) {
                        
                        for (int iatom=0; iatom<natoms; ++iatom) {
                            if (atom->fragment[istate*natoms + iatom] != ifrag || x0 != x || y0 != y || z0 != z) {
                                double mmq = atom->getCharge(iatom, istate);
                                //fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf 0.12446\n",
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
            
            fprintf(fs, "$end\n\n");
            
            fclose(fs);

        }
        ++index_mono;
    }
    
    // ***** Loop through dimers ***** //
    for (int idim=0; idim<run->n_dimers_sq; idim++) {
        
        //if (run->dimer_queue[idim] == 0) continue;
        if (run->dimer_list[idim] == -1) continue;
        
        if (ifrom_dim <= index_dim && index_dim < ito_dim) {
            
            run->dimer_proc[index_dim-ifrom_dim]=run->dimer_list[idim];


            int istate;
            int x,y,z;
            int ifrag,jfrag;
            
            //atom->getDimerIndices(idim,istate,x,y,z,ifrag,jfrag);
	    atom->getDimerIndices(run->dimer_list[idim],istate,x,y,z,ifrag,jfrag);

            //
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
            int midx = 0;
            int didx = 0;
            
            sprintf(snum, "%02d", istate);
            sprintf(state_directory, "state_%02d", istate);
            // Make the directory...
            //
            sprintf(make_directory, "mkdir -p %s/%s", run->scratch_dir,state_directory);
            int ierr = system(make_directory);
            
            // Get name of file to open
            char jobname[256];
            char filename[256];
            char inum[16];
            char jnum[16];
            
            // Get name of job
            sprintf(jobname, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d", snum, ifrag, jfrag, x+xa, y+xb, z+xc);

            sprintf(filename, "%s/%s/fmo_st%s_d%03d-%03d_cell.%d.%d.%d.in", run->scratch_dir, state_directory, snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            FILE *fs = fopen(filename, "w");
            if (fs == NULL) {
                char tmpstr[256];
                sprintf(tmpstr, "Failure to write Q-Chem input for file %s", filename);
                fmr->error(FLERR, tmpstr);
            }
            
            // Comment for labeling
            fprintf(fs, "$comment\n");
            fprintf(fs, "State %d Dimer %d %d Cell %d %d %d\n", istate, ifrag, jfrag, x+xa, y+xb, z+xc);
            fprintf(fs, "$end\n\n");
            
            // $rem section
            fprintf(fs, "$rem\n");
            if (jobtype == RUN_ENERGY)
                fprintf(fs, "jobtype sp\n");
            else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                fprintf(fs, "jobtype force\n");
            //fprintf(fs, "exchange pbe0\n");
            //fprintf(fs, "basis 6-31+G*\n");
            fprintf(fs, "basis %s\n", run->basis);
            //fprintf(fs, "basis cc-pvdz\n");
            //fprintf(fs, "aux_basis rimp2-cc-pvdz\n");
            fprintf(fs, "exchange %s\n", run->exchange);
            //fprintf(fs, "exchange hf\n");
            fprintf(fs, "correlation %s\n", run->correlation);
            //fprintf(fs, "correlation mp2\n");
            
            //fprintf(fs, "correlation rimp2\n");
            fprintf(fs, "scf_algorithm %s\n", run->algorithm);
            fprintf(fs, "scf_convergence 6\n");
            fprintf(fs, "max_scf_cycles 500\n");
            fprintf(fs, "qm_mm true\n");
            fprintf(fs, "print_input true\n");
            fprintf(fs, "sym_ignore true\n");
            fprintf(fs, "no_reorient true\n");
            fprintf(fs, "skip_charge_self_interact 1\n");
            fprintf(fs, "mem_static 160\n");
            //fprintf(fs, "gaussian_blur true\n");
            // Read previous step MO coeffs?
            //if (flag_read_MOs) fprintf(fs, "scf_guess read\n");
            fprintf(fs, "$end\n\n");
            
            // $molecule section
            fprintf(fs, "$molecule\n");
            if (ifrag == chgfrag && jfrag == chgfrag) {
                fprintf(fs, "2 1\n");
            } else if (ifrag == chgfrag || jfrag == chgfrag) {
                fprintf(fs, "1 1\n");
            } else {
                fprintf(fs, "0 1\n");
            }
            
            //zeroth unit cell monomer
            //int posi = getFragPosition(istate,0,0,0,ifrag);
            //int posj = getFragPosition(istate,x,y,z,jfrag)
            
            int ra = 2*xa+1;
            int rb = 2*xb+1;
            int rc = 2*xc+1;
            
            int totalatoms = ra*rb*rc*natoms;
            int state_start = istate*totalatoms;
            //int totalatoms = istate*ra*rb*rc*natoms;
            
            for (int iatom=state_start; iatom<state_start+totalatoms; ++iatom) {
                if (atom->AtomInFragment(iatom,jfrag,istate,x,y,z)) {
                    fprintf(fs, "%c %20.10lf %20.10lf %20.10lf\n",
                            atom->symbol[iatom%natoms],
                            atom->coord[3*(iatom%natoms)]   + x*cellA,
                            atom->coord[3*(iatom%natoms)+1] + y*cellB,
                            atom->coord[3*(iatom%natoms)+2] + z*cellC
                            );
                    
                }
                else if (atom->AtomInFragment(iatom,ifrag,istate,0,0,0)) {
                    fprintf(fs, "%c %20.10lf %20.10lf %20.10lf\n",
                            atom->symbol[iatom%natoms],
                            atom->coord[3*(iatom%natoms)],
                            atom->coord[3*(iatom%natoms)+1],
                            atom->coord[3*(iatom%natoms)+2]
                            );
                }
            }
            
            fprintf(fs, "$end\n\n");
            
            // $external_charges section
            fprintf(fs, "$external_charges\n");
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
            fprintf(fs, "$end\n\n");
            
            fclose(fs);

        }
        ++index_dim;
    }

    if (fmr->master_rank) printf("Done writing Q-Chem inputs.\n");

    // Hold up
    MPI_Barrier(fmr->world);
}


/*-----------------------------------------------------------------
 Perform all the FMO calculations
 -----------------------------------------------------------------*/
void Run::do_qchem_calculations_cutoff(int FORCE)
{
    // Check if we should branch to the env approximation
    if (EnvApprox) {
        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
        //do_nwchem_calculations_env();
        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
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
  
    int cellA    = fmr->atom->cellA;
    int cellB    = fmr->atom->cellB;
    int cellC    = fmr->atom->cellC;

    int nafield	 = 2*fmr->atom->afield + 1;
    int nbfield	 = 2*fmr->atom->bfield + 1;
    int ncfield	 = 2*fmr->atom->cfield + 1;
    
    int afield     = fmr->atom->afield;
    int bfield     = fmr->atom->bfield;
    int cfield     = fmr->atom->cfield;

    double comi[3],comj[3],comk[3];
    double massi,massj,massk,tmp,d2,d;
    
    bool python    = false;
    
    //n_monomers = nstates * nfragments * na*nb*nc;
    // Assuming all states have equal number of dimers, for now
    //n_dimers = nstates * (nf2 * (na*nb*nc-1) + (nfragments * (nfragments-1)) / 2);
    //n_dimers_sq = nstates * nf2 *na*nb*nc; // inclues self
    
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
    //if (monomer_energies == NULL) monomer_energies = new double [n_monomers];
    //if (dimer_energies == NULL)   dimer_energies   = new double [n_dimers_sq];
    for (int i=0; i<nstates; ++i)     fmo_energies[i]     = 0.0;
    //for (int i=0; i<n_monomers; ++i)  monomer_energies[i] = 0.0;
    //for (int i=0; i<n_dimers_sq; ++i) dimer_energies[i]   = 0.0;
    
    if (FORCE) {
        // ** Allocate gradients and zero ** //
        if (fmo_gradients == NULL)     fmo_gradients     = new double[nstates*3*natoms];
        //if (monomer_gradients == NULL) monomer_gradients = new double[n_monomers*3*natoms];
        //if (dimer_gradients == NULL)   dimer_gradients   = new double[n_dimers_sq*3*natoms];
        for (int i=0; i<nstates*3*natoms; ++i)     fmo_gradients[i]     = 0.0;
        //for (int i=0; i<n_monomers*3*natoms; ++i)  monomer_gradients[i] = 0.0;
        //for (int i=0; i<n_dimers_sq*3*natoms; ++i) dimer_gradients[i]   = 0.0;
    }
    
    // ** Determine load balance ** //
    int mdiv = n_monomers / fmr->world_size;
    int rem = n_monomers % fmr->world_size;
    int ifrom_mono, ito_mono;
    if (my_rank < rem) {
        ifrom_mono = my_rank*mdiv + my_rank;
        ito_mono   = ifrom_mono + mdiv + 1;
    } else {
        ifrom_mono = my_rank*mdiv + rem;
        ito_mono   = ifrom_mono + mdiv;
    }

    //
    //Only allocate monomer energies and gradients as needed
    if (monomer_energies == NULL) monomer_energies = new double [mdiv+1];
    for (int i=0; i<mdiv+1; ++i)  monomer_energies[i] = 0.0;
    if (FORCE) {
        if (monomer_gradients == NULL) monomer_gradients = new double[(mdiv+1)*3*natoms];
        for (int i=0; i<(mdiv+1)*3*natoms; ++i)  monomer_gradients[i] = 0.0;
    }

    int ddiv = n_dimers / fmr->world_size;
    rem = n_dimers % fmr->world_size;
    int ifrom_dim, ito_dim;
    if (my_rank < rem) {
        ifrom_dim = my_rank*ddiv + my_rank;
        ito_dim   = ifrom_dim + ddiv + 1;
    } else {
        ifrom_dim = my_rank*ddiv + rem;
        ito_dim   = ifrom_dim + ddiv;
    }

    //
    //Only allocate dimer energies and gradients as needed
    if (dimer_energies == NULL)   dimer_energies   = new double [ddiv+1];
    for (int i=0; i<ddiv+1; ++i) dimer_energies[i]   = 0.0;
    if (FORCE) {
        if (dimer_gradients == NULL) dimer_gradients = new double[(ddiv+1)*3*natoms];
        for (int i=0; i<(ddiv+1)*3*natoms; ++i)  dimer_gradients[i] = 0.0;
    }
    //
    
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
    
    // ***** Loop through monomers ***** //
    for (int imon=0; imon<n_monomers_tmp; imon++) {
        
        if (run->monomer_list[imon] == -1) continue;
	//if (run->monomer_queue[imon] == 0) continue;
        
        if (ifrom_mono <= index_mono && index_mono < ito_mono) {

            //** new index variable **//
            int index = index_mono - ifrom_mono;

            int istate;
            int x,y,z;
            int ifrag,jfrag;
            
            atom->getMonomerIndices(run->monomer_list[imon],istate,x,y,z,ifrag);
            //atom->getMonomerIndices(imon,istate,x,y,z,ifrag);

            char state_directory[256];
            char snum[16];
            char jobname[256];
            char filename[256];
            char inum[16];
            
            sprintf(snum, "%02d", istate);
            sprintf(state_directory, "state_%02d", istate);
            
            sprintf(inum,"%03d",ifrag);
            char cname[16];
            sprintf(cname,"cell.%d.%d.%d", x+xa, y+xb, z+xc);
            
            sprintf(jobname, "fmo_st%s_m%03d_%s", snum, ifrag, cname);
            sprintf(filename, "fmo_st%s_m%03d_%s",  snum, ifrag, cname);
            //printf("%s\n",filename);
            
            // change directory
            char directory[512];
            sprintf(directory, "%s/%s/", scratch_dir,state_directory);
            //chdir(directory);
    	
            sprintf(filename, "fmo_st%s_m%s_%s", snum, inum, cname);
            sprintf(command, "%s %s/%s/%s.in %s/%s/ > %s/%s/%s.out",
                    exec,
                    scratch_dir,
                    state_directory,
                    filename,
                    scratch_dir,
                    filename,
                    scratch_dir,
                    state_directory,
                    filename
                   );
            
            // ** The system call ** //
            ierr = system(command);
            //printf("error:%d\n",ierr);
            // ** Check for error ** //
            if (ierr) {
           //     if (ierr!=256 && ierr!=34304) {
                    printf("error:%d\n",ierr);
                    printf("Q-Chem run error on rank %d:\n", fmr->my_rank);
                    fmr->error(FLERR, command);
	    //	}
            }
            
            // ** Open output file and get the energy ** //
            char output_file[MAX_LENGTH];
            sprintf(output_file, "%s/%s/%s.in.energy", scratch_dir,state_directory, filename);
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
                    monomer_energies[index] = en;
                    //monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag] = en;
                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                }
            }
            fclose(fs);
            //printf("Rank %d: energy\n", my_rank);
            
            if (FORCE) {
                // ** Get gradient from file ** //
                sprintf(output_file, "%s/%s/%s.in.gradient", scratch_dir,state_directory, filename);
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
                        monomer_gradients[index*3*natoms + 3*(atnum%natoms)+0] = gx;
                        monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                        monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                        //monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                        //monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                        //monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                    }
                    // Increment atnum for the next round
                    atnum++;
                }
                fclose(fs);
                //printf("Rank %d: grad\n", my_rank);
                
                
                // ** Get field from file ** //
                sprintf(output_file, "%s/%s/%s.in.field", scratch_dir,state_directory, filename);
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
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)+0] = gx;
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                            //monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                            //monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                            //monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
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
    
    // ***** Loop through dimers ***** //
    for (int idim=0; idim<n_dimers_sq; idim++) {
       
        //if (run->dimer_queue[idim] == 0) continue;
        if (run->dimer_list[idim] == -1) continue;
        
        if (ifrom_dim <= index_dim && index_dim < ito_dim) {
        
            //** new index variable **//
            int index = index_dim - ifrom_dim;
            
            int istate;
            int x,y,z;
            int ifrag,jfrag;
            
            //
            //atom->getDimerIndices(idim,istate,x,y,z,ifrag,jfrag);
            atom->getDimerIndices(run->dimer_list[idim],istate,x,y,z,ifrag,jfrag);
            
            
            char state_directory[256];
            char snum[16];
            char jobname[256];
            char filename[256];
            
            sprintf(snum, "%02d", istate);
            sprintf(state_directory, "state_%02d", istate);
            sprintf(jobname, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            //printf("jobname %s %4d %2d\n",jobname,didx,my_rank);
            sprintf(filename, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            
            // change directory
            char directory[512];
            sprintf(directory, "%s/%s/", scratch_dir,state_directory);
            //chdir(directory);
            
            //sprintf(filename, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            sprintf(command, "%s %s/%s/%s.in %s/%s/ > %s/%s/%s.out",
                    exec,
                    scratch_dir,
                    state_directory,
                    filename,
                    scratch_dir,
                    filename,
                    scratch_dir,
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
            sprintf(output_file, "%s/%s/%s.in.energy", scratch_dir,state_directory, filename);
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
			dimer_energies[index] = en;
                        dimer_energies[index] = en;
                        //dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag] = en;
                        //dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag] = en;
                    } else {
                        dimer_energies[index] = en;
                        //dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag] = en;
                    }
                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                    
                }
	        }
	        fclose(fs);
            
            if (FORCE) {
                // ** Get gradient from file ** //
                sprintf(output_file, "%s/%s/%s.in.gradient", scratch_dir, state_directory, filename);
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

                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+0] = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+0] = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                            } else {
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+0] = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                            }
                        }
                    }
                    // Increment atnum for the next round
                    atnum++;
                }
                fclose(fs);
                
                // ** Get field from file ** //
                sprintf(output_file, "%s/%s/%s.in.field", scratch_dir, state_directory, filename);
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
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*jfrag + ifrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
                            } else {
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                //dimer_gradient_[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)]   = gx;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+1] = gy;
                                //dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*(atnum%natoms)+2] = gz;
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
    
    
    // Clock
    double FMO_end = MPI_Wtime();
    MPI_Barrier(fmr->world);
    if (fmr->master_rank) {
        printf("Finished with FMO calculations.\n");
        printf("Time for FMO: %.4f seconds\n\n", FMO_end - FMO_start);
    }

    // Clock
    MPI_Barrier(fmr->world);
    double comp_start = MPI_Wtime();
    
    chdir(fmr->home_dir);
    // *** Compute the FMO energies/forces for each state *** //
   
    if (fmr->master_rank) {
        printf("Compiling energies and gradients...\n");
    }
    int didx,midx;
    int idx,idxi,idxj;

    //** Distributed Compilation **//
    double den, men1, men2;

    FILE *fs = fopen("fmr_calc.log", "a");
    for (int istate=0; istate<nstates; ++istate) {
        //printf("----- State %4d -----\n", istate);
        if (fmr->print_level > 1) {
            fprintf(fs,"Rank %3d - Monomer | EI\n", fmr->my_rank);
        }
        
        double en_fmo1 = 0.0;
        for (int x=-xa; x<=xa; x++) {
            for (int y=-xb; y<=xb; y++) {
                for (int z=-xc; z<=xc; z++) {
                    
                    if (x==0 && y==0 && z==0) {

                        for (int ifrag=0; ifrag<nfragments; ++ifrag) {

                            midx=atom->getMonomerIndex(istate,x,y,z,ifrag);

                            men1=0.0;

                            for (int i=0; i<mdiv+1; i++) {
                                //				printf("proc:%d\n",monomer_proc[i]);
                                if (monomer_proc[i] == midx) { men1 = monomer_energies[i]; break; }
                            }


                            if (fmr->print_level > 1) {
                                fprintf(fs,"Rank %3d - FRAG I:%4d | %16.10f\n", fmr->my_rank,
                                        ifrag, men1);
                                        //ifrag, monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag]);
                            }
                            
                            en_fmo1 += men1;
//                            en_fmo1 += monomer_energies[nfragments*istate + ifrag];


                        }

                    }
                    
                }
            }
        }
        //printf("state:%d mon_en:%lf\n",istate,en_fmo1);

        if (fmr->print_level > 1) {
            fprintf(fs,"Rank %3d - Dimer | EIJ EI EJ (EIJ - EI - EJ)\n",fmr->my_rank);
        }
        double en_fmo2 = 0.0;
        
        for (int x=-xa; x<=xa; x++) {
            for (int y=-xb; y<=xb; y++) {
                for (int z=-xc; z<=xc; z++) {
                    
                    for (int ifrag=0; ifrag<nfragments; ++ifrag) {
                        for (int jfrag=0; jfrag<nfragments; ++jfrag) {
                            
                            if (x==0 && y==0 && z==0 && jfrag<=ifrag) continue;
                            
                            idx=atom->getDimerIndex(istate,x,y,z,ifrag,jfrag);
                            idxi=atom->getMonomerIndex(istate,0,0,0,ifrag);
                            idxj=atom->getMonomerIndex(istate,x,y,z,jfrag);

                            if (dimer_queue[idx] == 1) {

                                //printf("rank:%d  didx:%4d  idxi:%4d  idxj:%4d\n",fmr->my_rank,idx,idxi,idxj);
                                
                                //overcount for unit cell
                                double oc=0.5; if (x==0 && y==0 && z==0) oc=1.0;

                                //zero energies
                                den = men1 = men2 = 0.0;

                                for (int i=0; i<ddiv+1; i++) {
                                    if (dimer_proc[i]==idx) { 
                                        den = dimer_energies[i];
                                        break; }
                                }
                                for (int i=0; i<mdiv+1; i++) {
                                    if (monomer_proc[i] == idxi) men1 = monomer_energies[i];
                                    if (monomer_proc[i] == idxj) men2 = monomer_energies[i];
                                }


                                double etmp = den - men1 - men2;
                                //double etmp = dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag] -
                                //monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(xa) +   nc*nfragments*(xb) +   nfragments*(xc) +   ifrag] -
                                //monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag];
                                if (fmr->print_level > 1) {
                                    fprintf(fs,"Didx:%2d  idxi:%2d  idxj:%2d  State:%d Rank %3d - FRAG I:%02d--J:%02d  CELL x:%2d y:%2d z:%2d | %16.10f %16.10f %16.10f %16.10f\n", idx, idxi, idxj,istate, fmr->my_rank, ifrag, jfrag,x,y,z,
                                            den, men1, men2, etmp);


                                    //fprintf(fs,"Rank %3d - FRAG I:%02d--J:%02d  CELL x:%2d y:%2d z:%2d | %16.10f %16.10f %16.10f %16.10f\n", fmr->my_rank, ifrag, jfrag,x,y,z,
                                    //        dimer_energies[nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag],
                                   //         monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(xa) +   nc*nfragments*(xb) +   nfragments*(xc) +   ifrag],
                                   //         monomer_energies[nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag], etmp
                                   //         );
                                }
                                en_fmo2 += etmp*oc;
                            }
                        }
                    }
                    
                }
            }
        }
        fmo_energies[istate] = en_fmo1 + en_fmo2;

        // ** Do forces for this state ** //
        if (FORCE) {
            double gx, gy, gz;
            //gx = gy = gz = 0.0;
            // Monomer force
            for (int x=-xa; x<=xa; x++) {
                for (int y=-xb; y<=xb; y++) {
                    for (int z=-xc; z<=xc; z++) {
                        
                        if (x==0 && y==0 && z==0) {
                            for (int ifrag=0; ifrag<nfragments; ++ifrag) {

                                midx=atom->getMonomerIndex(istate,x,y,z,ifrag);



                                gx = gy = gz = 0.0;

                                for (int pi=0; pi<mdiv+1; pi++) {
                                    if (monomer_proc[pi] == midx) {
                                        for (int i=0; i<natoms; ++i) {
                                            gx = monomer_gradients[pi*3*natoms + 3*i];
                                            gy = monomer_gradients[pi*3*natoms + 3*i+1];
                                            gz = monomer_gradients[pi*3*natoms + 3*i+2];
                                            //printf("atomid %d: %f %f %f\n",i,gx,gy,gz);
                                            fmo_gradients[3*natoms*istate + 3*i]   += gx;
                                            fmo_gradients[3*natoms*istate + 3*i+1] += gy;
                                            fmo_gradients[3*natoms*istate + 3*i+2] += gz;
                                        }
                                        break;
                                    }
                                }
                                //for (int i=0; i<natoms; ++i) {
                                //    gx = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i];
                                //    gy = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i+1];
                                //    gz = monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + ifrag)*3*natoms + 3*i+2];
                                //    fmo_gradients[3*natoms*istate + 3*i]   += gx;
                                //    fmo_gradients[3*natoms*istate + 3*i+1] += gy;
                                //    fmo_gradients[3*natoms*istate + 3*i+2] += gz;
                                //}
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
                                
                                idx=atom->getDimerIndex(istate,x,y,z,ifrag,jfrag);
                                idxi=atom->getMonomerIndex(istate,0,0,0,ifrag);
                                idxj=atom->getMonomerIndex(istate,x,y,z,jfrag);
                                
                                if (dimer_queue[idx] == 1) {
                                    
                                    //overcount for unit cell
                                    double oc=1.0; if (x==0 && y==0 && z==0) oc=1.0;
                                    
                                    gx = gy = gz = 0.0;
                                    
                                    for (int pi=0; pi<ddiv+1; pi++) {
                                        if (dimer_proc[pi] == idx) {
                                            
                                            for (int i=0; i<natoms; ++i) {
                                                gx = dimer_gradients[pi*3*natoms + 3*i+0];
                                                gy = dimer_gradients[pi*3*natoms + 3*i+1];
                                                gz = dimer_gradients[pi*3*natoms + 3*i+2];
                                                //printf("%lf %lf %lf\n",gx,gy,gz);
                                                
                                                fmo_gradients[3*natoms*istate + 3*i+0] += gx*oc;
                                                fmo_gradients[3*natoms*istate + 3*i+1] += gy*oc;
                                                fmo_gradients[3*natoms*istate + 3*i+2] += gz*oc;
                                            }
                                            //printf("\n");
                                            break;
                                        }
                                    }
                                    
                                    for (int pi=0; pi<mdiv+1; pi++) {
                                        if (monomer_proc[pi] == idxi) {
                                            for (int i=0; i<natoms; ++i) {
                                                gx = monomer_gradients[pi*3*natoms + 3*i+0];
                                                gy = monomer_gradients[pi*3*natoms + 3*i+1];
                                                gz = monomer_gradients[pi*3*natoms + 3*i+2];
                                                
                                                fmo_gradients[3*natoms*istate + 3*i+0] -= gx*oc;
                                                fmo_gradients[3*natoms*istate + 3*i+1] -= gy*oc;
                                                fmo_gradients[3*natoms*istate + 3*i+2] -= gz*oc;
                                                //printf("%lf %lf %lf\n",monomer_gradients[pi*3*natoms + 3*i+0],monomer_gradients[pi*3*natoms + 3*i+1],monomer_gradients[pi*3*natoms + 3*i+2]);
                                            }
                                            //printf("\n");
                                        }
                                        if (monomer_proc[pi] == idxj) {
                                            for (int i=0; i<natoms; ++i) {
                                                gx = monomer_gradients[pi*3*natoms + 3*i+0];
                                                gy = monomer_gradients[pi*3*natoms + 3*i+1];
                                                gz = monomer_gradients[pi*3*natoms + 3*i+2];
                                                
                                                fmo_gradients[3*natoms*istate + 3*i+0] -= gx*oc;
                                                fmo_gradients[3*natoms*istate + 3*i+1] -= gy*oc;
                                                fmo_gradients[3*natoms*istate + 3*i+2] -= gz*oc;
                                                //printf("%lf %lf %lf\n",monomer_gradients[pi*3*natoms + 3*i+0],monomer_gradients[pi*3*natoms + 3*i+1],monomer_gradients[pi*3*natoms + 3*i+2]);
                                            }
                                            //printf("\n");
                                        }
                                    }
                                    
//                                    for (int i=0; i<natoms; ++i) {
//                                    
//                                        gx = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i];
//                                        gy = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+1];
//                                        gz = dimer_gradients[(nf2*na*nb*nc*istate + nb*nc*nf2*(x+xa) + nc*nf2*(y+xb) + nf2*(z+xc) + nfragments*ifrag + jfrag)*3*natoms + 3*i+2];
//                                        
//                                        gx -= monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i] +
//                                        monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i];
//                                        gy -= monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i+1] +
//                                        monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i+1];
//                                        gz -= monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(0+xa) + nc*nfragments*(0+xb) + nfragments*(0+xc) + ifrag)*3*natoms + 3*i+2] +
//                                        monomer_gradients[(nfragments*na*nb*nc*istate + nb*nc*nfragments*(x+xa) + nc*nfragments*(y+xb) + nfragments*(z+xc) + jfrag)*3*natoms + 3*i+2];
//                                        
//                                        fmo_gradients[3*natoms*istate + 3*i]   += gx*oc;
//                                        fmo_gradients[3*natoms*istate + 3*i+1] += gy*oc;
//                                        fmo_gradients[3*natoms*istate + 3*i+2] += gz*oc;
//                                        
//                                    }
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
    
    // *** Reduce the energies from the parallel Q-Chem calls *** //
    double *rbuffer;
    //if (FORCE) rbuffer = new double [MAX_SIZE];
    //else       rbuffer = new double [MAX_SIZE];
    if (FORCE) rbuffer = new double [nstates*natoms*3];
    else       rbuffer = new double [nstates];
    //printf("Initialized rbuffer arrays...\n");
    
    // FMO energies
    //printf("Reduced FMO energies...\n");
    for (int i=0; i<nstates; ++i) rbuffer[i] = 0.0;
    MPI_Allreduce(fmo_energies, rbuffer, nstates, MPI_DOUBLE, MPI_SUM, fmr->world);
    for (int i=0; i<nstates; ++i) fmo_energies[i] = rbuffer[i];
    
    // FMO gradients
    if (FORCE) {
        //printf("Reduced FMO gradients...\n");
        for (int i=0; i<nstates*natoms*3; ++i) rbuffer[i] = 0.0;
        MPI_Allreduce(fmo_gradients, rbuffer, nstates*natoms*3, MPI_DOUBLE, MPI_SUM, fmr->world);
        for (int i=0; i<nstates*natoms*3; ++i) fmo_gradients[i] = rbuffer[i];
    }
    
    delete [] rbuffer;
    
    // *** Compute the FMO energies/forces for each state *** //
    // Broadcast the FMO energy/force to worker ranks
    MPI_Bcast(fmo_energies, nstates, MPI_DOUBLE, MASTER_RANK, fmr->world);
    if (FORCE) {
        MPI_Bcast(fmo_gradients, nstates*3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
    }
    
#ifdef FMR_DEBUG
    /*
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
    */
#endif
    
    // Clock
    MPI_Barrier(fmr->world);
    double comp_end = MPI_Wtime();
    if (fmr->master_rank) {
        printf("Time for FMO compilation: %.4f seconds\n\n", comp_end - comp_start);
    }
}

