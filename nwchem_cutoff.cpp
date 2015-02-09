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
void State::write_nwchem_inputs_cutoff(int jobtype)
{
    // Writes a separate input file for all monomers and all dimers
    // Master rank does all the work here
    
    if (fmr->master_rank) { printf("Writing NWChem inputs.\n"); printf("Read MOs: %d\n", flag_read_MOs); }
    
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
    for (int istate=0; istate<nstates; ++istate) {
        for (int ifrag=0; ifrag<nfragments; ++ifrag) {
            pos=atom->getMonomerIndex(istate,0,0,0,ifrag);
            run->monomer_list[istate*nmonomers + ifrag] = pos;
        }
    }

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
                            
                            
                            //                            statedimers++;
                            //                            idx =atom->getDimerIndex(0,x,y,z,ifrag,jfrag);
                            for (int state=0; state<nstates; ++state) {
                                idx=atom->getDimerIndex(state,x,y,z,ifrag,jfrag);
                                run->dimer_list[statedimers] = idx;
                                run->dimer_queue[idx] = 1;
                                statedimers++;
                                
                            }
                            
                            // include j monomer in queue list
                            if (x != 0 || y != 0 || z != 0) {
                                for (int state=0; state<nstates; ++state) {
                                    
                                    idxj=atom->getMonomerIndex(state,x,y,z,jfrag);
                                    run->monomer_list[statemonomers] = idxj;
                                    run->monomer_queue[idxj] = 1;
                                    statemonomers++;
                                    
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Assuming all states have equal number of dimers and monomers, for now
    run->n_monomers = statemonomers;//nstates * statemonomers;
    run->n_dimers   = statedimers;//nstates * statedimers;

    // ** Determine load balance ** //
    int div = run->n_monomers / fmr->world_size;
    int rem = run->n_monomers % fmr->world_size;

    int ifrom_mono, ito_mono;
    if (my_rank < rem) {
        ifrom_mono = my_rank*div + my_rank;
        ito_mono   = ifrom_mono + div + 1;
    } else {
        ifrom_mono = my_rank*div + rem;
        ito_mono   = ifrom_mono + div;
    }

    //** Distributed Processors **//
    if (run->monomer_proc == NULL) run->monomer_proc = new int [div+1];
    for (int i=0; i<div+1; ++i) run->monomer_proc[i] = -1;

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
    for (int i=0; i<div+1; ++i) run->dimer_proc[i] = -1;

    int index_mono = 0;
    int index_dim  = 0;
    
    // ***** Loop through monomers ***** //
    for (int imon=0; imon<run->n_monomers_tmp; imon++) {

        if (run->monomer_list[imon] == -1) continue;
        
        if (ifrom_mono <= index_mono && index_mono < ito_mono) {

            run->monomer_proc[index_mono-ifrom_mono]=run->monomer_list[imon];

            int istate;
            int x,y,z;
            int ifrag;
            
            atom->getMonomerIndices(run->monomer_list[imon],istate,x,y,z,ifrag);
            printf("testing:%5d rank:%2d\n", imon, fmr->my_rank);
            printf("state:%d ifrag:%d x:%d y:%d z:%d\n", istate, ifrag, x, y, z);
            
            // Determine the charged reactive fragment for this state
            int chgfrag = 0;
            for (int i=0; i<natoms; ++i) {
                if (atom->reactive[istate*natoms + i]) {
                    chgfrag = atom->fragment[istate*natoms + i];
                    break;
                }
            }
            
	    int ionfrag = 0;
            for (int i=0; i<natoms; ++i) {
	        if (atom->symbol[i] == 'L') {
		    ionfrag = atom->fragment[istate*natoms + i];
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
            //            sprintf(make_directory, "mkdir -p %s", state_directory);
            sprintf(make_directory, "mkdir -p %s/%s", run->scratch_dir,state_directory);
            int ierr = system(make_directory);
            
            // Get name of file to open
            char jobname[256];
            char filename[256];

            // Get name of job
            sprintf(jobname, "fmo_st%s_m%03d_cell.%d.%d.%d", snum, ifrag, x+xa, y+xb, z+xc);
            sprintf(filename, "%s/%s/fmo_st%s_m%03d_cell.%d.%d.%d.nw", run->scratch_dir,state_directory, snum, ifrag, x+xa, y+xb, z+xc);

            FILE *fs = fopen(filename, "w");
            if (fs == NULL) {
                char tmpstr[256];
                sprintf(tmpstr, "Failure to write NWChem input for file %s", filename);
                fmr->error(FLERR, tmpstr);
            }
            
            // Comment for labeling
            fprintf(fs, "start grad_%s\n", jobname);
            fprintf(fs, "title \"State %d Monomer %d Cell %d %d %d\"\n\n", istate, ifrag, x+xa, y+xb, z+xc);
            
            // geometry section
            fprintf(fs, "geometry nocenter noautoz units angstrom\n");
            fprintf(fs, "symmetry c1\n");
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
            fprintf(fs, "end\n\n");
            
            //determine COM of fragment I
            comi[0] = comi[1] = comi[2] = massi = 0.0;
            for (int i=0; i<natoms; ++i) {
                if (atom->fragment[i] == ifrag) {
                    comi[0] += atom->mass[i] * atom->coord[3*i+0];
                    comi[1] += atom->mass[i] * atom->coord[3*i+1];
                    comi[2] += atom->mass[i] * atom->coord[3*i+2];
                    massi += atom->mass[i];
                }
            }
            comi[0] /= massi;
            comi[1] /= massi;
            comi[2] /= massi;
            
            // bq section
            fprintf(fs, "bq units angstrom\n");
            for (int x0=-afield; x0<=afield; ++x0) {
                for (int y0=-bfield; y0<=bfield; ++y0) {
                    for (int z0=-cfield; z0<=cfield; ++z0) {
                        
                        for (int iatom=0; iatom<natoms; ++iatom) {
                            if (atom->fragment[istate*natoms + iatom] != ifrag || x0 != x || y0 != y || z0 != z) {


			        //BUGBUGBUGBUGBUGUBUGUBUGBUGBUGUBUGUBUGUBUG
			        //if statement for distance criteria between ifrag and loop frag
			        //if distance is greater than cutoff skip charges
			        //else do normal execution
			        //BUGBUGBUGBUGBUGUBUGUBUGBUGBUGUBUGUBUGUBUG

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

/*
            fprintf(fs, "bq units angstrom\n");
            for (int x0=-afield; x0<=afield; ++x0) {
                for (int y0=-bfield; y0<=bfield; ++y0) {
                    for (int z0=-cfield; z0<=cfield; ++z0) {
                        
                        for (int jfrag=0; jfrag<nfragments; ++jfrag) {
                            
                            if (jfrag != ifrag || x0 != x || y0 != y || z0 != z) {
                                
                                //determine COM of fragment J
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
                                
                                
                                //COM distance from fragment I to J
                                d2 = 0;
                                for (int l=0; l<3; ++l) {
                                    tmp = comj[l]-comi[l];
                                    d2 += tmp*tmp;
                                }
                                d=sqrt(d2);
                                
                                if (d>25.0) continue;
                                
                                for (int jatom=0; jatom<natoms; ++jatom) {
                                    if (atom->fragment[jatom] == jfrag) {
                                        double mmq = atom->getCharge(jatom, istate);
                                        fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                                                atom->coord[3*jatom+0] + x0*cellA,
                                                atom->coord[3*jatom+1] + y0*cellB,
                                                atom->coord[3*jatom+2] + z0*cellC,
                                                mmq
                                                );
                                    }
                                }
                            }
                        }
                        
                    }
                }
            }
*/
            fprintf(fs, "end\n\n");
            
            // scf section
            fprintf(fs, "scf\n");
            fprintf(fs, "print low\n");
            fprintf(fs, "singlet\n");
            fprintf(fs, "direct\n");
            fprintf(fs, "thresh 1e-6\n");
            fprintf(fs, "sym off\n");
            //if (flag_read_MOs)
            //  fprintf(fs, "vectors input %s.movecs output %s.movecs\n",jobname,jobname);
            //else
            //  fprintf(fs, "vectors output %s.movecs\n",jobname);
            fprintf(fs, "end\n\n");
            
            if (strcmp(run->correlation,"mp2") == 0) {
                fprintf(fs, "mp2\n");
                fprintf(fs, "print low\n");
                fprintf(fs, "end\n");
            }
            
            // charge section
            if (ifrag == chgfrag) {
                fprintf(fs, "charge 1\n\n");
            } else if (ifrag == ionfrag) {
		fprintf(fs, "charge -1\n\n");
            } else {
                fprintf(fs, "charge 0\n\n");
            }
            
            
            // scratch section
            char scratch[256];
            sprintf(scratch, "%s/%s/", run->scratch_dir,state_directory);
            fprintf(fs, "scratch_dir %s\n", scratch);
            fprintf(fs, "permanent_dir %s\n\n",scratch);
            // Make the scratch directory...
            sprintf(make_directory, "mkdir -p %s", scratch);
            ierr = system(make_directory);
            
            
            // basis set section
            fprintf(fs, "basis\n");
            fprintf(fs, "* library %s\n", run->basis);
            fprintf(fs, "end\n\n");
            
            if (python) {
                
                fprintf(fs, "python\n");
                fprintf(fs, "  abc=task_gradient('mp2')\n");
                fprintf(fs, "  fener=open('%s.nw.energy','w')\n",jobname);
                fprintf(fs, "  fener.write('%%15.10f'%%(abc[0]))\n");
                fprintf(fs, "  fener.close()\n");
                fprintf(fs, "  fgrad=open('%s.nw.gradient','w')\n",jobname);
                fprintf(fs, "  for i in range(0,len(abc[1]),3):\n");
                fprintf(fs, "    fgrad.write('%%15.10f %%15.10f %%15.10f\\n'%%(abc[1][i+0],abc[1][i+1],abc[1][i+2]))\n");
                fprintf(fs, "  fgrad.close()\n");
                fprintf(fs, "end\n\n");
                
                fprintf(fs, "task python\n\n");
                
            } else {
                
                // task section
                if (jobtype == RUN_ENERGY)
                    fprintf(fs, "task %s %s\n\n", run->correlation, "energy");
                else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                    fprintf(fs, "task %s %s\n\n", run->correlation, "gradient");
                
            }
            /*
             *
             */
            
            // Comment for labeling
            fprintf(fs, "\n\nstart field_%s\n", jobname);
            fprintf(fs, "print none\n\n");
            
            // geometry section
            fprintf(fs, "geometry nocenter noautoz units angstrom\n");
            fprintf(fs, "symmetry c1\n");
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
            fprintf(fs, "end\n\n");
            
            // bq section
            fprintf(fs, "bq units angstrom\n");
            for (int x0=-afield; x0<=afield; ++x0) {
                for (int y0=-bfield; y0<=bfield; ++y0) {
                    for (int z0=-cfield; z0<=cfield; ++z0) {
                        
                        for (int iatom=0; iatom<natoms; ++iatom) {
                            if (atom->fragment[istate*natoms + iatom] != ifrag || x0 != x || y0 != y || z0 != z) {


                                //BUGBUGBUGBUGBUGUBUGUBUGBUGBUGUBUGUBUGUBUG
                                //if statement for distance criteria between ifrag and loop frag
                                //if distance is greater than cutoff skip charges
                                //else do normal execution
                                //BUGBUGBUGBUGBUGUBUGUBUGBUGBUGUBUGUBUGUBUG

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
/*
            fprintf(fs, "bq units angstrom\n");
            fprintf(fs, "force %s.nw.field\n", jobname);
            for (int x0=-afield; x0<=afield; ++x0) {
                for (int y0=-bfield; y0<=bfield; ++y0) {
                    for (int z0=-cfield; z0<=cfield; ++z0) {
                        
                        for (int jfrag=0; jfrag<nfragments; ++jfrag) {
                            
                            if (jfrag != ifrag || x0 != x || y0 != y || z0 != z) {
                                
                                //determine COM of fragment J
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
                                
                                
                                //COM distance from fragment I to J
                                d2 = 0;
                                for (int l=0; l<3; ++l) {
                                    tmp = comj[l]-comi[l];
                                    d2 += tmp*tmp;
                                }
                                d=sqrt(d2);
                                
                                if (d>25.0) continue;
                                
                                for (int jatom=0; jatom<natoms; ++jatom) {
                                    if (atom->fragment[jatom] == jfrag) {
                                        double mmq = atom->getCharge(jatom, istate);
                                        fprintf(fs, "%20.10lf %20.10lf %20.10lf %16.4lf\n",
                                                atom->coord[3*jatom+0] + x0*cellA,
                                                atom->coord[3*jatom+1] + y0*cellB,
                                                atom->coord[3*jatom+2] + z0*cellC,
                                                mmq
                                                );
                                    }
                                }
                            }
                        }
                        
                    }
                }
            }
*/
            fprintf(fs, "end\n\n");
            
            // scf section
            fprintf(fs, "scf\n");
            fprintf(fs, "print none\n");
            fprintf(fs, "singlet\n");
            fprintf(fs, "direct\n");
            fprintf(fs, "thresh 1e-6\n");
            fprintf(fs, "sym off\n");
            //fprintf(fs, "vectors input %s.movecs\n",jobname);
            fprintf(fs, "end\n\n");
            
            // charge section
            if (ifrag == chgfrag) {
                fprintf(fs, "charge 1\n\n");
            } else if (ifrag == ionfrag) {
                fprintf(fs, "charge -1\n\n");
            } else {
                fprintf(fs, "charge 0\n\n");
            }
            
            // scratch section
            fprintf(fs, "scratch_dir %s\n", scratch);
            fprintf(fs, "permanent_dir %s\n\n", scratch);

            //fprintf(fs, "memory 2 gb\n");
            //fprintf(fs, "memory hardfail\n");
            //fprintf(fs, "memory noverify\n\n");
            // basis set section
            fprintf(fs, "basis\n");
            fprintf(fs, "* library %s\n", run->basis);
            fprintf(fs, "end\n\n");
            
            // task section
            if (jobtype == RUN_ENERGY)
                fprintf(fs, "task scf energy\n\n");
            else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                fprintf(fs, "task scf gradient\n\n");
            
            
            fclose(fs);
            
            
        }
        ++index_mono;
    }
    
    // ***** Loop through dimers ***** //
    for (int idim=0; idim<run->n_dimers_sq; idim++) {
        
        if (run->dimer_queue[idim] == 0) continue;
        
        if (ifrom_dim <= index_dim && index_dim < ito_dim) {
            
            run->dimer_proc[index_dim-ifrom_dim]=idim;
            //printf("rank: %d idim: %3d\n",fmr->my_rank,idim);

            int istate;
            int x,y,z;
            int ifrag,jfrag;
            
            atom->getDimerIndices(idim,istate,x,y,z,ifrag,jfrag);
            //atom->getDimerIndices(run->dimer_list[idim],istate,x,y,z,ifrag,jfrag);
            
            //printf("testing:%2d rank:%d state:%d ifrag:%d jfrag;%d x:%d y:%d z:%d\n", idim, fmr->my_rank,istate, ifrag, jfrag, x, y, z);
            //
            // Determine the charged reactive fragment for this state
            int chgfrag = 0;
            for (int i=0; i<natoms; ++i) {
                if (atom->reactive[istate*natoms + i]) {
                    chgfrag = atom->fragment[istate*natoms + i];
                    break;
                }
            }
           
            int ionfrag = 0;
            for (int i=0; i<natoms; ++i) {
                if (atom->symbol[i] == 'L') {
                    ionfrag = atom->fragment[istate*natoms + i];
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
            //            sprintf(make_directory, "mkdir -p %s", state_directory);
            sprintf(make_directory, "mkdir -p %s/%s", run->scratch_dir,state_directory);
            int ierr = system(make_directory);
            
            // Get name of file to open
            char jobname[256];
            char filename[256];
            
            // Get name of job
            sprintf(jobname, "fmo_st%s_d%03d-%03d_cell.%d.%d.%d", snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            //printf(" %s %s %03d %03d %d %d %d\n",jobname,snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            sprintf(filename, "%s/%s/fmo_st%s_d%03d-%03d_cell.%d.%d.%d.nw", run->scratch_dir, state_directory, snum, ifrag, jfrag, x+xa, y+xb, z+xc);
            FILE *fs = fopen(filename, "w");
            if (fs == NULL) {
                char tmpstr[256];
                sprintf(tmpstr, "Failure to write NWChem input for file %s", filename);
                fmr->error(FLERR, tmpstr);
            }
            
            // Comment for labeling
            fprintf(fs, "start grad_%s\n", jobname);
            fprintf(fs, "title \"State %d Dimer %d %d Cell %d %d %d\"\n\n", istate, ifrag, jfrag, x+xa, y+xb, z+xc);
            
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
            fprintf(fs, "end\n\n");
            
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
            fprintf(fs, "print low\n");
            fprintf(fs, "singlet\n");
            fprintf(fs, "direct\n");
            fprintf(fs, "thresh 1e-6\n");
            fprintf(fs, "sym off\n");
            //if (flag_read_MOs)
            //  fprintf(fs, "vectors input %s.movecs output %s.movecs\n",jobname,jobname);
            //else
            //  fprintf(fs, "vectors output %s.movecs\n",jobname);
            fprintf(fs, "end\n\n");
            
            if (strcmp(run->correlation,"mp2") == 0) {
                fprintf(fs, "mp2\n");
                fprintf(fs, "print low\n");
                fprintf(fs, "end\n");
            }
            
            // charge section
            if (ifrag == chgfrag && jfrag == chgfrag) {
                fprintf(fs, "charge 2\n\n");
            } else if (ifrag == ionfrag && jfrag == ionfrag) {
                fprintf(fs, "charge -2\n\n");
            } else if (ifrag == chgfrag && jfrag == ionfrag) {
		fprintf(fs, "charge 0\n\n");
	    } else if (ifrag == ionfrag && jfrag == chgfrag) {
		fprintf(fs, "charge 0\n\n");
	    } else if (ifrag == ionfrag || jfrag == ionfrag) {
		fprintf(fs, "charge -1\n\n");
	    } else if (ifrag == chgfrag || jfrag == chgfrag) {
		fprintf(fs, "charge 1\n\n");
	    } else {	
                fprintf(fs, "charge 0\n\n");
            }
            
            // scratch section
            char scratch[256];
            sprintf(scratch, "%s/%s/", run->scratch_dir,state_directory);
            fprintf(fs, "scratch_dir %s\n", scratch);
            fprintf(fs, "permanent_dir %s\n\n", scratch);
            
            // Make the scratch directory...
            sprintf(make_directory, "mkdir -p %s", scratch);
            ierr = system(make_directory);
            
            // basis set section
            fprintf(fs, "basis\n");
            fprintf(fs, "* library %s\n", run->basis);
            fprintf(fs, "end\n\n");
            
            if (python) {
                
                fprintf(fs, "python\n");
                fprintf(fs, "  abc=task_gradient('mp2')\n");
                fprintf(fs, "  fener=open('%s.nw.energy','w')\n",jobname);
                fprintf(fs, "  fener.write('%%15.10f'%%(abc[0]))\n");
                fprintf(fs, "  fener.close()\n");
                fprintf(fs, "  fgrad=open('%s.nw.gradient','w')\n",jobname);
                fprintf(fs, "  for i in range(0,len(abc[1]),3):\n");
                fprintf(fs, "    fgrad.write('%%15.10f %%15.10f %%15.10f\\n'%%(abc[1][i+0],abc[1][i+1],abc[1][i+2]))\n");
                fprintf(fs, "  fgrad.close()\n");
                fprintf(fs, "end\n\n");
                fprintf(fs, "task python\n\n");
                
            } else {
                
                // task section
                if (jobtype == RUN_ENERGY)
                    fprintf(fs, "task %s energy\n\n", run->correlation);
                else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                    fprintf(fs, "task %s gradient\n\n", run->correlation);
                
            }
            /*
             *
             */
            
            
            // Comment for labeling
            fprintf(fs, "start field_%s\n", jobname);
            fprintf(fs, "title \"State %d Dimer %d %d Cell %d %d %d\"\n", istate, ifrag, jfrag, x+xa, y+xb, z+xc);
            fprintf(fs, "print none\n\n");
            
            // geometry section
            fprintf(fs, "geometry nocenter noautoz units angstrom\n");
            fprintf(fs, "symmetry c1\n");
            
            //int ra = 2*xa+1;
            //int rb = 2*xb+1;
            //int rc = 2*xc+1;
            //int totalatoms = ra*rb*rc*natoms;
            //int state_start = istate*totalatoms;
            //int totalatoms = istate*ra*rb*rc*natoms;
            
            for (int iatom=state_start; iatom<state_start+totalatoms; ++iatom) {
                if (atom->AtomInFragment(iatom,jfrag,istate,x,y,z)) {
                    fprintf(fs, "  %c %20.10lf %20.10lf %20.10lf\n",
                            atom->symbol[iatom%natoms],
                            atom->coord[3*(iatom%natoms)]   + x*cellA,
                            atom->coord[3*(iatom%natoms)+1] + y*cellB,
                            atom->coord[3*(iatom%natoms)+2] + z*cellC
                            );
                    
                }
                else if (atom->AtomInFragment(iatom,ifrag,istate,0,0,0)) {
                    fprintf(fs, "  %c %20.10lf %20.10lf %20.10lf\n",
                            atom->symbol[iatom%natoms],
                            atom->coord[3*(iatom%natoms)],
                            atom->coord[3*(iatom%natoms)+1],
                            atom->coord[3*(iatom%natoms)+2]
                            );
                }
            }
            fprintf(fs, "end\n\n");
            
            // bq section
            fprintf(fs, "bq units angstrom\n");
            fprintf(fs, "force %s.nw.field\n", jobname);
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
            fprintf(fs, "print none\n");
            //fprintf(fs, "vectors input %s.movecs\n",jobname);
            fprintf(fs, "end\n\n");
            
            // charge section
            // charge section
            if (ifrag == chgfrag && jfrag == chgfrag) {
                fprintf(fs, "charge 2\n\n");
            } else if (ifrag == chgfrag && jfrag == ionfrag) {
                fprintf(fs, "charge 0\n\n");
            } else if (ifrag == ionfrag && jfrag == chgfrag) {
                fprintf(fs, "charge 0\n\n");
            } else if (ifrag == ionfrag || jfrag == ionfrag) {
                fprintf(fs, "charge -1\n\n");
            } else if (ifrag == chgfrag || jfrag == chgfrag) {
                fprintf(fs, "charge 1\n\n");    
            } else {        
                fprintf(fs, "charge 0\n\n");
            }

            // scratch section
            fprintf(fs, "scratch_dir %s\n",scratch);
            fprintf(fs, "permanent_dir %s\n\n",scratch);

//            fprintf(fs, "memory 2 gb\n");
//	    fprintf(fs, "memory noverify\n\n");
            
            // basis set section
            fprintf(fs, "basis\n");
            fprintf(fs, "* library %s\n", run->basis);
            fprintf(fs, "end\n\n");
            
            
            // task section
            if (jobtype == RUN_ENERGY)
                fprintf(fs, "task scf energy\n\n");
            else if (jobtype == RUN_FORCE || jobtype == RUN_MOLDYN)
                fprintf(fs, "task scf gradient\n\n");
            
            
            fclose(fs);
            
        }
        ++index_dim;
    }
    
    if (fmr->master_rank) printf("Done writing NWChem inputs.\n");
    
    // Hold up
    MPI_Barrier(fmr->world);
}


/*-----------------------------------------------------------------
 Perform all the FMO calculations
 -----------------------------------------------------------------*/
void Run::do_nwchem_calculations_cutoff(int FORCE)
{
    // Check if we should branch to the env approximation
    if (EnvApprox) {
        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
        //do_nwchem_calculations_env();
        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
        return;
    }
    
    // Divide the list up across MPI ranks
    // Then, run the FMO calculations with call to serial NWChem in parallel
    
    
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
        
        if (ifrom_mono <= index_mono && index_mono < ito_mono) {

            //** new index variable **//
            int index = index_mono - ifrom_mono;

            int istate;
            int x,y,z;
            int ifrag,jfrag;
            
            atom->getMonomerIndices(run->monomer_list[imon],istate,x,y,z,ifrag);

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
            
            // change directory
            char directory[512];
            //            sprintf(directory, "%s/", state_directory);
            sprintf(directory, "%s/%s/", scratch_dir,state_directory);
            chdir(directory);
            
            //printf("%s %2d\n",jobname,fmr->my_rank);
            sprintf(command, "%s %s.nw > %s.nwout",
                    exec,
                    jobname,
                    jobname
                    );
            
            // ** The system call ** //
            ierr = system(command);
            
            // ** Check for error ** //
            //printf("ierr:%d\n",ierr);
            if (ierr) {
                printf("NWChem run error on rank %d:\n", fmr->my_rank);
                fmr->error(FLERR, command);
            }
            
            sprintf(command, "rm %s/%s/field_%s*",scratch_dir,state_directory,jobname);
            ierr = system(command);
            
            if (python) {
                // ** Open output file and get the energy ** //
                char output_file[MAX_LENGTH];
                sprintf(output_file, "%s.nw.energy",  jobname);
                FILE *fs = fopen(output_file, "r");
                if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                }
                char line[MAX_LENGTH];
                double en;
                
                
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    if ( sscanf(line, "%lf", &en) == 1 ) {
                        
                        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                        monomer_energies[index] = en;
                        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                        
                    }
                }
                fclose(fs);
                //printf("Rank %d: energy\n", my_rank);
                
                if (FORCE) {
                    // ** Get gradient from file ** //
                    sprintf(output_file, "%s.nw.gradient", jobname);
                    fs = fopen(output_file, "r");
                    if (fs == NULL) {
                        char tmpstr[MAX_LENGTH];
                        sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
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
                        if ( sscanf(line, "%lf %lf %lf", &gx, &gy, &gz) == 3 ) {
                            //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                            //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                        }
                        // Increment atnum for the next round
                        atnum++;
                    }
                    fclose(fs);
                    //printf("Rank %d: grad\n", my_rank);
                }
                
            } else {
                
                // ** Open output file and get the energy ** //
                char output_file[MAX_LENGTH];
                sprintf(output_file, "%s.nwout",  jobname);
                FILE *fs = fopen(output_file, "r");
                if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                }
                char line[MAX_LENGTH];
                double en;
                double gw, gx, gy, gz;
                
                char tmp0[16],tmp1[16],tmp2[16],tmp3[16],tmp4[16];
                int atnum;
                atnum = 0; // index of QM atom for storing gradient
                
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    
                    if (strcmp(run->correlation,"mp2") == 0) {
                        if ( strstr(line, "mp2 ENERGY GRADIENTS") ) {
                            
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            
                            for (int iatom=0; iatom<natoms; ++iatom) {
                                
                                fgets(line, MAX_LENGTH, fs);
                                if (strcmp(line,"\n") == 0) break;
                                while ( !fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z) ) {
                                    atnum++;
                                }
                                
                                if ( sscanf(line, "%lf %s %s %s %s %lf %lf %lf", &gw, tmp0, tmp1, tmp2, tmp3, &gx, &gy, &gz) == 8 ) {
                                    //printf("%s %d",line,atnum);
                                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                    monomer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                    monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                    monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                    
                                }
                                atnum++;
                            }
                            
                        } else if ( strstr(line, "Total MP2 energy") ) {
                            
                            if ( sscanf(line, "%s %s %s %lf", tmp0, tmp1, tmp2, &en) == 4 ) {
                                
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                monomer_energies[index] = en;
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                
                            }
                        }
                        
                    } else if (strcmp(run->correlation,"scf") == 0) {
                        
                        if ( strstr(line, "RHF ENERGY GRADIENTS") ) {
                            
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            
                            for (int iatom=0; iatom<natoms; ++iatom) {
                                
                                fgets(line, MAX_LENGTH, fs);
                                if (strcmp(line,"\n") == 0) break;
                                while ( !fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z) ) {
                                    atnum++;
                                }
                                
                                if ( sscanf(line, "%lf %s %s %s %s %lf %lf %lf", &gw, tmp0, tmp1, tmp2, tmp3, &gx, &gy, &gz) == 8 ) {
                                    //printf("%s %d",line,atnum);
                                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                    monomer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                    monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                    monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                    //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                    
                                }
                                atnum++;
                            }
                            
                        } else if ( strstr(line, "Total SCF energy") ) {
                            
                            if ( sscanf(line, "%s %s %s %s %lf", tmp0, tmp1, tmp2, tmp3, &en) == 5 ) {
                                
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                monomer_energies[index] = en;
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                
                            }
                        }
                    }
                }
                fclose(fs);
                
            }
            if (FORCE) {
                // ** Get field from file ** //
                char output_file[MAX_LENGTH];
                sprintf(output_file, "%s.nw.field", jobname);
                FILE *fs = fopen(output_file, "r");
                if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                }
                char line[MAX_LENGTH];
                double en;
                double gw, gx, gy, gz;
                char tmp0[16],tmp1[16],tmp2[16],tmp3[16],tmp4[16];
                int atnum;

                atnum = 0; // index of non-QM atom for storing gradient
                fgets(line, MAX_LENGTH, fs); // skip initial comment line
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    // Advance atnum until it matches as a non-QM atom index for this monomer fragment
                    while ( fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z, afield, bfield, cfield) ) {
                        //while ( fmr->atom->AtomInFragment(atnum, ifrag, istate, x, y, z) ) {
                        atnum++;
                    }

                    //BUGBUGBUGBUGBUGUBUGUBUGBUGBUGUBUGUBUGUBUG
                    //if statement for distance criteria between ifrag and loop frag
                    //if distance is greater than cutoff skip charges
                    //else do normal execution
                    //BUGBUGBUGBUGBUGUBUGUBUGBUGBUGUBUGUBUGUBUG

                    //printf("atnum:%d ifrag:%d istate:%d istate:%d x:%d y:%d z:%d iatom:%d gx:%d gy:%d gz:%d \n",atnum,ifrag,istate,x,y,z,iatom,gx,gy,gz);

                    if ( sscanf(line, "%lf %lf %lf", &gx, &gy, &gz) == 3 ) {
                        
                        if (fmr->atom->AtomInCell(atnum,istate,0,0,0, afield, bfield, cfield)) {
                            // gx,gy,gz = the electric field
                            // multiply by charge to get force (i.e. negative gradient) on atom
                            //double mmq = fmr->atom->getCharge(atnum%natoms, istate);
                            //gx *= -mmq;
                            //gy *= -mmq;
                            //gz *= -mmq;
                            //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
			    monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                            monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                            //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                        }
		    }
                    // Increment atnum for the next round
                    atnum++;
		}
                fclose(fs);
                         
/*                
                //determine COM of fragment I
                comi[0] = comi[1] = comi[2] = massi = 0.0;
                for (int i=0; i<natoms; ++i) {
                    if (atom->fragment[i] == ifrag) {
                        comi[0] += atom->mass[i] * (atom->coord[3*i+0] + x*cellA);
                        comi[1] += atom->mass[i] * (atom->coord[3*i+1] + y*cellB);
                        comi[2] += atom->mass[i] * (atom->coord[3*i+2] + z*cellC);
                        massi += atom->mass[i];
                    }
                }
                comi[0] /= massi;
                comi[1] /= massi;
                comi[2] /= massi;
                
                atnum = 0; // index of non-QM atom for storing gradient
                fgets(line, MAX_LENGTH, fs); // skip initial comment line
                
                for (int x0=-afield; x0<=afield; ++x0) {
                    for (int y0=-bfield; y0<=bfield; ++y0) {
                        for (int z0=-cfield; z0<=cfield; ++z0) {
                            
                            for (int jfrag=0; jfrag<nfragments; ++jfrag) {
                                
                                if (jfrag != ifrag || x0 != x || y0 != y || z0 != z) {
                                    
                                    //determine COM of fragment J
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
                                    
                                    
                                    //COM distance from fragment I to J
                                    d2 = 0;
                                    for (int l=0; l<3; ++l) {
                                        tmp = comj[l]-comi[l];
                                        d2 += tmp*tmp;
                                    }
                                    d=sqrt(d2);
                                    
                                    if (d>25.0) continue;
                                    
                                    for (int jatom=0; jatom<natoms; ++jatom) {
                                        if (atom->fragment[jatom] == jfrag) {
                                            
                                            fgets(line, MAX_LENGTH, fs);
                                            if ( sscanf(line, "%lf %lf %lf", &gx, &gy, &gz) == 3 ) {
                                                
                                                // gx,gy,gz = the electric field
                                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                                monomer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                                monomer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                                monomer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                                
                                            }
                                        }
                                        atnum++;
                                    }
                                }
                            }
                        }
                    }

                }
                fclose(fs);
*/
            }
            chdir("../");
        }
        ++index_mono;
    }
    
    // ***** Loop through dimers ***** //
    for (int idim=0; idim<n_dimers_sq; idim++) {
        
        if (run->dimer_queue[idim] == 0) continue;
        //if (run->dimer_list[idim] == -1) continue;
        
        if (ifrom_dim <= index_dim && index_dim < ito_dim) {
            
            //** new index variable **//
            int index = index_dim - ifrom_dim;
            
            int istate;
            int x,y,z;
            int ifrag,jfrag;
            
            //atom->getDimerIndices(run->dimer_list[idim],istate,x,y,z,ifrag,jfrag);
            atom->getDimerIndices(idim,istate,x,y,z,ifrag,jfrag);
            
            //printf("testing2:%2d rank:%d state:%d ifrag:%d jfrag;%d x:%d y:%d z:%d\n", idim, fmr->my_rank,istate, ifrag, jfrag, x, y, z);
            
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
            sprintf(directory, "%s/", state_directory);
            sprintf(directory, "%s/%s/", scratch_dir,state_directory);
            chdir(directory);
            
            //printf("%s %2d %4d %4d\n",jobname,fmr->my_rank,index_dim,idim);
            sprintf(command, "%s %s.nw > %s.nwout",
                    exec,
                    jobname,
                    jobname
                    );
            
            // ** The system call ** //
            ierr = system(command);
            
            // ** Check for error ** //
            //printf("ierr:%d\n",ierr);
            if (ierr) {
                printf("NWChem run error on rank %d:\n", fmr->my_rank);
                fmr->error(FLERR, command);
            }
            
            sprintf(command, "rm %s/%s/field_%s*",scratch_dir,state_directory,jobname);
            ierr = system(command);
            
            if (python) {
                // ** Open output file and get the energy ** //
                char output_file[MAX_LENGTH];
                //sprintf(output_file, "%s/%s.out", state_directory, filename);
                sprintf(output_file, "%s.nw.energy", jobname);
                FILE *fs = fopen(output_file, "r");
                if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
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
                        } else {
                            dimer_energies[index] = en;
                        }
                        //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                        
                    }
                }
                fclose(fs);
                
                if (FORCE) {
                    // ** Get gradient from file ** //
                    sprintf(output_file, "%s.nw.gradient", jobname);
                    fs = fopen(output_file, "r");
                    if (fs == NULL) {
                        char tmpstr[MAX_LENGTH];
                        sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
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
                        
                        if ( sscanf(line, "%lf %lf %lf", &gx, &gy, &gz) == 3 ) {
                            
                            if (fmr->atom->AtomInCell(atnum,istate,0,0,0)) {
                                //if (fmr->atom->AtomInFragment(atnum, jfrag, istate, 0, 0, 0) || fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0)) {
                                
                                // store symmetrically for zeroth unit cell
                                if (x==0 && y==0 && z==0) {
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                } else {
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                    dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                }
                            }
                        }
                        // Increment atnum for the next round
                        atnum++;
                    }
                    fclose(fs);
                }
            } else {
                
                // ** Open output file and get the energy ** //
                char output_file[MAX_LENGTH];
                sprintf(output_file, "%s.nwout",  jobname);
                FILE *fs = fopen(output_file, "r");
                if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                }
                char line[MAX_LENGTH];
                double en;
                double gw, gx, gy, gz;
                
                char tmp0[16],tmp1[16],tmp2[16],tmp3[16],tmp4[16];
                int atnum;
                atnum = 0; // index of QM atom for storing gradient
                
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    //printf("%s",line);
                    if (strcmp(run->correlation,"mp2") == 0) {
                        if ( strstr(line, "mp2 ENERGY GRADIENTS") ) {
                            
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            
                            for (int iatom=0; iatom<natoms; ++iatom) {
                                
                                fgets(line, MAX_LENGTH, fs);
                                if (strcmp(line,"\n") == 0) break;
                                
                                while ( !(fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0) ||
                                          fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z)) ) {
                                    atnum++;
                                }
                                //printf("%s",line);
                                if ( sscanf(line, "%lf %s %s %s %s %lf %lf %lf", &gw, tmp0, tmp1, tmp2, tmp3, &gx, &gy, &gz) == 8 ) {
                                    
                                    //printf(line);
                                    if (fmr->atom->AtomInCell(atnum,istate,0,0,0)) {
                                        //if (fmr->atom->AtomInFragment(atnum, jfrag, istate, 0, 0, 0) || fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0)) {
                                        
                                        // store symmetrically for zeroth unit cell
                                        if (x==0 && y==0 && z==0) {
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                        } else {
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                        }
                                    }
                                }
                                atnum++;
                            }
                            
                        } else if ( strstr(line, "Total MP2 energy") ) {
                            
                            //printf("%s",line);
                            if ( sscanf(line, "%s %s %s %lf", tmp0, tmp1, tmp2, &en) == 4 ) {
                                
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                // save symmetrized for zeroth unit cell
                                if (x==0 && y==0 && z==0) {
                                    dimer_energies[index] = en;
                                    dimer_energies[index] = en;
                                } else {
                                    dimer_energies[index] = en;
                                }
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                
                            }
                        }
                        
                    } else if (strcmp(run->correlation,"scf") == 0) {
                        if ( strstr(line, "RHF ENERGY GRADIENTS") ) {
                            
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            fgets(line, MAX_LENGTH, fs);
                            
                            for (int iatom=0; iatom<natoms; ++iatom) {
                                
                                fgets(line, MAX_LENGTH, fs);
                                if (strcmp(line,"\n") == 0) break;
                                
                                while ( !(fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0) ||
                                          fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z)) ) {
                                    atnum++;
                                }
                                //printf("%s",line);
                                if ( sscanf(line, "%lf %s %s %s %s %lf %lf %lf", &gw, tmp0, tmp1, tmp2, tmp3, &gx, &gy, &gz) == 8 ) {
                                    
                                    //printf(line);
                                    if (fmr->atom->AtomInCell(atnum,istate,0,0,0)) {
                                        //if (fmr->atom->AtomInFragment(atnum, jfrag, istate, 0, 0, 0) || fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0)) {
                                        
                                        // store symmetrically for zeroth unit cell
                                        if (x==0 && y==0 && z==0) {
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                        } else {
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                            dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                        }
                                    }
                                }
                                atnum++;
                            }
                            
                        } else if ( strstr(line, "Total SCF energy") ) {
                            
                            if ( sscanf(line, "%s %s %s %s %lf", tmp0, tmp1, tmp2, tmp3, &en) == 5 ) {
                                
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                // save symmetrized for zeroth unit cell
                                if (x==0 && y==0 && z==0) {
                                    dimer_energies[index] = en;
                                    dimer_energies[index] = en;
                                } else {
                                    dimer_energies[index] = en;
                                }
                                //BUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUGBUG
                                
                            }
                        }
                    }
                    
                }
                fclose(fs);
                
            }
            
            if (FORCE) {
                
                // ** Get field from file ** //
                char output_file[MAX_LENGTH];
                sprintf(output_file, "%s.nw.field", jobname);
                FILE *fs = fopen(output_file, "r");
                if (fs == NULL) {
                    char tmpstr[MAX_LENGTH];
                    sprintf(tmpstr, "Failure to read NWChem output file: %s", output_file);
                    fmr->error(FLERR, tmpstr);
                }
                char line[MAX_LENGTH];
                double gw, gx, gy, gz;
                int atnum = 0; // index of non-QM atom for storing gradient
                
                fgets(line, MAX_LENGTH, fs); // skip initial comment line
                while ( fgets(line, MAX_LENGTH, fs) != NULL ) {
                    // Advance atnum until it matches as a non-QM atom index for this dimer fragment
                    //while ( (fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0) ||
                    //         fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z)) ) {
                    while ( (fmr->atom->AtomInFragment(atnum, ifrag, istate, 0, 0, 0, afield, bfield, cfield) ||
                             fmr->atom->AtomInFragment(atnum, jfrag, istate, x, y, z, afield, bfield, cfield)) ) {
                        atnum++;
                    }
                    
                    if ( sscanf(line, "%lf %lf %lf", &gx, &gy, &gz) == 3 ) {
                        
                        if (fmr->atom->AtomInCell(atnum,istate,0,0,0,afield,bfield,cfield)) {
                            // gx,gy,gz = the electric field
                            // multiply by charge to get force (i.e. negative gradient) on atom
                            //double mmq = fmr->atom->getCharge(atnum%natoms, istate);
                            //gx *= -mmq;
                            //gy *= -mmq;
                            //gz *= -mmq;
                            // store symmetrically for zeroth unit cell
                            if (x==0 && y==0 && z==0) {
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                            } else {
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)]   = gx;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+1] = gy;
                                dimer_gradients[index*3*natoms + 3*(atnum%natoms)+2] = gz;
                            }
                        }
                    }
                    // Increment atnum for the next round
                    atnum++;
                }
                fclose(fs);
            }
            chdir("../");
            
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
                            //if (monomer_queue[midx] == 1) {
                            //			    printf("monomer_idx:%8d\n",midx);
                            men1=0.0;

                            for (int i=0; i<mdiv+1; i++) {
                                //				printf("proc:%d\n",monomer_proc[i]);
                                if (monomer_proc[i] == midx) { men1 = monomer_energies[i]; break; }
                            }

                            if (fmr->print_level > 1) {
                                fprintf(fs,"Rank %3d - FRAG I:%4d | %16.10f\n", fmr->my_rank,
                                        ifrag, men1);
                            }
                            
                            en_fmo1 += men1;
                            //}
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
                                    if (dimer_proc[i]==idx) { den = dimer_energies[i];
                                        //printf("ddiv:%d state:%d rank:%d  didx:%4d  idxi:%4d  idxj:%4d en:%f\n",i,istate,fmr->my_rank,idx,idxi,idxj,den);
                                        break; }
                                }
                                for (int i=0; i<mdiv+1; i++) {
                                    if (monomer_proc[i] == idxi) men1 = monomer_energies[i];
                                    if (monomer_proc[i] == idxj) men2 = monomer_energies[i];
                                }


                                double etmp = den - men1 - men2;
                                
                                if (fmr->print_level > 1) {
                                    fprintf(fs,"Didx:%2d  idxi:%2d  idxj:%2d  State:%d Rank %3d - FRAG I:%02d--J:%02d  CELL x:%2d y:%2d z:%2d | %16.10f %16.10f %16.10f %16.10f\n", idx, idxi, idxj,istate, fmr->my_rank, ifrag, jfrag,x,y,z,
                                            den, men1, men2, etmp);
                                }

                                en_fmo2 += etmp*oc;
                            }

                        }
                    }
                    
                }
            }
        }
        fmo_energies[istate] = en_fmo1 + en_fmo2;
        //printf("state: %d FMO ener BUGBUG: %f\n",istate,fmo_energies[istate]);
        
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
                                //printf("midx:%d\n",midx);
                                //if (monomer_queue[midx] == 1) {
                                
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
                                    
                                    //for (int i=0; i<natoms; ++i) {
                                        
                                        //printf("%lf %lf %lf\n",fmo_gradients[3*natoms*istate + 3*i+0],fmo_gradients[3*natoms*istate + 3*i+1],fmo_gradients[3*natoms*istate + 3*i+2]);
                                        //fmo_gradients[3*natoms*istate + 3*i+0] += gx*oc;
                                        //fmo_gradients[3*natoms*istate + 3*i+1] += gy*oc;
                                        //fmo_gradients[3*natoms*istate + 3*i+2] += gz*oc;
                                        //printf("%lf %lf %lf\n",gx,gy,gz);
                                        
                                        
                                    //}
                                    
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
    
    // *** Reduce the energies from the parallel NWChem calls *** //
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
    /*
     if (fmr->my_rank==1) {
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
    
    // Clock
    MPI_Barrier(fmr->world);
    double comp_end = MPI_Wtime();
    if (fmr->master_rank) {
        printf("Time for FMO compilation: %.4f seconds\n\n", comp_end - comp_start);
    }
}

