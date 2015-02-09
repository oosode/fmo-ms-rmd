/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "electrostatics.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Electrostatics::Electrostatics(FMR *fmr) : Pointers(fmr)
{
/*
    // Basic initializations
    qsum_coc   = NULL;
    r_coc      = NULL;
*/
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Electrostatics::~Electrostatics()
{
/*
    int nstates = fmr->atom->nstates;
    
    if (r_coc != NULL) {
        for (int i=0; i<nstates; ++i)
            delete [] r_coc[i];
        delete [] r_coc;
    }

    if (qsum_coc != NULL) delete [] qsum_coc;
*/
}

/*-----------------------------------------------------------------
  State search algorithm 
-----------------------------------------------------------------*/
void Electrostatics::coulomb()
{

    Run *run       = fmr->run;
    Atom *atom     = fmr->atom;
    int natoms     = atom->natoms;
    int nstates    = fmr->atom->nstates;

    int cellA      = fmr->atom->cellA;
    int cellB      = fmr->atom->cellB;
    int cellC      = fmr->atom->cellC;

    int xa         = fmr->atom->na;
    int xb         = fmr->atom->nb;
    int xc         = fmr->atom->nc;

    int afield     = fmr->atom->afield;
    int bfield     = fmr->atom->bfield;
    int cfield     = fmr->atom->cfield;

    double fxi,fyi,fzi,duij,uij,xij,yij,zij,rij,rinv,qi,qj;

    if (fmr->master_rank) {

        for (int istate=0; istate<nstates; ++istate) {
  
            for (int x=-15; x<=15; ++x) {
                for (int y=-15; y<=15; ++y) {
                    for (int z=-15; z<=15; ++z) {

	    	        if ( xa <= afield && xb <= bfield && xc <= cfield ) continue;
 
                        for (int i=0; i<natoms; ++i) {
                            for (int j=0; j<natoms; ++j) {
                           
                                qi = atom->getCharge(i, istate);
                                qj = atom->getCharge(j, istate);
			
  	    	    	        xij = atom->coord[3*i+0] - atom->coord[3*j+0] + x*cellA;
			        yij = atom->coord[3*i+1] - atom->coord[3*j+1] + y*cellB;
			        zij = atom->coord[3*i+2] - atom->coord[3*j+2] + z*cellC;

                                rij = sqrt( xij*xij + yij*yij + zij*zij );
                                rinv = 1.0/rij;

                                uij = qi*qj*rinv;
				std::cout << uij << std::endl;
		
			        run->fmo_energies[istate] += uij;
			
			        duij = -uij*rinv;

                                fxi = - duij*xij*rinv;
                                fyi = - duij*yij*rinv;
                                fzi = - duij*zij*rinv;

			        run->fmo_gradients[3*natoms*istate + 3*i+0] += fxi;
			        run->fmo_gradients[3*natoms*istate + 3*i+1] += fyi;
			        run->fmo_gradients[3*natoms*istate + 3*i+2] += fzi;

                                run->fmo_gradients[3*natoms*istate + 3*j+0] -= fxi;
                                run->fmo_gradients[3*natoms*istate + 3*j+1] -= fyi;
                                run->fmo_gradients[3*natoms*istate + 3*j+2] -= fzi;
			    
			    }
		        }
		    }
	        }
	    }
        }
    }
}

void Electrostatics::ewald()
{



}

