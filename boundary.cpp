/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "boundary.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Boundary::Boundary(FMR *fmr) : Pointers(fmr)
{

    do_boundary_conditions = 0;
    
    k[0] = k[1] = k[2] = 0.0;
    
    center[0] = center[1] = center[2] = 0.0;
    
    radius[0] = radius[1] = radius[2] = 0.0;
    vector[0] = vector[1] = vector[2] = 0.0;
     
    bound_coord  = 0;
    
}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Boundary::~Boundary()
{


    
}

void Boundary::setup()
{
    
    Atom *atom         = fmr->atom;
    
    int natoms         = atom->natoms;
    //int nstates        = atom->nstates;
    
    if (fmr->master_rank)  printf("Setting up boundary conditions...\n");
    
    // center of mass of i fragment
    double mass = 0.0;
    for (int i=0; i<natoms; ++i) {
        center[0] += atom->mass[i]/fmr->math->amu2au * atom->coord[3*i];
        center[1] += atom->mass[i]/fmr->math->amu2au * atom->coord[3*i+1];
        center[2] += atom->mass[i]/fmr->math->amu2au * atom->coord[3*i+2];
    }
    center[0] /= atom->totalMass;
    center[1] /= atom->totalMass;
    center[2] /= atom->totalMass;

    if (bound_coord==COORD_SPHERICAL) {
        
        k[1] = k[0];
        k[2] = k[0];
        
        radius[1] = radius[0];
        radius[2] = radius[0];
        
        
    } else if ( bound_coord==COORD_CYLINDRICAL) {
        
	
        k[2] = k[0];

        radius[2] = radius[0];       
 
    }
}

void Boundary::compute()
{
    
    Atom *atom	       = fmr->atom;
    Matrix *matrix     = fmr->matrix;
    
    double GSEnergy    = matrix->GSEnergy;
    double *GSGradient = fmr->atom->force;
    double *coord      = atom->coord;
    int natoms         = atom->natoms;
    int nstates        = atom->nstates;
    
    if (fmr->master_rank) printf("Computing boundary conditions...\n");
    
    virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;
    
    dx[0] = dx[1] = dx[2] = 0.0;
    ff[0] = ff[1] = ff[2] = 0.0;
    
    if (bound_coord==COORD_SPHERICAL)
    {
        for (int iatom=0; iatom<natoms; ++iatom) {
            
            double energy = 0.0;
            double sum = 0.0;
            double distance = 0.0;
            
            for (int i=0; i<3; i++) {
                di[i] = coord[3*iatom+i]-center[i];
            }
		
            distance = sqrt(di[0]*di[0] + di[1]*di[1] + di[2]*di[2]);
            if (distance > radius[0]) {
                
                diff = distance - radius[0];
                diff2 = diff*diff;
                energy = k[0]*diff2;
                
                for (int i=0; i<3; i++) {
                    
                    di[i] /= distance;
                    
                    f[0][i] = (2 * k[0] * diff) * di[i];
                    f[1][i] = -f[0][i];
                    
                }
//		if (fmr->master_rank) printf("d: %lf, diff: %lf, %lf %lf %lf\n",distance,diff,di[0],di[1],di[2]);
            } else {

		energy = 0.0;
                f[0][0] = f[0][1] = f[0][2] = f[1][0] = f[1][1] = f[1][2] = 0.0;
	    }

            if (fmr->master_rank) {
                GSEnergy += energy;
                GSGradient[3*iatom + 0] += f[0][0];
                GSGradient[3*iatom + 1] += f[0][1];
                GSGradient[3*iatom + 2] += f[0][2];
//		printf("%2d %lf %lf %lf %lf\n",iatom,energy,f[0][0],f[0][1],f[0][2]);
            }
        }
        
    } else if (bound_coord==COORD_CYLINDRICAL) {
        
        
    }
    
	MPI_Bcast(&GSEnergy, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);
	MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
    
}
