/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "cec.h"
#include "umbrella.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Umbrella::Umbrella(FMR *fmr) : Pointers(fmr)
{
    // Basic initializations
    deltaE     = NULL;
    koppa      = NULL;
    upsilon    = NULL;
    
    umb_typ    = 0;

}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Umbrella::~Umbrella()
{
    int nstates = fmr->atom->nstates;
    
    if (deltaE != NULL) delete [] deltaE;
    
    if (koppa != NULL) {
	    for (int i=0; i<nstates; ++i) {
            delete [] koppa[i];
        }
        delete [] koppa;
    }
    if (upsilon != NULL) {
	    for (int i=0; i<nstates; ++i) {
            delete [] upsilon[i];
        }
        delete [] upsilon;
    }

}

void Umbrella::decompose_force(double* force)
{
    
    Atom *atom	       = fmr->atom;
    Matrix *matrix     = fmr->matrix;
    
    int natoms         = atom->natoms;
    int nstates        = atom->nstates;

    double *GSGradient = fmr->atom->force;
    double *GSCoeffs   = fmr->matrix->GSCoeffs;
    
    /******************************************************************/
    /*** Calculate derivitive of [qCOC(i)] ****************************/
    /******************************************************************/
    if (fmr->master_rank) {
        for(int istate=0; istate<nstates; ++istate)
        {
            double C2 = GSCoeffs[istate]*GSCoeffs[istate];
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (atom->reactive[istate*natoms + iatom]) {

                    for(int l=0; l<3; l++) {
                       GSGradient[3*iatom + l] += force[l]*C2;//*(atom->getCharge(iatom,istate)/qsum_coc[istate]);
                    }      
                }
            } // loop all the atoms in the coc
        } // loop all the coc(s)
    }
    MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world); 

    /******************************************************************/
    /*** Calculate derivitive of [C(i)^2]  ****************************/
    /******************************************************************/
    partial_C_N2(force);    
}

void Umbrella::partial_C_N2(double *force)
{
    /******************************************************************/
    /*** JPCB, 112, 2349 Eq 24-26 *************************************/
    /******************************************************************/
    
    Atom *atom	       = fmr->atom;
    Matrix *matrix     = fmr->matrix;
    Cec *cec           = fmr->cec;
    
    int natoms         = atom->natoms;
    int nstates        = atom->nstates;
    int prev_nstates   = atom->prev_nstates;
    
    double *GSGradient = fmr->atom->force;
    double **Evecs     = matrix->Evecs;
    double *Evals      = matrix->Evals;
    double ***HX       = matrix->HX;
    double **r_coc     = cec->r_coc;
    
    int ground = 0; // convention retained from matrix diagonalization
    double factor;
   
    // ** Allocate and assign deltaE array ** //
    if (deltaE != NULL) delete [] deltaE;
    deltaE = new double [nstates];
    for (int i=0; i<nstates; ++i) deltaE[i] = Evals[ground]-Evals[i];
    
    // ** Allocate and assign summation vectors array ** //
    if (koppa != NULL || upsilon != NULL) {
        for (int i=0; i<prev_nstates; ++i) {
            delete [] koppa[i];
            delete [] upsilon[i];
        }
        delete [] koppa;
        delete [] upsilon;
    }
    koppa   = new double*[nstates];
    upsilon = new double*[nstates];
    for (int i=0; i<nstates; ++i) {
        koppa[i]   = new double[3];
        upsilon[i] = new double[3];
        for (int k=0; k<3; ++k) {
            koppa[i][k]   = 0.0;
            upsilon[i][k] = 0.0;
        }
    }
    
    /********** Eq. 24 ***************/
    if (fmr->master_rank) {
        
        for(int j=0; j<nstates; j++)
        {
            if(j==ground) continue;
            
            for(int i=0; i<nstates; i++)
            {
                factor = Evecs[i][j]*Evecs[i][ground];
                koppa[j][0] += 2*force[0]*(factor*r_coc[i][0]);
                koppa[j][1] += 2*force[1]*(factor*r_coc[i][1]);
                koppa[j][2] += 2*force[2]*(factor*r_coc[i][2]);
            }
        }
        
        /********** Eq. 25 ***************/
        
        for(int l=0; l<nstates; l++)
        {
            for(int j=0; j<nstates; j++)
            {
                if(j==ground) continue;
                
                factor = Evecs[l][j]/deltaE[j];
                upsilon[l][0] += factor*koppa[j][0];
                upsilon[l][1] += factor*koppa[j][1];
                upsilon[l][2] += factor*koppa[j][2];
            }
            
        }
        
        /*********** Eq. 26 ******************/
        
        for (int l=0; l<nstates; l++) {
            for (int m=0; m<nstates; m++) {
                
                double prefactor[3];
                prefactor[0] = Evecs[m][ground]*upsilon[l][0];
                prefactor[1] = Evecs[m][ground]*upsilon[l][1];
                prefactor[2] = Evecs[m][ground]*upsilon[l][2];
                
                for (int iatom=0; iatom<natoms; iatom++) {
                    
                    for (int x=0; x<3; x++) {
                        GSGradient[3*iatom + x] -= HX[l][m][3*iatom + x]*prefactor[x];
                    }
                }
            }
        }
    }
    MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

void Umbrella::decompose_energy(double energy)
{
    Matrix *matrix  = fmr->matrix;
    
    double GSEnergy = matrix->GSEnergy;  
  
    if (fmr->master_rank) GSEnergy += energy;
            
    MPI_Bcast(&GSEnergy, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

void Umbrella::compute()
{
    Atom *atom	     = fmr->atom;
    Cec *cec	     = fmr->cec;

    int natoms       = atom->natoms;
    int nstates      = atom->nstates;

    double *r_cec    = cec->r_cec;
    
    if (fmr->master_rank) {
        
        di[0]=di[1]=di[2]=1;
        center[0]=center[1]=center[2]=0.0;
        ref[0]=ref[1]=ref[2]=1.0/sqrt(3.0);
        k[0]=k[1]=k[2]=1.0*fmr->math->kcal2au;
        
        energy = 0.0;
        f[0][0] = f[0][1] = f[0][2] = f[1][0] = f[1][1] = f[1][2] = 0.0;
        virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;
        
        dx[0] = dx[1] = dx[2] = 0.0;
        ff[0] = ff[1] = ff[2] = 0.0;
        
        if(umb_typ==COORD_CART)
        {
            printf("\nCOORD_CART\n");
            for(int i=0; i<3; i++) if (di[i]) dx[i] = center[i]-r_cec[i];
            printf("dx %f %f %f\n",dx[0],dx[1],dx[2]);
            //VECTOR_PBC(dx);
            //for(int i=0; i<3; i++) if (di[i]) dx[i] = dx[i]-ref[i];
            //VECTOR_PBC(dx);
            
            for(int i=0; i<3; i++) if (di[i])
            {
                ff[i] = -k[i];
                f[0][i] = - k[i] * dx[i]; f[1][i] = -f[0][i];
                dx2[i] = dx[i]*dx[i];
                energy += 0.5 * k[i] * dx2[i];
            }
            
            for(int i=0; i<3; i++) if (di[i]) dx[i] = center[i]-r_cec[i];
            //VECTOR_PBC(dx);
        }
        
        // Virial calculation
        virial[0] += dx[0]*dx[0]*ff[0];
        virial[1] += dx[1]*dx[1]*ff[1];
        virial[2] += dx[2]*dx[2]*ff[2];
        
    }
    decompose_energy(energy);
    decompose_force(f[0]);

    //MPI_Bcast(&fmr->matrix->GSEnergy,  1, MPI_DOUBLE, MASTER_RANK, fmr->world);
    //MPI_Bcast(fmr->atom->force, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

void Umbrella::write_log()
{

}

