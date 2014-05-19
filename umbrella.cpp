/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "cec.h"
#include "umbrella.h"
#include "run.h"
#include "matrix.h"
#include "dynamics.h"
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
    beta       = NULL;
    omega      = NULL;
    
    k[0]       = 0.0;
    k[1]       = 0.0;
    k[2]       = 0.0;
    
    center[0]  = 0.0;
    center[1]  = 0.0;
    center[2]  = 0.0;
    
    di[0]      = 0;
    di[1]      = 0;
    di[2]      = 0;

    ref[0]     = 0.0;
    ref[1]     = 0.0;
    ref[2]     = 0.0;
    
    umb_coord  = 0;
    do_umbrella_sampling = 0;
    
    sprintf(umb_typ, "%s", "null");
    sprintf(umbFile, "%s", "fmr_umb.log"); // default file name

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
    if (omega != NULL) {
	    for (int i=0; i<nstates; ++i) {
            delete [] omega[i];
        }
        delete [] omega;
    }
    if (beta != NULL) {
	    for (int i=0; i<nstates; ++i) {
            delete [] beta[i];
        }
        delete [] beta;
    }
    

}

void Umbrella::decompose_force(double* force)
{
    
    Atom *atom	       = fmr->atom;
    Matrix *matrix     = fmr->matrix;
    Cec *cec           = fmr->cec;
    
    int natoms         = atom->natoms;
    int nstates        = atom->nstates;

    double *GSGradient = atom->force;
    double *GSCoeffs   = matrix->GSCoeffs;
    
    double *qsum_coc   = cec->qsum_coc;
    
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
                       //GSGradient[3*iatom + l] += force[l]*C2;//*(atom->getCharge(iatom,istate)/qsum_coc[istate]);
                       GSGradient[3*iatom + l] += force[l]*C2*fabs(atom->getCharge(iatom,istate))/qsum_coc[istate];
                    }      
                }
            } // loop all the atoms in the coc
        } // loop all the coc(s)
    }
    MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world); 

    /******************************************************************/
    /*** Calculate derivitive of [C(i)^2]  ****************************/
    /******************************************************************/
    //partial_C_N4(force);
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

void Umbrella::partial_C_N3(double *force)
{
    /******************************************************************/
    /*** JPCB, 112, 2349 Eq 21-23 *************************************/
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
    if (beta != NULL || omega != NULL) {
        for (int i=0; i<prev_nstates; ++i) {
            delete [] beta[i];
            delete [] omega[i];
        }
        delete [] beta;
        delete [] omega;
    }
    beta   = new double*[nstates];
    omega = new double*[nstates];
    for (int i=0; i<nstates; ++i) {
        beta[i]   = new double[3];
        omega[i] = new double[3];
        for (int k=0; k<3; ++k) {
            beta[i][k]   = 0.0;
            omega[i][k] = 0.0;
        }
    }
    
    /********** Eq. 21 ***************/
    if (fmr->master_rank) {
        
        for (int j=0; j<nstates; j++)
        {
            if (j==ground) continue;
            
            for (int l=0; l<nstates; l++)
            {
                for (int m=0; m<nstates; m++)
                {
                    factor = Evecs[l][j]*Evecs[m][ground];
                    for (int iatom=0; iatom<natoms; iatom++) {
                        
                        beta[j][0] += 2*force[0]*(factor*HX[m][l][3*iatom + 0]);
                        beta[j][1] += 2*force[1]*(factor*HX[m][l][3*iatom + 1]);
                        beta[j][2] += 2*force[2]*(factor*HX[m][l][3*iatom + 2]);
                    }
                }
            }
        }
        
        /********** Eq. 22 ***************/
        
        for(int i=0; i<nstates; i++)
        {
            for(int j=0; j<nstates; j++)
            {
                if(j==ground) continue;
                
                factor = Evecs[i][j]/deltaE[j];
                omega[i][0] += factor*beta[j][0];
                omega[i][1] += factor*beta[j][1];
                omega[i][2] += factor*beta[j][2];
            }
        }
        
        /*********** Eq. 23 ******************/
        
        for (int i=0; i<nstates; i++)
        {
            double prefactor[3];
            prefactor[0] = Evecs[i][ground]*r_coc[i][0]*omega[i][0];
            prefactor[1] = Evecs[i][ground]*r_coc[i][1]*omega[i][1];
            prefactor[2] = Evecs[i][ground]*r_coc[i][2]*omega[i][2];
            
            for (int iatom=0; iatom<natoms; iatom++) {
                
                for (int x=0; x<3; x++) {
                    GSGradient[3*iatom + x] -= prefactor[x];
                }
            }
        }
    }
    MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world);
}

void Umbrella::partial_C_N4(double *force)
{
    /******************************************************************/
    /*** JPCB, 112, 2349 Eq 20 ****************************************/
    /******************************************************************/
    
    Atom *atom	       = fmr->atom;
    Matrix *matrix     = fmr->matrix;
    Cec *cec           = fmr->cec;
    
    int natoms         = atom->natoms;
    int nstates        = atom->nstates;
    
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
    
    /********** Eq. 20 ***************/
    if (fmr->master_rank) {
        
        for (int i=0; i<nstates; i++)
        {
            for (int j=0; j<nstates; j++)
            {
                if (j==ground) continue;
                
                for (int l=0; l<nstates; l++)
                {
                    for (int m=0; m<nstates; m++)
                    {
                        factor = Evecs[m][ground]*Evecs[l][j]*Evecs[i][j]*Evecs[i][ground]/deltaE[j];
                        for (int iatom=0; iatom<natoms; iatom++) {
                            
                            GSGradient[3*iatom + 0] -= 2*force[0]*HX[m][l][3*iatom + 0]*factor*r_coc[i][0];
                            GSGradient[3*iatom + 1] -= 2*force[1]*HX[m][l][3*iatom + 1]*factor*r_coc[i][1];
                            GSGradient[3*iatom + 2] -= 2*force[2]*HX[m][l][3*iatom + 2]*factor*r_coc[i][2];
                        }
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

void Umbrella::setup()
{

    Atom *atom         = fmr->atom;
    Cec *cec           = fmr->cec;
    Dynamics *dynamics = fmr->dynamics;

    int natoms         = atom->natoms;
    int nstates        = atom->nstates;
    
    double sum         = 0.0;

    if (fmr->master_rank) printf("Setting up umbrella sampling...\n");
    
    if (umb_coord==COORD_SPHERICAL) {
        
        k[1] = k[0];
        k[2] = k[0];
        
        for (int i=0; i<3; i++) {
            sum   += ref[i];
            di[i]  = 1;
        }
        for (int i=0; i<3; i++) ref[i] /= sum;

    }
    else if (umb_coord==COORD_CARTESIAN) {
        
        for (int i=0; i<3; i++) {
            
            sum += ref[i];
            if (k[i]) di[i] = 1;
            
        }
        for (int i=0; i<3; i++) ref[i] /= sum;


    }
    else if (umb_coord==COORD_CYLINDRICAL) {
        
        for (int i=0; i<3; i++) sum += ref[i];
        for (int i=0; i<3; i++) ref[i] /= sum;
    }

}

void Umbrella::compute()
{
    Atom *atom	       = fmr->atom;
    Cec *cec	       = fmr->cec;
    Dynamics *dynamics = fmr->dynamics;

    int natoms         = atom->natoms;
    int nstates        = atom->nstates;

    double *r_cec      = cec->r_cec;
    
    int iCurrentStep   = dynamics->iCurrentStep;
    
    if (fmr->master_rank) printf("Computing umbrella sampling...\n");
    
    if (fmr->master_rank) {

        energy = 0.0;
        diff   = 0.0;
        
        f[0][0]   = f[0][1]   = f[0][2]   = f[1][0]   = f[1][1]   = f[1][2]   = 0.0;
        virial[0] = virial[1] = virial[2] = virial[3] = virial[4] = virial[5] = 0.0;
        
        dx[0] = dx[1] = dx[2] = 0.0;
        ff[0] = ff[1] = ff[2] = 0.0;
        
        if (umb_coord==COORD_CARTESIAN)
        {
            for(int i=0; i<3; i++) if (di[i]) dx[i] = center[i]-r_cec[i];
            //VECTOR_PBC(dx);
            
            for (int i=0; i<3; i++) if (di[i])
            {
                ff[i]    = -k[i];
                f[0][i]  = -k[i] * dx[i];
                f[1][i]  = -f[0][i];
                dx2[i]   = dx[i] * dx[i];
                energy  += 0.5 * k[i] * dx2[i];
                diff    += dx[i] * ref[i]; 
            }
        }
        
        // Virial calculation
        virial[0] += dx[0]*dx[0]*ff[0];
        virial[1] += dx[1]*dx[1]*ff[1];
        virial[2] += dx[2]*dx[2]*ff[2];
        
    }
    
    if (di[0] || di[0] || di[2]) {
       
      decompose_energy(energy);
      decompose_force(f[0]);

      writeStepUmb(iCurrentStep);
    }
}

/*-----------------------------------------------------------------
 Report information about umbrella sampling for this step
 -----------------------------------------------------------------*/
void Umbrella::writeStepUmb(int mode)
{
    Atom *atom	       = fmr->atom;
    Matrix *matrix     = fmr->matrix;
    Cec *cec           = fmr->cec;
        
    double *r_cec     = cec->r_cec;
    
    if (fmr->master_rank) {
        
        FILE *fs;
        
        if (mode == 0) fs = fopen(umbFile, "w");
        else           fs = fopen(umbFile, "a");
        if (fs == NULL) {
            char tmpstr[256];
            sprintf(tmpstr, "Failure to write to umbrella file %s", umbFile);
            fmr->error(FLERR, tmpstr);
        }
        
        // Print header
        if (mode == 0) fprintf(fs, "# %10s %15s %15s %15s %15s %15s %15s %15s %15s\n",
                                   "StepNumber","Ener[kcal]","Diff[ang.]",
                                   "DX[ang.]","X2[ang.]",
                                   "DY[ang.]","Y2[ang.]",
                                   "DZ[ang.]","Z2[ang.]");
        
        fprintf(fs, "  %10d %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f %15.10f\n",
                mode, energy*fmr->math->au2kcal, diff,
                dx[0],r_cec[0],
                dx[1],r_cec[1],
                dx[2],r_cec[2]);
        
        fclose(fs);
    }
}

