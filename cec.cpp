/* AWGL */
#include "fmr.h"
#include "fmr_math.h"
#include "atom.h"
#include "state.h"
#include "cec.h"
#include "run.h"
#include "matrix.h"
#include <vector>
#include <string>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Cec::Cec(FMR *fmr) : Pointers(fmr)
{
    // Basic initializations
    qsum_coc   = NULL;
    r_coc      = NULL;
    deltaE     = NULL;
    koppa      = NULL;
    upsilon    = NULL;
    
    umb_typ    = 0;

}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Cec::~Cec()
{
    int nstates = fmr->atom->nstates;
    
    if (r_coc != NULL) {
        for (int i=0; i<nstates; ++i)
            delete [] r_coc[i];
        delete [] r_coc;
    }

    if (qsum_coc != NULL) delete [] qsum_coc;
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

/*-----------------------------------------------------------------
  State search algorithm 
-----------------------------------------------------------------*/
void Cec::compute_coc()
{
    
    Atom *atom	 = fmr->atom;
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;
    int prev_nstates = atom->prev_nstates;
    
    // ** Allocate qsum_coc array ** //
    // If qsum_coc is already allocated from previous step, de-allocate
    if (qsum_coc != NULL) delete [] qsum_coc;
    
	// Allocate qsum_coc array based on number of states
    qsum_coc = new double [nstates];
    for (int i=0; i<nstates; ++i) qsum_coc[i] = 0.0;
    
    // ** Allocate r_coc array ** //
    // If r_coc is already allocated from previous step, de-allocate
	if (r_coc != NULL) {
	    for (int i=0; i<prev_nstates; ++i) {
            delete [] r_coc[i];
        }
        delete [] r_coc;
    }
    
    // Allocate r_coc array based on number of states for this step
    r_coc = new double*[nstates];
    for (int i=0; i<nstates; ++i) {
        r_coc[i] = new double[3];
        r_coc[i][0] = 0.0;
        r_coc[i][1] = 0.0;
        r_coc[i][2] = 0.0;
    }
    
    if (fmr->master_rank) {
        
        double ref[3];
        /********************************************/
        /*** Assign COC information  ****************/
        /********************************************/
        for (int istate=0; istate<nstates; ++istate) {
            int atoms_per_ireactive = 0;
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (atom->reactive[istate*natoms + iatom]) {
                    
                    if (atom->symbol[iatom] == 'O') {
                        ref[0] = atom->coord[3*iatom];
                        ref[1] = atom->coord[3*iatom + 1];
                        ref[2] = atom->coord[3*iatom + 2];
                    }
                    atoms_per_ireactive++;
                    qsum_coc[istate] += atom->getCharge(iatom,istate);
                    
                }
            }
            
            //printf("reference state %d: %f %f %f\n",istate,ref[0],ref[1],ref[2]);
            //natom_coc[istate] = atoms_per_ireactive;
            for (int iatom=0; iatom<natoms; ++iatom) {
                if (atom->reactive[istate*natoms + iatom]) {
                    
                    double rr[3];
                    double dr[3];
                    
                    rr[0] = atom->coord[3*iatom];
                    rr[1] = atom->coord[3*iatom + 1];
                    rr[2] = atom->coord[3*iatom + 2];
                    
                    //printf("current position: %f %f %f\n",rr[0],rr[1],rr[2]);
                    
                    VECTOR_SUB(dr,rr,ref);
                    //VECTOR_PBC(dr);
                    
                    //printf("charge: %f/%f = %f\n",atom->getCharge(iatom,istate),qsum_coc[istate],atom->getCharge(iatom,istate)/qsum_coc[istate]);
                    VECTOR_SCALE(dr,atom->getCharge(iatom,istate)/qsum_coc[istate]);
                    VECTOR_ADD(r_coc[istate],r_coc[istate],dr);
                    //printf("coc: %f %f %f\n",r_coc[istate][0],r_coc[istate][1],r_coc[istate][2]);
                    
                }
            }
            VECTOR_ADD(r_coc[istate],r_coc[istate],ref);
            printf("COC position for state %2d: %15.10f %15.10f %15.10f\n",istate,r_coc[istate][0],r_coc[istate][1],r_coc[istate][2]);
            
        }
    }
}

void Cec::compute_cec()
{

  double *GSCoeffs = fmr->matrix->GSCoeffs;
  int natoms       = atom->natoms;
  int nstates      = atom->nstates;

  if (fmr->master_rank) {

    r_cec[0]=r_cec[1]=r_cec[2]=0.0;
    
    double ref[3];
    ref[0] = r_coc[0][0];
    ref[1] = r_coc[0][1];
    ref[2] = r_coc[0][2];

    for (int istate=1; istate<nstates; istate++) {

      double dr[3];
      double C2 = GSCoeffs[istate]*GSCoeffs[istate];

      VECTOR_SUB(dr,r_coc[istate],ref);
      //VECTOR_PBC(dr);
      VECTOR_SCALE_ADD(r_cec,dr,C2);

    }
    VECTOR_ADD(r_cec,r_cec,ref);

    printf("CEC position:              %15.10f %15.10f %15.10f\n",r_cec[0],r_cec[1],r_cec[2]);

  }
}


void Cec::decompose_force2(double* force)
{
    
    Atom *atom	     = fmr->atom;
    Matrix *matrix   = fmr->matrix;
    
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;

    double *GSGradient = fmr->atom->force;
    double *GSCoeffs   = fmr->matrix->GSCoeffs;
    
    //EVB_Matrix* matrix;
    //if(evb_engine->ncomplex==1) matrix = (EVB_Matrix*)(evb_engine->full_matrix);
    //else matrix = (EVB_Matrix*)(evb_engine->all_matrix[cplx->id-1]);
    
    /******************************************************************/
    /*** Calculate derivitive of [qCOC(i)] ****************************/
    /******************************************************************/
    //printf("Calculate derivitive of [qCOC(i)]\n"); 
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
    
    /******************************************************************/
    /*** Calculate derivitive of [C(i)^2]  ****************************/
    /******************************************************************/
    
    //partial_C_N3(force);
    //printf("Calculate derivitive of [C(i)^2]\n");
    partial_C_N22(force);
    
    MPI_Bcast(GSGradient, 3*natoms, MPI_DOUBLE, MASTER_RANK, fmr->world); 

}

void Cec::partial_C_N22(double *force)
{
    /******************************************************************/
    /*** JPCB, 112, 2349 Eq 24-26 *************************************/
    /******************************************************************/
    
    Atom *atom	     = fmr->atom;
    Matrix *matrix   = fmr->matrix;
    
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;
    int prev_nstates = atom->prev_nstates;
    
    double *GSGradient = fmr->atom->force;
    //double *GSCoeffs   = matrix->GSCoeffs;
    double **Evecs     = matrix->Evecs;
    double *Evals      = matrix->Evals;
    double ***HX       = matrix->HX;
    
    int ground = 0; // convention retained from matrix diagonalization
    double factor;
   
    //int nall = atom->nlocal + atom->nghost;
    //EVB_Matrix* matrix = evb_matrix;
    
    //int *parent = cplx->parent_id;
    //double ***diagonal = matrix->f_diagonal;
    //double ***off_diagonal = matrix->f_off_diagonal;
    //double ***extra_coupl = matrix->f_extra_coupling;
   
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
    
    printf("decompose_force %f %f %f\n\n",force[0],force[1],force[2]);
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

void Cec::decompose_energy2(double energy)
{
    Atom *atom	     = fmr->atom;
    Matrix *matrix   = fmr->matrix;
    
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;
    
    double GSEnergy   = fmr->matrix->GSEnergy;
    
    if (fmr->master_rank) {
     
        GSEnergy += energy;
        
    }
    MPI_Bcast(&GSEnergy, 1, MPI_DOUBLE, MASTER_RANK, fmr->world);
    
}

void Cec::compute2()
{
    Atom *atom	     = fmr->atom;
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;
    
    if (fmr->master_rank) {
        
        di[0]=di[1]=di[2]=1;
        center[0]=center[1]=center[2]=0.0;
        ref[0]=ref[1]=ref[2]=1.0/sqrt(3.0);
        k[0]=k[1]=k[2]=100.0*fmr->math->kcal2au;
        
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
            //printf("dx %f %f %f\n",dx[0],dx[1],dx[2]);
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
        
        decompose_energy2(energy);
        decompose_force2(f[0]);
        
    }
}

