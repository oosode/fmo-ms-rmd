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

    if (qsum_coc    != NULL) delete [] qsum_coc;

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
    }
    
	// ** Initialize r_coc array to zero ** //
    for (int i=0; i<nstates; ++i)
    {
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

void Cec::compute()
{

  double *C2     = fmr->matrix->GSCoeffs;
  int natoms     = atom->natoms;
  int nstates    = atom->nstates;

  if (fmr->master_rank) {

    r_cec[0]=r_cec[1]=r_cec[2]=0.0;
    
    double ref[3];
    ref[0] = r_coc[0][0];
    ref[1] = r_coc[0][1];
    ref[2] = r_coc[0][2];

    for (int istate=1; istate<nstates; istate++) {

      double dr[3];
      VECTOR_SUB(dr,r_coc[istate],ref);
      //VECTOR_PBC(dr);
      VECTOR_SCALE_ADD(r_cec,dr,C2[istate]);

    }
    VECTOR_ADD(r_cec,r_cec,ref);

    printf("CEC position:              %15.10f %15.10f %15.10f\n",r_cec[0],r_cec[1],r_cec[2]);

  }
}


void Cec::decompose_force(double* force)
{
    
    Atom *atom	     = fmr->atom;
    Matrix *matrix   = fmr->matrix;
    
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;

    double *GSGradient = fmr->matrix->GSGradient;
    double *GSCoeffs   = fmr->matrix->GSCoeffs;
    
    //EVB_Matrix* matrix;
    //if(evb_engine->ncomplex==1) matrix = (EVB_Matrix*)(evb_engine->full_matrix);
    //else matrix = (EVB_Matrix*)(evb_engine->all_matrix[cplx->id-1]);
    
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
                       GSGradient[3*iatom + l] += force[l]*C2*(atom->getCharge(iatom,istate)/qsum_coc[istate]);
                    }      
                }
            } // loop all the atoms in the coc
        } // loop all the coc(s)
    }
    
    /******************************************************************/
    /*** Calculate derivitive of [C(i)^2]  ****************************/
    /******************************************************************/
    
    partial_C_N3(force);

}

void Cec::partial_C_N2(double *force)
{
    /******************************************************************/
    /*** JPCB, 112, 2349 Eq 24-26 *************************************/
    /******************************************************************/
    
    Atom *atom	     = fmr->atom;
    Matrix *matrix   = fmr->matrix;
    
    int natoms       = atom->natoms;
    int nstates      = atom->nstates;
    
    double *GSGradient = matrix->GSGradient;
    double *GSCoeffs   = matrix->GSCoeffs;
    double **Evecs     = matrix->Evecs;
    double *Evals      = matrix->Energy;
    
    int ground = 0; // convention retained from matrix diagonalization
    double factor;


    int nall = atom->nlocal + atom->nghost;
    EVB_Matrix* matrix = evb_matrix;

    int *parent = cplx->parent_id;
    double ***diagonal = matrix->f_diagonal;
    double ***off_diagonal = matrix->f_off_diagonal;
    double ***extra_coupl = matrix->f_extra_coupling;

    // ** Allocate and assign deltaE array ** //
    if (deltaE != NULL) delete [] deltaE;
    deltaE = new double [nstates];
    for (int i=0; i<nstates; ++i) deltaE[i] = E[ground]-E[i];
    
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
    }
    
	// ** Initialize r_coc array to zero ** //
    for (int i=0; i<nstates; ++i)
    {
        r_coc[i][0] = 0.0;
        r_coc[i][1] = 0.0;
        r_coc[i][2] = 0.0;
    }
    
    /********** Eq. 24 ***************/
    
    for(int j=0; j<nstate; j++)
    {
        if(j==ground) continue;
        
        for(int i=0; i<nstate; i++)
        {
            factor = C[i][j]*C[i][ground];
            array1[j][0] += 2*force[0]*(factor*r_coc[i][0]);
            array1[j][1] += 2*force[1]*(factor*r_coc[i][1]);
            array1[j][2] += 2*force[2]*(factor*r_coc[i][2]);
        }
    }
    
    /********** Eq. 25 ***************/
    
    for(int l=0; l<nstate; l++)
    {
        for(int j=0; j<nstate; j++)
        {
            if(j==ground) continue;
            
            factor = C[l][j]/deltaE[j];
            array2[l][0] += factor*array1[j][0];
            array2[l][1] += factor*array1[j][1];
            array2[l][2] += factor*array1[j][2];
        }
        
    }
    
    /*********** Eq. 26 ******************/
    
    for(int l=0; l<nstate; l++)
    {
        double prefactor[3];
        prefactor[0] = C[l][ground]*array2[l][0];
        prefactor[1] = C[l][ground]*array2[l][1];
        prefactor[2] = C[l][ground]*array2[l][2];
        
        for(int i=0; i<nall; i++)
        {
            f[i][0] -= (diagonal[l][i][0])*prefactor[0];
            f[i][1] -= (diagonal[l][i][1])*prefactor[1];
            f[i][2] -= (diagonal[l][i][2])*prefactor[2];
        }
    }
    
    for(int k=1; k<nstate; k++)
    {
        int l = k;
        int m = parent[k];
        
        double prefactor[3];
        prefactor[0] = C[m][ground]*array2[l][0]+C[l][ground]*array2[m][0];
        prefactor[1] = C[m][ground]*array2[l][1]+C[l][ground]*array2[m][1];
        prefactor[2] = C[m][ground]*array2[l][2]+C[l][ground]*array2[m][2];
        
        for(int i=0; i<nall; i++)
        {
            f[i][0] -= off_diagonal[k-1][i][0]*prefactor[0];
            f[i][1] -= off_diagonal[k-1][i][1]*prefactor[1];
            f[i][2] -= off_diagonal[k-1][i][2]*prefactor[2];
        }
    }
    
    for(int k=0; k<cplx->nextra_coupling; k++)
    {
        int l = cplx->extra_i[k];
        int m = cplx->extra_j[k];
        
        double prefactor[3];
        prefactor[0] = C[m][ground]*array2[l][0]+C[l][ground]*array2[m][0];
        prefactor[1] = C[m][ground]*array2[l][1]+C[l][ground]*array2[m][1];
        prefactor[2] = C[m][ground]*array2[l][2]+C[l][ground]*array2[m][2];
        
        for(int i=0; i<nall; i++)
        {
            f[i][0] -= extra_coupl[k][i][0]*prefactor[0];
            f[i][1] -= extra_coupl[k][i][1]*prefactor[1];
            f[i][2] -= extra_coupl[k][i][2]*prefactor[2];
        }
    }
}


