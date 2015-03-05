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
#include <algorithm>
#include <complex>
#include <math.h>

using namespace FMR_NS;

/*-----------------------------------------------------------------
  Constructor
-----------------------------------------------------------------*/
Electrostatics::Electrostatics(FMR *fmr) : Pointers(fmr)
{

    // Basic initializations
    eigax = NULL;
    eigay = NULL;
    eigaz = NULL;
    eigbx = NULL;
    eigby = NULL;
    eigbz = NULL;
    eigcx = NULL;
    eigcy = NULL;
    eigcz = NULL;


}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/
Electrostatics::~Electrostatics()
{
    int atom = fmr->atom->natoms;

    if (eigax != NULL) {
	for (int i=0; i<natoms; ++i) {
	    delete [] eigax[i];
            delete [] eigay[i];
	    delete [] eigaz[i];
            delete [] eigbx[i];
            delete [] eigby[i];
            delete [] eigbz[i];
            delete [] eigcx[i];
            delete [] eigcy[i];
            delete [] eigcz[i];
	}
        delete [] eigax;
        delete [] eigay;
        delete [] eigaz;
        delete [] eigbx;
        delete [] eigby;
        delete [] eigbz;
        delete [] eigcx;
        delete [] eigcy;
        delete [] eigcz;	    
    }

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

void Electrostatics::setup()
{

    Run *run       = fmr->run;  
    Atom *atom     = fmr->atom; 
    int natoms     = atom->natoms;
    int nstates    = fmr->atom->nstates;
                                
    int cellA      = fmr->atom->cellA;
    int cellB      = fmr->atom->cellB;
    int cellC      = fmr->atom->cellC;

    int *nbox_ewald = new int [3];
    for (int l=0; l<3; ++l) nbox_ewald[l] = 0.0;
	
    int *lmax_ewald = new int [3];
    for (int l=0; l<3; ++l) lmax_ewald[l] = 0.0; 
                                
    eps_ewald   = 1.0e-8;
    ratio_ewald = 4.0; 

    double snew = 0.0;
    double sold;
    for (int i=0; i<1000; ++i) {
        sold = snew;
        snew = exp(-snew+1.0)/(eps_ewald + exp(-snew));
	sdif = abs(sold/snew-1.0);
	if (sdif < 1e-15) break;
    }

    double s_ewald = sqrt(snew);   
    double alpha_ewald = ((ratio_ewald+natom*pi)**3/volume**2)**(1.0/6.0);
    double rcut_ewald = s_ewald/alpha_ewald;
    
    absx = cellA;
    absy = cellB;
    absz = cellC; 

    lmax_ewald[0] = (int)(s_ewald*absx*alpha_ewald/M_PI) + 1;
    lmax_ewald[1] = (int)(s_ewald*absy*alpha_ewald/M_PI) + 1;
    lmax_ewald[2] = (int)(s_ewald*absz*alpha_ewald/M_PI) + 1;

    absa = 1.0/(cellA*cellA);           
    absb = 1.0/(cellB*cellB);
    absc = 1.0/(cellC*cellC);

    nbox_ewald[0] = (int)(2.0*rcut_ewald*absa) + 1;
    nbox_ewald[1] = (int)(2.0*rcut_ewald*absb) + 1;
    nbox_ewald[2] = (int)(2.0*rcut_ewald*absc) + 1;

}

void Electrostatics::rs_ewald() 
{

    Run *run       = fmr->run;
    Atom *atom     = fmr->atom;
    int natoms     = atom->natoms;
    int nstates    = atom->nstates;
    int nfragments = atom->nfragments;
    int my_rank    = fmr->my_rank;
    int nf2        = nfragments*nfragments;


    double rcut_ewald2 = rcut_ewald*rcut_ewald;

    for (istate=0; istate<nstates; ++istate) {

        for (int iatom=0; iatom<natoms; ++iatom) {
            for (int jatom=0; jatom<natoms; ++jatom) {

                qi = atom->getCharge(iatom, istate):
 		qj = atom->getCharge(jatom, istate);

                qiqj = qi*qj;

                if (qiqj == 0.0) continue;

                for (int jx=2; jx<nbox_ewald[0]; ++jx) {
		    for (int jy=2; jy<nbox_ewald[1]; ++jy) {
			for (int jz=2; jz<nbox_ewald[2]; ++jz) { 
                	    
			    double j2 = jx*jx + jy*jy + jz*jz;
			
			    double xij = atom->coord[3*iatom] - atom->coord[3*jatom]; 
			    double yij = atom->coord[3*iatom] - atom->coord[3*jatom];
			    double zij = atom->coord[3*iatom] - atom->coord[3*jatom];

			    xij = xij - cellA*jx;
			    yij = yij - cellB*jy;
			    zij = zij - cellC*jz;

			    aij = cellA*nbox_ewald[0]*xij;
			    bij = cellB*nbox_ewald[1]*yij;
			    cij = cellC*nbox_ewald[2]*zij;

			    aij -= round(aij);
			    bij -= round(bij);
			    cij -= round(cij);

			    xij = cellA*nbox_ewald[0]*aij;
			    yij = cellB*nbox_ewald[1]*bij;
			    zij = cellC*nbox_ewald[2]*cij;

			    double r2 = xij*xij + yij*yij + zij*zij;
				
			    if (r2>rcut_ewald2) continue;
	
			    double r = sqrt(r2);
			
			    double rinv  = 1.0/r;
			    double rinv2 = rinv*rinv;
			    double rinv3 = rinv*rinv2;
	
			    double ar = alpha_ewald*r;
			    double invpi = 1.0/sqrt(M_PI);		    
	
			    pot += qiqj*erfc(ar)*rinv

			    double exponent =  -alpha_ewald*alpha_ewald * xij*xij
			    factor = erfc(ar)*rinv3 + 2*alpha_ewald*invpi*exp(exponent);

                            factor = qi*qj*factor;
			
			    fxi = factor*xij;
			    fyi = factor*yij;
			    fzi = factor*zij;

                            fmo_gradients[3*natoms*istate + 3*iatom + 0] += fxi;
                            fmo_gradients[3*natoms*istate + 3*iatom + 1] += fyi;
                            fmo_gradients[3*natoms*istate + 3*iatom + 2] += fzi;
		 
 			}
		    }
		}
	    }
	}
	
        fmo_energies[istate] += 0.5*pot;
    }
/*    
    for (int istate=0; istate<nstate; ++istate) {
        for (int i=0; i<natoms; ++i) {
	    for (int j=0; j<natoms; ++j) {
		
	        double qi = atom->getCharge(iatom, istate):
                double qj = atom->getCharge(jatom, istate);

                double xij = atom->coord[3*iatom] - atom->coord[3*jatom];
                double yij = atom->coord[3*iatom] - atom->coord[3*jatom];
                double zij = atom->coord[3*iatom] - atom->coord[3*jatom];

	        double rij = sqrt(xij*xij + yij*yij + zij*zij);
	
		double rinv = 1.0/rij;
		
		double uij = (factor - 1.0)*qi*qj*rinv;

		fmo_energies[istate] += uij;
		
		double duij = -uij*rinv;
		double fxi  = -duij*xij*rinv;
		double fyi  = -duij*yij*rinv;
		double fzi  = -duij*zij*rinv;
		
		fmo_gradients[3*natoms*istate + 3*i + 0] += fxi;
		fmo_gradients[3*natoms*istate + 3*i + 1] += fyi;
		fmo_gradients[3*natoms*istate + 3*i + 2] += fzi;
	
                fmo_gradients[3*natoms*istate + 3*j + 0] -= fxi;
                fmo_gradients[3*natoms*istate + 3*j + 1] -= fyi;
                fmo_gradients[3*natoms*istate + 3*j + 2] -= fzi;

	    }
	}
    }
*/
}

void Electrostatics::fs_ewald() 
{

    Run *run       = fmr->run;
    Atom *atom     = fmr->atom;
    double *xx     = atom->coord;
    int natoms     = atom->natoms;
    int nstates    = atom->nstates;
    int my_rank    = fmr->my_rank;

    int cellA      = fmr->atom->cellA;
    int cellB      = fmr->atom->cellB;
    int cellC      = fmr->atom->cellC;

    if (eigax != NULL) {
        for (int i=0; i<natoms; ++i) {
            delete [] eigax[i];
            delete [] eigay[i];
            delete [] eigaz[i];
            delete [] eigbx[i];
            delete [] eigby[i];
            delete [] eigbz[i];
            delete [] eigcx[i];
            delete [] eigcy[i];
            delete [] eigcz[i];
        }
        delete [] eigax;
        delete [] eigay;
        delete [] eigaz;
        delete [] eigbx;
        delete [] eigby;
        delete [] eigbz;
        delete [] eigcx;
        delete [] eigcy;
        delete [] eigcz;           
    }

    // Allocate eig based on number of atoms
    eigax = new complex<double>*[natoms];
    eigay = new complex<double>*[natoms];
    eigaz = new complex<double>*[natoms];
    eigbx = new complex<double>*[natoms];
    eigby = new complex<double>*[natoms];
    eigbz = new complex<double>*[natoms];
    eigcx = new complex<double>*[natoms];
    eigcy = new complex<double>*[natoms];
    eigcz = new complex<double>*[natoms];

    int lmax0 = lmax_ewald[0]*2 + 1;
    int lmax1 = lmax_ewald[1]*2 + 1;
    int lmax2 = lmax_ewald[2]*2 + 1;

    for (int i=0; i<natoms; ++i) {
        eigax = new complex<double>[lmax0];
        eigay = new complex<double>[lmax0];
        eigaz = new complex<double>[lmax0];
        eigbx = new complex<double>[lmax1];
        eigby = new complex<double>[lmax1];
        eigbz = new complex<double>[lmax1];
        eigcx = new complex<double>[lmax2];
        eigcy = new complex<double>[lmax2];
        eigcz = new complex<double>[lmax2];	
    }

    for (int istate=0; istate<nstates; ++istate) {

        double ax = 2.0*M_PI*(1.0/cellA);
	double ay = 0.0;
	double az = 0.0;
	
	double bx = 0.0;
	double by = 2.0*M_PI*(1.0/cellB);
	double bz = 0.0;
	
	double cx = 0.0;
	double cy = 0.0;
	double cz = 2.0*M_PI*(1.0/cellC);

	a2 = ax*ax + ay*ay + az*az;
	b2 = bx*bx + by*by + bz*bz;
	c2 = cx*cx + cy*cy + cz*cz;

	al2 = a2*lmax_ewald[0]*lmax_ewald[0];
	bl2 = b2*lmax_ewald[1]*lmax_ewald[1];
	cl2 = c2*lmax_ewald[2]*lmax_ewald[2];

	g2max = std::min(al2, std::min(bl2, cl2));

	for (int k=0; k<natoms; ++k) {

	    eigax[k][lmax0] = complex<double>(1.0, 0.0);
            eigax[k][lmax0] = complex<double>(1.0, 0.0); 
            eigax[k][lmax0] = complex<double>(1.0, 0.0); 
            eigbx[k][lmax1] = complex<double>(1.0, 0.0); 
            eigbx[k][lmax1] = complex<double>(1.0, 0.0); 
            eigbx[k][lmax1] = complex<double>(1.0, 0.0); 
            eigcx[k][lmax2] = complex<double>(1.0, 0.0); 
            eigcx[k][lmax2] = complex<double>(1.0, 0.0); 
            eigcx[k][lmax2] = complex<double>(1.0, 0.0);

            eigax[k][lmax0+1] = complex<double>(cos(ax*xx[3*k+0]), sin(ax*xx[3*k+0]));
            eigax[k][lmax0+1] = complex<double>(cos(ax*xx[3*k+1]), sin(ax*xx[3*k+1]));
            eigax[k][lmax0+1] = complex<double>(cos(ax*xx[3*k+2]), sin(ax*xx[3*k+2]));
            eigax[k][lmax1+1] = complex<double>(cos(ax*xx[3*k+0]), sin(ax*xx[3*k+0]));
            eigax[k][lmax1+1] = complex<double>(cos(ax*xx[3*k+1]), sin(ax*xx[3*k+1]));
            eigax[k][lmax1+1] = complex<double>(cos(ax*xx[3*k+2]), sin(ax*xx[3*k+2]));
            eigax[k][lmax2+1] = complex<double>(cos(ax*xx[3*k+0]), sin(ax*xx[3*k+0]));
            eigax[k][lmax2+1] = complex<double>(cos(ax*xx[3*k+1]), sin(ax*xx[3*k+1]));
            eigax[k][lmax2+1] = complex<double>(cos(ax*xx[3*k+2]), sin(ax*xx[3*k+2]));

            eigax[k][lmax0-1] = complex<double>(cos(ax*xx[3*k+0]), -sin(ax*xx[3*k+0]));
            eigax[k][lmax0-1] = complex<double>(cos(ax*xx[3*k+1]), -sin(ax*xx[3*k+1]));
            eigax[k][lmax0-1] = complex<double>(cos(ax*xx[3*k+2]), -sin(ax*xx[3*k+2]));
            eigax[k][lmax1-1] = complex<double>(cos(ax*xx[3*k+0]), -sin(ax*xx[3*k+0]));
            eigax[k][lmax1-1] = complex<double>(cos(ax*xx[3*k+1]), -sin(ax*xx[3*k+1]));
            eigax[k][lmax1-1] = complex<double>(cos(ax*xx[3*k+2]), -sin(ax*xx[3*k+2]));
            eigax[k][lmax2-1] = complex<double>(cos(ax*xx[3*k+0]), -sin(ax*xx[3*k+0]));
            eigax[k][lmax2-1] = complex<double>(cos(ax*xx[3*k+1]), -sin(ax*xx[3*k+1]));
            eigax[k][lmax2-1] = complex<double>(cos(ax*xx[3*k+2]), -sin(ax*xx[3*k+2]));

	    for (int l=2; l<=lmax_ewald[0]; ++l) {
		idx0 = lmax_ewald[0] + l;
                idx1 = lmax_ewald[0] - l;

		eigax[k][idx0] = eigax[k][idx0-1]*eigax[k][lmax_ewald[0]+1];
		eigay[k][idx0] = eigay[k][idx0-1]*eigay[k][lmax_ewald[0]+1];
		eigaz[k][idx0] = eigaz[k][idx0-1]*eigaz[k][lmax_ewald[0]+1];
		eigax[k][idx1] = std::conj(eigax[k][idx0]);
		eigay[k][idx1] = std::conj(eigay[k][idx0]);
		eigaz[k][idx1] = std::conj(eigaz[k][idx0]);
	    }
            for (int l=2; l<=lmax_ewald[1]; ++l) {
                idx0 = lmax_ewald[1] + l;
                idx1 = lmax_ewald[1] - l;
                
                eigbx[k][idx0] = eigbx[k][idx0-1]*eigbx[k][lmax_ewald[1]+1];
                eigby[k][idx0] = eigby[k][idx0-1]*eigby[k][lmax_ewald[1]+1];
                eigbz[k][idx0] = eigbz[k][idx0-1]*eigbz[k][lmax_ewald[1]+1];
                eigbx[k][idx1] = std::conj(eigbx[k][idx0]);
                eigby[k][idx1] = std::conj(eigby[k][idx0]);
                eigbz[k][idx1] = std::conj(eigbz[k][idx0]);
            }
            for (int l=2; l<=lmax_ewald[2]; ++l) {
                idx0 = lmax_ewald[2] + l;
                idx1 = lmax_ewald[2] - l;
                
                eigcx[k][idx0] = eigcx[k][idx0-1]*eigcx[k][lmax_ewald[2]+1];
                eigcy[k][idx0] = eigcy[k][idx0-1]*eigcy[k][lmax_ewald[2]+1];
                eigcz[k][idx0] = eigcz[k][idx0-1]*eigcz[k][lmax_ewald[2]+1];
                eigcx[k][idx1] = std::conj(eigcx[k][idx0]);
                eigcy[k][idx1] = std::conj(eigcy[k][idx0]);
                eigcz[k][idx1] = std::conj(eigcz[k][idx0]);
            }	
	}	    		

	double volume = cellA*cellB*cellC;

        double factor_1 = (4.0*M_PI)/(2.0*volume);

        for (int la=0; la<lmax_ewald[0]; ++la) {
            for (int lb=0; lb<lmax_ewald[1]; ++lb) {
                for (int lc=0; lc<lmax_ewald[2]; ++lc) {
		
		    double l2 = la*la + lb*lb + lc*lc;
 		
		    if (l2 == 0) continue;
		
		    if (la == 0) factor_2 = 1.0;
		    else factor_2 = 2.0;
		
		    double gx = ax*la + bx*lb + cx*lc;
		    double gy = ay*la + by*lb + cy*lc;
		    double gz = az*la + bz*lb + cz*lc;

		    doubl g2 = gx*gx + gy*gy + gz*gz; 

		    if (g2 > g2max) continue;
		
		    double factor_3 = exp(-g2/(4.0*alpha_ewald*alpha_ewald))/g2;

		    double qcos = 0.0;
		    double qsin = 0.0;

		    for (int k=0; k<natoms; ++k) {
		    
		        double cos_gxyz = std::real((eigax[k][la]*eigbx[k][lb]*eigcx[k][lc]*
		   				     eigay[k][la]*eigby[k][lb]*eigcy[k][lc]*
			    			     eigaz[k][la]*eigbz[k][lb]*eigcz[k][lc]));
		
 		        double sin_gxyz = std::imag((eigax[k][la]*eigbx[k][lb]*eigcx[k][lc]*
                                                     eigay[k][la]*eigby[k][lb]*eigcy[k][lc]*
                                                     eigaz[k][la]*eigbz[k][lb]*eigcz[k][lc]));

		        qcos += qk*cos_gxyz;
		        qsin += qk*sin_gxyz;
		    }
		
		    double qexp2 = qcos*qcos + qsin*qsin;
		
		    fmo_energies[istate] += factor_1*factor_2*factor_3*qexp2;
	
  		    for (int k=0; k<natoms; ++k) {
		    
                        double cos_gxyz = std::real((eigax[k][la]*eigbx[k][lb]*eigcx[k][lc]*
                                                     eigay[k][la]*eigby[k][lb]*eigcy[k][lc]*
                                                     eigaz[k][la]*eigbz[k][lb]*eigcz[k][lc]));

                        double sin_gxyz = std::imag((eigax[k][la]*eigbx[k][lb]*eigcx[k][lc]*
                                                     eigay[k][la]*eigby[k][lb]*eigcy[k][lc]*
                                                     eigaz[k][la]*eigbz[k][lb]*eigcz[k][lc]));
		
		        double factor_4 = sin_gxyz*qcos - cos_gxyz*qsin;
		        double factor_5 = 2.0*qk*factor_1*factor_2*factor_3*factor4;

		        double fxi = factor_5*gx;
		        double fyi = factor_5*gy;
		        double fzi = factor_5*gz;

		        fmo_gradients[3*natoms*istate + 3*k + 0] += fxi;
		        fmo_gradients[3*natoms*istate + 3*k + 1] += fyi;
		        fmo_gradients[3*natoms*istate + 3*k + 2] += fzi;
		    }
	        }
	    }
        }
    }
}

void Electrostatics::charge_ewald()
{

    Run *run       = fmr->run;
    Atom *atom     = fmr->atom;
    double *xx     = atom->coord;
    int natoms     = atom->natoms;
    int nstates    = atom->nstates;
    int my_rank    = fmr->my_rank;

    int cellA      = fmr->atom->cellA;
    int cellB      = fmr->atom->cellB;
    int cellC      = fmr->atom->cellC;

    double volume = cellA*cellB*cellC;
    double factor = M_PI/(2.0*volume*alpha_ewald*alpha_ewald);


    for (int istate=0; istate<nstates; ++istate) {

	double qsum = 0.0;
	for (int k=0; k<natoms; ++k) {
	    
	    double qk = atom->getCharge(k, istate);  
	    qsum += qk;

	}
	fmo_energies[istate] += qsum*qsum*factor;

    }
} 

void Electrostatics::self_ewald()
{

    Run *run       = fmr->run;
    Atom *atom     = fmr->atom;
    double *xx     = atom->coord;
    int natoms     = atom->natoms;
    int nstates    = atom->nstates;
    int my_rank    = fmr->my_rank;

    int cellA      = fmr->atom->cellA;
    int cellB      = fmr->atom->cellB;
    int cellC      = fmr->atom->cellC;

    double factor = alpha_ewald/sqrt(M_PI);
    
    for (int istate=0; istate<nstates; ++istate) {

	double q2sum = 0.0;
	for (int k=0; k<natoms; k++) {
	
	    double qk = atom->getCharge(k, istate);
            q2sum += qk;
	}
	fmo_energies[istate] += q2sum;
    }
}


void Electrostatics::ewald()
{



}

