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

Atom::Atom(FMR *fmr) : Pointers(fmr)
{
  // Basic initializations
  natoms       		= 0;
  nfragments   		= 0;
  nstates      		= 0;
  prev_nstates 		= 0;
  ireactive    		= 0;
  coord        		= NULL;
  force        		= NULL;
  veloc        		= NULL;
  mass         		= NULL;
  symbol       		= NULL;
  name			= NULL;
  fragment     		= NULL;
  reactive     		= NULL;
  available    		= NULL;
  hop          		= NULL;
  environment  		= NULL;
  totalMass    		= 1.0;

  // periodic boundary conditions	
  cellA 		= 1000.0;  
  cellB			= 1000.0;
  cellC			= 1000.0;
  na			= 0;
  nb			= 0;
  nc 			= 0;
  afield		= 0;
  bfield		= 0;  
  cfield		= 0;

  // MM charges
  qO_SPCE = -0.8476;
  qH_SPCE = -(qO_SPCE) * 0.5;
  qO_hydronium = -0.5;
  qH_hydronium = (1.0 - qO_hydronium) / 3.0;
  qCl = -1;

  

}

/*-----------------------------------------------------------------
  Destructor
-----------------------------------------------------------------*/

Atom::~Atom()
{
  if (coord     != NULL) delete [] coord; 
  if (force     != NULL) delete [] force; 
  if (veloc     != NULL) delete [] veloc; 
  if (mass      != NULL) delete [] mass; 
  if (symbol    != NULL) delete [] symbol;
  if (name      != NULL) delete [] name; 
  if (fragment  != NULL) delete [] fragment; 
  if (reactive  != NULL) delete [] reactive; 
  if (available != NULL) delete [] available; 
  if (hop       != NULL) delete [] hop; 
}

void Atom::setGlobalBox() 
{
  xprd = cellA;
  yprd = cellB;
  zprd = cellC;

//  h[0] = xprd;
//  h[1] = yprd;
//  h[2] = zprd;
//  h_inv[0] = 1.0/h[0];
//  h_inv[1] = 1.0/h[1];
//  h_inv[2] = 1.0/h[2];

  xprd_half = 0.5*xprd;
  yprd_half = 0.5*yprd;
  zprd_half = 0.5*zprd;
}

void Atom::minimum_image(double &dx, double &dy, double &dz) {
  if (na) {
    if (fabs(dx) > xprd_half) {
      if (dx < 0.0) dx += xprd;
      else dx -= xprd;
    }
  }
  if (nb) {
    if (fabs(dy) > yprd_half) {
      if (dy < 0.0) dy += yprd;
      else dy -= yprd;
    }
  }
  if (nc) {
    if (fabs(dz) > zprd_half) {
      if (dz < 0.0) dz += zprd;
      else dz -= zprd;
    }
  }
}

/*-----------------------------------------------------------------
  Set all the atom masses to a.u. 
-----------------------------------------------------------------*/
void Atom::setAtomMasses()
{
  // Also compute the total mass
  totalMass = 0.0;
  for (int i=0; i<natoms; ++i) {
    mass[i] = getAtomMass(i) * fmr->math->amu2au;
    totalMass += mass[i];
  }
}
