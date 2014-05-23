/* AWGL */
#ifndef FMR_RUN_H
#define FMR_RUN_H

#include "pointers.h"

// Run type definitions
#define RUN_ENERGY   0
#define RUN_FORCE    1
#define RUN_MOLDYN   2
#define RUN_TYPE_MIN 0
#define RUN_TYPE_MAX 2

namespace FMR_NS {


class Run : protected Pointers {
 public:

   // ** Constructor/destructor ** //
   Run(FMR *);
   ~Run();

   // ** Variables ** //
   int run_type;              // Run type index
   int FMO_only;              // Run with only FMO, no Multi-state FMO
   int EnvApprox;             // Environment approximation on/off

   double cut_dimer;          // radius of dimer cut

   char scratch_dir[256];     // Full path of Q-Chem scratch directory 
   char exec[256];            // Full path of Q-Chem exectuable

   char gamess_version[256];  // Gamess version number
   int  gamess_ncores;        // Number of cores to run FMO with GDDI

   char correlation[256];     // correlation method
   char exchange[256];	      // exchange method
   char basis[256];	      // basis set name
   char algorithm[256];       // scf algorithm

   int n_monomers;            // # monomers within range
   int n_monomers_tmp;        // (# fragments) * (# states)
   int n_dimers;              // # fragments within range
   int n_dimers_tmp;          // ((# fragments) * (# fragments - 1) / 2 ) * (# states)
   int n_dimers_sq;           // ((# fragments) * (# fragments) ) * (# states)
   int    *monomer_queue;     // list of monomers to calculate
   int    *dimer_queue;       // list of dimers to calculate
   double *fmo_energies;      // FMO energies for each state
   double *monomer_energies;  // FMO monomer energies
   double *dimer_energies;    // FMO dimer energies
   double *fmo_gradients;     // FMO gradients for each state
   double *monomer_gradients; // FMO monomer gradients
   double *dimer_gradients;   // FMO dimer gradients

   int *monomer_list;
   int *dimer_list;

   // ** Functions ** //
   void run_calculation();     // General run calculation
   void setup();               // Setup calculation ...
   void post_force();
   void calculate_energy();    // Computes energy only
   void calculate_force();     // Computes energy + force
   void calculate_moldyn();    // Computes molecular dynamics
   void do_qchem_calculations(int); // Performs the all Q-Chem FMO calculations in parallel
   void do_qchem_calculations_env();  // Performs the FMO calculations in parallel, using env approximation
   void do_nwchem_calculations(int);  // Performs all the NWChem calculations in parallel
   void do_nwchem_calculations_cutoff(int);
   void do_nwchem_calculations_cutoff_old(int);
   void do_nwchem_calculations_env(); // Performs the FMO calculations in parallel, using env approximation
   void do_gamess_calculations(int);  // Performs all the Gamess calculations in parallel
   void do_gamess_calculations_env(); // Performs the FMO calculations in parallel, using env approximation
   void perturb_coords();      // Just to perturb the coordinates slightly randomly

};

}

#endif
