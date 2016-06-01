
#include <iostream>
#include "../lib/trajectory.h"              // TRAJECTORY class
#include "../lib/association.h"             // ASSOCIATION class
#include "../lib/handle_association.h"      // CHAIN class
#include "../lib/potential.h"               // POTENTIAL_SET class
#include "../lib/file_IO.h"                 // RECORD_DATA class, COND class

#include <string>

#include "repulsive_brownian.h"             // EQUILIBRATION simulation
#include "stochastic_HEUR_flowers.h"        // stochastic simulation for HEUR flowers
using namespace std;

int help()
{
  cout << "USAGE of Brownian Dynamics with Quadratic Repulsive Potential\n";
  cout << "argv[1] == condition files including variables:\n";
  cout << "\tLONG INTEGER format: N_dimension, Nt, Np\n";
  cout << "\tDOUBLE FLOAT format: dt, repulsion_coefficient, effective_distance\n";
  cout << "\tSTRING format: filename_trajectory, filename_energy, filename_energy_detail\n";
  return 0;
}


int main(int argc, char* argv[])
{
  if(argc==1)
    {
      help();
      return 0;
    }
  else
    {
      COND given_condition(argv[1]);
      // if(atoi(given_condition("N_THREADS_BD").c_str())==1)
      //   {
      //     std::cout << "ERR: N_THREADS_BD = 1 has bug at this moment\n";
      //     return -1;
      //   }
      given_condition.cond_print();

      MKL_LONG N_basic = 2;
      if (given_condition("Integrator") == "Euler") 
        {
          // at the moment, the only method available is simple Euler integrator
          // for this reason, this condition will not be varied during condition test
          // however, the explicit description should be needed for condition data
          // that is because to prevent potential errors when we implement higher order methods.
          // Note that the only needed information for Euler is the previous time step. Therefore, we need previous and now which is in 2.
          N_basic = 2;

        }
      TRAJECTORY TRAJ(given_condition, N_basic);

      POTENTIAL_SET POTs;
      RECORD_DATA DATA(given_condition);
      if(given_condition("Method") == "NAPLE_ASSOCIATION")
        {
          if(given_condition("Step") == "EQUILIBRATION")
            {
              FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
              REPULSIVE_BROWNIAN::main_EQUILIBRATION(TRAJ, POTs, DATA, given_condition);
            }
          else
            {
              ASSOCIATION CONNECT(given_condition);
              FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTs, given_condition);
              CHAIN_HANDLE CHAIN(given_condition, CONNECT);
              // if(given_condition("CONTINUATION_CONNECTION") == "TRUE")
              //   {
              //     if(given_condition("tracking_individual_chain") == "TRUE")
              //       {
              //         if(given_condition("CONTINUATION_CHAIN")=="TRUE")
              //           {
              //             printf("CONTINUATION the chain information is not tested for 'tracking individual chain' mode\n");
              //             return -1;
              //           }
              //       }
              //     // CHAIN.allocate_existing_bridges(CONNECT);
              //   }
              stochastic_simulation_HEUR_flowers(TRAJ, POTs, CONNECT, CHAIN, DATA, given_condition);
            }
        }
      else
        {
          printf("Method except NAPLE_ASSOCIATION is not clearly defined. The given Method input is %s", given_condition("Method").c_str());
        }
    }
  return 0;
}



