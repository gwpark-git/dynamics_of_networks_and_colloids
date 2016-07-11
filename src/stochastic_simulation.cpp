
#include <iostream>
#include "../lib/trajectory.h"		// TRAJECTORY class
#include "../lib/association.h"		// ASSOCIATION class
#include "../lib/handle_association.h"	// CHAIN class
#include "../lib/potential.h"		// POTENTIAL_SET class
#include "../lib/file_IO.h"		// RECORD_DATA class, COND class
#include "../lib/analysis.h" // related with analysis and post-processing


#include <string>

#include "brownian.h"			// basic Brownian motion
#include "repulsive_brownian.h"		// EQUILIBRATION simulation
#include "stochastic_HEUR_flowers.h"	// stochastic simulation for HEUR flowers
#include "dumbbell_model.h"		// dumbbell model based on Brownian motion
using namespace std;

int
help()
{
  cout << "USAGE of Brownian Dynamics with Quadratic Repulsive Potential\n";
  cout << "argv[1] == condition files including variables:\n";
  cout << "\tLONG INTEGER format: N_dimension, Nt, Np\n";
  cout << "\tDOUBLE FLOAT format: dt, repulsion_coefficient, effective_distance\n";
  cout << "\tSTRING format: filename_trajectory, filename_energy, filename_energy_detail\n";
  return 0;
}


int
main
(int argc, char* argv[])
{
  if(argc==1)
    {
      help();
      return 0;
    }
  else if(argc > 2)
    {
      MATRIX ener_data;
      ener_data.read_exist_data(argv[2]);
      printf("READ_FILE_DONE\n");
      MATRIX result(ener_data.rows/2, 7, 0.);
      for(MKL_LONG t=0; t<result.rows; t++)
        result(t, 0) = ener_data(t, 0);
      MKL_LONG N_THREADS = 6;

      using namespace POST_PROCESSING;
      compute_autocorrelation_OMP(ener_data, 21, N_THREADS, result, 1);
      printf("PROCESS_21\n");
      compute_autocorrelation_OMP(ener_data, 22, N_THREADS, result, 2);
      printf("PROCESS_22\n");      
      compute_autocorrelation_OMP(ener_data, 23, N_THREADS, result, 3);
      printf("PROCESS_23\n");      
      compute_autocorrelation_OMP(ener_data, 27, N_THREADS, result, 4);
      printf("PROCESS_27\n");
      compute_autocorrelation_OMP(ener_data, 28, N_THREADS, result, 5);
      printf("PROCESS_28\n");
      compute_autocorrelation_OMP(ener_data, 29, N_THREADS, result, 6);
      printf("PROCESS_29\n");      

      ofstream FILE;
      FILE.open(argv[3]);
      result.fprint(FILE);
      FILE.close();
      return 0;
      // for(MKL_LONG i=0; i<data.rows; i++)
      //   {
      //     for(MKL_LONG j=0; j<data.cols; j++)
      //       {
      //         printf("%3.2lf\t", data(i,j));
      //       }
      //     printf("\n");
      //   }
      
    }
  else
    {
        
      // inp_cnt started with index 1 since 0 is itself
      // single input file case, argc==2, argv[0] == execution, argv[1] == inputfile
      for(MKL_LONG inp_cnt = 1; inp_cnt < argc; inp_cnt ++)
        {
          /*
            This for-loop is deisgn to perform sequence of jobs using the given input file.
            The sequence is exactly the same with the given index for arguments.
            Note that the all the controls should be properly set with the given input files, which means the inheritance from the previous computation in the given input file is exactly matched with the output file name from the previous input file.
           */
          
          COND given_condition(argv[inp_cnt]);
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
              // note that the time scale for Langevin dynamics on here is affected by repulsive coefficient, C_rep, which have differ non-dimensional form compared with pure Brownian motion
              
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

		  HEUR::stochastic_simulation_HEUR_flowers(TRAJ, POTs, CONNECT, CHAIN, DATA, given_condition);
                    // }
                }
            }
          else if(given_condition("Method") == "REPULSIVE_BROWNIAN")
            {
              FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTs, given_condition);
              REPULSIVE_BROWNIAN::main_EQUILIBRATION(TRAJ, POTs, DATA, given_condition);
            }
          else if(given_condition("Method") == "BROWNIAN")
            {
              
              // As mentioned in repulsive Brownian motion, here the time scale is used as the basic Brownian motion. Therefore, the conditional file has time step dt/tauB instead of dt/tauR.
              FORCE::BROWNIAN::MAP_potential_set(POTs, given_condition); // null allocator
              BROWNIAN::main_PURE_BROWNIAN(TRAJ, POTs, DATA, given_condition);
            }
          else if(given_condition("Method") == "DUMBBELL")
            {
              // CONNECTIVITY CONNECT(given_condition); // note that
              MKL_LONG N_max_connection = 1; // note that dumbbell only have one permernant connection between a pir of micelle
	      printf("test:init\n");
              CONNECTIVITY CONNECT(TRAJ.Np, N_max_connection);
	      printf("test:aft:connect\n");
	      DUMBBELL::generate_dumbbell_connectivity(CONNECT); // dumbbell
	      printf("test:aft:dumbbell_connect\n");
              FORCE::DUMBBELL::MAP_potential_set(POTs, given_condition);
	      printf("test:aft:potential_set\n");
              DUMBBELL::main_DUMBBELL(TRAJ, CONNECT, POTs, DATA, given_condition);
            }
          else
            {
              printf("Method except NAPLE_ASSOCIATION is not clearly defined. The given Method input is %s", given_condition("Method").c_str());
            }
        }
    }
  return 0;
}



