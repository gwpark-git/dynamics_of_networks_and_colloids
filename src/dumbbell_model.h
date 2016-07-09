
#ifndef DUMBBELL_MODEL_H
#define DUMBBELL_MODEL_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"
#include "brownian.h"
#include "omp.h"

namespace DUMBBELL
{
  struct DUMBBELL_VARIABLE; // tell compiler to expect struct def.
  
  MKL_LONG main_DUMBBELL(TRAJECTORY& TRAJ, CONNECTIVITY& CONNECT, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition);
  MKL_LONG generate_dumbbell_connectivity(CONNECTIVITY& CONNECT);

  double record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, MATRIX& energy, const MKL_LONG index_t_now);
  double report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, DUMBBELL_VARIABLE& VAR);
  /* double OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, RDIST& R_boost, MKL_LONG* tmp_index_vec, const MKL_LONG N_THREADS_BD); */
  double OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, CONNECTIVITY& CONNECT, POTENTIAL_SET& POTs, MATRIX* force_random, MATRIX* force_spring, RNG_BOOST& RNG, RDIST& R_boost, const MKL_LONG N_THREADS_BD, COND& given_condition, DUMBBELL_VARIABLE& VAR);
}

struct DUMBBELL::DUMBBELL_VARIABLE : BROWNIAN::BROWNIAN_VARIABLE
{
  // the variables on here is supporting for interface
  // do not use the core variable inside this structure
  // for designing purpose, it is of importance to use individual variable definition outside of this structure

  // additional member variables
  double time_LV_force_connector;
  MATRIX *force_spring;
  
  DUMBBELL_VARIABLE()
  {
    std::cout << "ERR: Basic constructor for BROWNIAN_VARIABLE structure is not supported\n";
  }
 DUMBBELL_VARIABLE(COND& given_condition, MKL_LONG given_N_basic) : BROWNIAN::BROWNIAN_VARIABLE(given_condition, given_N_basic)
    {
      time_LV_force_connector = 0.;
      force_spring = new MATRIX [Np];

      MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
      for(MKL_LONG i=0; i<Np; i++)
	{
	  force_spring[i].initial(N_dimension, 1, 0.);
	}
      /* all the necessary information stored in the construction for BROWNIAN::BROWNIAN_VARIABLE */
    }
  
  virtual ~DUMBBELL_VARIABLE()
    {
      if(INITIALIZATION)
	delete[] force_spring;
    }
};

#endif
