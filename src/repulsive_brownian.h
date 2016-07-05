
#ifndef REPULSIVE_BROWNIAN_H
#define REPULSIVE_BROWNIAN_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"

#include "omp.h"

namespace REPULSIVE_BROWNIAN
{
  struct TEMPORAL_VARIABLE; // tell compiler to expect struct def.
  
  MKL_LONG main_EQUILIBRATION(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition);

  double record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, MATRIX& energy, const MKL_LONG index_t_now);
  double report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, TEMPORAL_VARIABLE& VAR);
  double OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, RDIST& R_boost, MKL_LONG* tmp_index_vec, const MKL_LONG N_THREADS_BD);
  double OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, POTENTIAL_SET& POTs, RDIST& R_boost, MATRIX* vec_boost_Nd_parallel, MATRIX* force_repulsion, MATRIX* force_random, RNG_BOOST& RNG, const MKL_LONG N_THREADS_BD, COND& given_condition, TEMPORAL_VARIABLE& VAR);
}

struct REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE
{
  // the variables on here is supporting for interface
  // do not use the core variable inside this structure
  // for designing purpose, it is of importance to use individual variable definition outside of this structure
  MKL_LONG N_THREADS_BD;
  MKL_LONG *tmp_index_vec;
  MATRIX *vec_boost_Nd_parallel;
  MATRIX *force_repulsion;
  MATRIX *force_random;
  MKL_LONG N_skip;
  MKL_LONG Nt;
  MKL_LONG Np;
  MKL_LONG N_basic;
  double time_LV, time_DIST, time_file, time_AN, time_RECORDED;
  double time_LV_init, time_LV_force, time_LV_update;
  double simulation_time;

  // related with shear flow
  bool SIMPLE_SHEAR;
  double Wi_tau_R;
  MKL_LONG shear_axis;
  MKL_LONG shear_grad_axis;
  double shear_PBC_shift;
  
  bool INITIALIZATION;
  // member functions
  /* double construct(TEMPORAL_VARIABLE& VAR, COND& given_condition, MKL_LONG given_N_basic); */
  /* double destruct(TEMPORAL_VARIABLES& VAR); */
  TEMPORAL_VARIABLE()
  {
    std::cout << "ERR: Basic constructor for TEMPORAL_VARIABLE structure is not supported\n";
  }
  TEMPORAL_VARIABLE(COND& given_condition, MKL_LONG given_N_basic);
  virtual ~TEMPORAL_VARIABLE();

};

#endif
