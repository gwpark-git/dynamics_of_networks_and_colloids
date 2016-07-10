
#ifndef BROWNIAN_H
#define BROWNIAN_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"

#include "omp.h"

namespace BROWNIAN
{
  struct BROWNIAN_VARIABLE; // tell compiler to expect struct def.

  MKL_LONG
  main_PURE_BROWNIAN(TRAJECTORY& TRAJ,
		     POTENTIAL_SET& POTs,
		     RECORD_DATA& DATA,
		     COND& given_condition);
  double
  report_simulation_info(TRAJECTORY& TRAJ,
			 MATRIX& energy,
			 BROWNIAN_VARIABLE& VAR);
  
  double
  OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next,
			   POTENTIAL_SET& POTs, MATRIX* force_random,
			   RNG_BOOST& RNG,
			   const MKL_LONG N_THREADS_BD,
			   COND& given_condition, BROWNIAN_VARIABLE& VAR);

  // inline functions
  double
  sum_virial_components(MATRIX& energy);
}

struct
BROWNIAN::
  BROWNIAN_VARIABLE
{
  // the variables on here is supporting for interface
  // do not use the core variable inside this structure
  // for designing purpose, it is of importance to use individual variable definition outside of this structure
  MKL_LONG N_THREADS_BD;

  double volume_PBC_box;
  
  MKL_LONG *tmp_index_vec;
  MATRIX *force_random;
  MKL_LONG N_skip_ener, N_skip_file;
  MKL_LONG Nt;
  MKL_LONG Np;
  MKL_LONG N_basic;
  double time_LV, time_DIST, time_file, time_AN, time_RECORDED;
  double time_LV_init, time_LV_force, time_LV_update;
  double time_LV_force_random;
  double simulation_time;

  // the following are related with virial stress
  // the form of single variable is because of reduction in OpenMP is not allowed array-type
  // hence, the form cannot be simplified
  double RF_random_xx, RF_random_yy, RF_random_zz;
  double RF_random_xy, RF_random_xz, RF_random_yz;
  MKL_LONG N_components_energy;

  // member functions for virials
  MKL_LONG
  virial_initial();
  double
  record_virial_into_energy_array(MATRIX& energy);
  
  // related with shear flow
  bool SIMPLE_SHEAR;
  double Wi_tau_B;
  MKL_LONG shear_axis;
  MKL_LONG shear_grad_axis;
  double shear_PBC_shift;
  bool INITIALIZATION;


  
  // member functions
  BROWNIAN_VARIABLE()
  {
    std::cout << "ERR: Basic constructor for BROWNIAN_VARIABLE structure is not supported\n";
  }
  BROWNIAN_VARIABLE(COND& given_condition,
		    MKL_LONG given_N_basic);
  virtual
  ~BROWNIAN_VARIABLE();
};

// for inline functions
inline MKL_LONG
BROWNIAN::BROWNIAN_VARIABLE::
virial_initial()
{
  RF_random_xx = 0.; RF_random_yy = 0.; RF_random_zz = 0.;
  RF_random_xy = 0.; RF_random_xz = 0.; RF_random_yz = 0.;
  return 0;
}

inline double
BROWNIAN::BROWNIAN_VARIABLE::
record_virial_into_energy_array(MATRIX& energy)
{
  double time_st = dsecnd();
  energy(12) = RF_random_xx/(2.*volume_PBC_box);
  energy(13) = RF_random_yy/(2.*volume_PBC_box);
  energy(14) = RF_random_zz/(2.*volume_PBC_box);

  energy(15) = RF_random_xy/(2.*volume_PBC_box);
  energy(16) = RF_random_xz/(2.*volume_PBC_box);
  energy(17) = RF_random_yz/(2.*volume_PBC_box);
  return time_st - dsecnd();
}

inline double
BROWNIAN::
sum_virial_components(MATRIX& energy)
{
  double time_st = dsecnd();
  // it is not importance for brownian motion
  // but, will be useful when there are other potential contributions
  MKL_LONG index_st = 6;
  MKL_LONG number_of_components = 6;
  for(MKL_LONG i=6; i<12; i++)
    energy(i) = energy(i + index_st);

  return dsecnd() - time_st;
}


#endif
