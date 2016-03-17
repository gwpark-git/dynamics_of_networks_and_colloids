
#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

#include "trajectory.h"
#include "geometry.h"
#include "potential.h"
#include "fstream"
#include "association.h"
#include <mkl.h>

namespace INTEGRATOR
{
  namespace EULER
  {
    MKL_LONG simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
    MKL_LONG cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index);
    MKL_LONG cal_repulsion_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG index_i);
    MKL_LONG cal_repulsion_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG index_i, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost);
    // note that the vec_boost_Nd_parallel in main function has MATRIX* form
    // In the parallel region, the given variable is transfered by vec_boost_Nd_parallel[i] that is MATRIX
    // Therefore, here should have MATRIX& for call-by-reference rather than MATRIX*.
    MKL_LONG cal_random_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t);
    MKL_LONG cal_random_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, gsl_rng* r_boost);
    
  }
  namespace EULER_ASSOCIATION
  {
    MKL_LONG cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index);
    MKL_LONG cal_connector_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost);
  }

}

namespace ANALYSIS
{
  MKL_LONG GET_dCDF_POTENTIAL_target(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG& index_particle, MKL_LONG& index_target, double& INDEX_dCDF_U_ij, double& dCDF_U_ij, MATRIX& vec_boost_ordered_pdf_ij);

  MKL_LONG CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& mat_energy, MKL_LONG index_t);
  
  double cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
  
  double cal_kinetic_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
  double cal_total_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
  MKL_LONG cal_detail_repulsion(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, const char* fn, MKL_LONG index_t);


  MKL_LONG GET_dCDF_POTENTIAL(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG index_particle, MATRIX& INDEX_dCDF_U, MATRIX& dCDF_U, MATRIX& R_minimum_distance_boost_particle);
  MKL_LONG GET_ORDERED_dCDF_POTENTIAL(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG index_particle, MATRIX& INDEX_dCDF_U, MATRIX& dCDF_U, MATRIX& vec_boost_ordered_pdf, double& time_dCDF, double& time_SORT);
  namespace ANAL_ASSOCIATION
  {
    MKL_LONG CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, MKL_LONG index_t, MATRIX& tmp_vec);
    
    double cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t);
    double cal_potential_energy_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t, MATRIX& tmp_vec);
    

  }
  
}

#endif
