
#ifndef LIB_EVOLUTION_H
#define LIB_EVOLUTION_H
#include "lib_traj.h"
#include "lib_geometry.h"
#include "lib_potential.h"
#include "fstream"
#include "lib_association.h"
#include <mkl.h>

namespace INTEGRATOR
{
  namespace EULER
  {
    long simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t);
    long cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t, long given_index);
    long cal_repulsion_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t, long index_i);
    long cal_repulsion_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t, long index_i, MATRIX& vec_boost_Nd);
    // note that the vec_boost_Nd_parallel in main function has MATRIX* form
    // In the parallel region, the given variable is transfered by vec_boost_Nd_parallel[i] that is MATRIX
    // Therefore, here should have MATRIX& for call-by-reference rather than MATRIX*.
    long cal_random_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, long index_t);
  }
  namespace EULER_ASSOCIATION
  {
    long cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, long index_t, long given_index);
    long cal_connector_force_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, long index_t, long given_index, MATRIX& vec_boost_Nd);
  }

}

namespace ANALYSIS
{
  long GET_dCDF_POTENTIAL_target(TRAJECTORY& TRAJ, long index_t, POTENTIAL_SET& POTs, long& index_particle, long& index_target, double& INDEX_dCDF_U_ij, double& dCDF_U_ij, MATRIX& vec_boost_ordered_pdf_ij);

  long CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& mat_energy, long index_t);
  
  double cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t);
  
  double cal_kinetic_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t);
  double cal_total_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, long index_t);
  long cal_detail_repulsion(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, const char* fn, long index_t);
  long GET_dCDF_POTENTIAL(TRAJECTORY& TRAJ, long index_t, POTENTIAL_SET& POTs, long index_particle, MATRIX& INDEX_dCDF_U, MATRIX& dCDF_U, MATRIX& vec_boost_ordered_pdf);
  long GET_ORDERED_dCDF_POTENTIAL(TRAJECTORY& TRAJ, long index_t, POTENTIAL_SET& POTs, long index_particle, MATRIX& INDEX_dCDF_U, MATRIX& dCDF_U, MATRIX& vec_boost_ordered_pdf, double& time_dCDF, double& time_SORT);
  namespace ANAL_ASSOCIATION
  {
    long CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, long index_t, MATRIX& tmp_vec);
    
    double cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long index_t);
    double cal_potential_energy_boost(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, long index_t, MATRIX& tmp_vec);
    

  }
  
}

#endif
