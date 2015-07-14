
#ifndef LIB_EVOLUTION_H
#define LIB_EVOLUTION_H
#include "lib_traj.h"
#include "lib_geometry.h"
#include "lib_potential.h"
/* #include "mkl_vsl.h" */
#include "fstream"
/* #include "potential_definition.h" */
#include "lib_association.h"

namespace INTEGRATOR
{
  namespace EULER
  {
    MKL_LONG simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
    MKL_LONG cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index);
    MKL_LONG cal_repulsion_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG index_i);
    MKL_LONG cal_random_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_t);
  }
  namespace EULER_ASSOCIATION
  {
    /* MKL_LONG simple_Euler(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t); */
    MKL_LONG cal_connector_force(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG index_t, MKL_LONG given_index);
  }

}

namespace ANALYSIS
{

  MKL_LONG CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MATRIX& mat_energy, MKL_LONG index_t);
  double cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
  double cal_kinetic_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
  double cal_total_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, MKL_LONG index_t);
  MKL_LONG cal_detail_repulsion(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, const char* fn, MKL_LONG index_t);
  MKL_LONG GET_PDF_POTENTIAL(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG index_particle, MATRIX& INDEX_PDF_U, MATRIX& PDF_U);
  MKL_LONG GET_ORDERED_PDF_POTENTIAL(TRAJECTORY& TRAJ, MKL_LONG index_t, POTENTIAL_SET& POTs, MKL_LONG index_particle, MATRIX& INDEX_PDF_U, MATRIX& PDF_U);
  namespace ANAL_ASSOCIATION
  {
    MKL_LONG CAL_ENERGY(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, MKL_LONG index_t);
    double cal_potential_energy(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG index_t);


  }
  
}

#endif
