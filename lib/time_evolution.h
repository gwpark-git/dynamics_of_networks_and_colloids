
#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H

#include "geometry.h"
#include "potential.h"
#include "fstream"
#include "association.h"
#include <mkl.h>

namespace INTEGRATOR
{
  namespace EULER
  {
    double
      cal_repulsion_force_R_boost
      (POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_particle, RDIST& R_boost);

    double
    cal_repulsion_force_R_boost_with_RF
    (POTENTIAL_SET& POTs, MATRIX& given_vec, MKL_LONG index_particle, RDIST& R_boost,
     double& RF_repulsion_xx, double& RF_repulsion_yy, double& RF_repulsion_zz,
     double& RF_repulsion_xy, double& RF_repulsion_xz, double& RF_repulsion_yz);

    
    double
      cal_random_force_boost
      (POTENTIAL_SET& POTs, MATRIX& given_vec, gsl_rng* r_boost);
    double
      cal_random_force_boost_simplified
      (POTENTIAL_SET& POTs, MATRIX& given_vec, gsl_rng* r_boost);

    double
      cal_connector_force_boost
      (POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, MATRIX& given_vec, MKL_LONG given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost);

    double
    cal_connector_force_boost_with_RF
    (POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, MATRIX& given_vec, MKL_LONG given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost,
     double& RF_connector_xx, double& RF_connector_yy, double& RF_connector_zz,
     double& RF_connector_xy, double& RF_connector_xz, double& RF_connector_yz);

    
  }
  namespace EULER_ASSOCIATION
  {
    double
      cal_connector_force_boost
      (POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost);

    double
    cal_connector_force_boost_with_RF
    (POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& given_vec, MKL_LONG given_index, MATRIX** R_minimum_vec_boost, MATRIX* R_minimum_distance_boost,
     double& RF_connector_xx, double& RF_connector_yy, double& RF_connector_zz,
     double& RF_connector_xy, double& RF_connector_xz, double& RF_connector_yz);

    
  }

}

namespace ANALYSIS
{
  double
    CAL_ENERGY_BROWNIAN
    (POTENTIAL_SET& POTs, MATRIX& mat_energy, double time);
  double
    CAL_ENERGY_R_boost
    (POTENTIAL_SET& POTs, MATRIX& mat_energy, double time, RDIST& R_boost);

  double
    cal_potential_energy_R_boost
    (POTENTIAL_SET& POTs, RDIST& R_boost);
  double
    cal_total_energy_R_boost
    (POTENTIAL_SET& POTs, RDIST& R_boost);
  namespace DUMBBELL
  {
    double
      CAL_ENERGY_R_boost(POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, MATRIX& mat_energy, double time, RDIST& R_boost); // note that this functionality is subect to re-define. It is due to the fact that RDIST is not necessary for the dumbbell model, but it is used for using the same interface with stochastic HEUR simulation
    double
      cal_potential_energy_R_boost(POTENTIAL_SET& POTs, CONNECTIVITY& CONNECT, RDIST& R_boost);
  }
  namespace ANAL_ASSOCIATION
  {
    double
      CAL_ENERGY_R_boost
      (POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MATRIX& mat_energy, double time, RDIST& R_boost);
    double
      cal_potential_energy_R_boost
      (POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, RDIST& R_boost);
  }
}

#endif
