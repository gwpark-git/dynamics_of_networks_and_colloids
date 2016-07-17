
#ifndef DUMBBELL_MODEL_H
#define DUMBBELL_MODEL_H

#include <iomanip>

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"
#include "brownian.h"
#include "repulsive_brownian.h"
#include "omp.h"

namespace DUMBBELL
{
  struct
  DUMBBELL_VARIABLE; // tell compiler to expect struct def.
    
  MKL_LONG
  main_DUMBBELL(TRAJECTORY& TRAJ,
                CONNECTIVITY& CONNECT,
                POTENTIAL_SET& POTs,
                RECORD_DATA& DATA,
                COND& given_condition);



  MKL_LONG
  generate_dumbbell_connectivity(CONNECTIVITY& CONNECT);

  double
  OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now,
                    RDIST& R_boost, MKL_LONG* tmp_index_vec,
                    CONNECTIVITY& CONNECT,
                    const MKL_LONG N_THREADS_BD);

  double
  record_RDIST(ofstream &file_DIST,
	       RDIST& R_boost, CONNECTIVITY& CONNECT);
  
  // double
  //   record_simulation_data(RECORD_DATA& DATA,
  // 			   TRAJECTORY& TRAJ, const MKL_LONG index_t_now,
  // 			   MATRIX& energy);

  double
  report_simulation_info(TRAJECTORY& TRAJ,
                         MATRIX& energy,
                         DUMBBELL_VARIABLE& VAR);
  /* double OMP_compute_RDIST(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, RDIST& R_boost, MKL_LONG* tmp_index_vec, const MKL_LONG N_THREADS_BD); */

  double
  OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next,
                           CONNECTIVITY& CONNECT,
                           POTENTIAL_SET& POTs, MATRIX* force_random, MATRIX* force_spring,
                           RNG_BOOST& RNG,
                           RDIST& R_boost,
                           const MKL_LONG N_THREADS_BD,
                           COND& given_condition,
                           DUMBBELL_VARIABLE& VAR);

  // inline functions
  double
  sum_virial_components(MATRIX& energy);
}

struct
DUMBBELL::
  DUMBBELL_VARIABLE
  : BROWNIAN::BROWNIAN_VARIABLE
{
  // the variables on here is supporting for interface
  // do not use the core variable inside this structure
  // for designing purpose, it is of importance to use individual variable definition outside of this structure

  // additional member variables
  double time_LV_force_connector;
  MATRIX *force_spring;


  // related with virials
  double RF_connector_xx, RF_connector_yy, RF_connector_zz;
  double RF_connector_xy, RF_connector_xz, RF_connector_yz;

  double energy_elastic_potential;
  
  // member functions for virials
  MKL_LONG
  virial_initial();
  double
  record_virial_into_energy_array(MATRIX& energy);

  
  DUMBBELL_VARIABLE()
  {
    std::cout << "ERR: Basic constructor for BROWNIAN_VARIABLE structure is not supported\n";
  }
  DUMBBELL_VARIABLE(COND& given_condition, MKL_LONG given_N_basic)
    : BROWNIAN::BROWNIAN_VARIABLE(given_condition, given_N_basic)
  {
    time_LV_force_connector = 0.;
    force_spring = new MATRIX [Np];

    MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
    for(MKL_LONG i=0; i<Np; i++)
      {
        force_spring[i].initial(N_dimension, 1, 0.);
      }

    virial_initial();
    N_components_energy = 24;
    /* all the necessary information stored in the construction for BROWNIAN::BROWNIAN_VARIABLE */
  }
  
  virtual
  ~DUMBBELL_VARIABLE()
  {
    if(INITIALIZATION)
      delete[] force_spring;
  }
};


inline MKL_LONG
DUMBBELL::DUMBBELL_VARIABLE::
virial_initial()
{

  BROWNIAN::BROWNIAN_VARIABLE::
    virial_initial();

  RF_connector_xx = 0.; RF_connector_yy = 0.; RF_connector_zz = 0.;
  RF_connector_xy = 0.; RF_connector_xz = 0.; RF_connector_yz = 0.;

  energy_elastic_potential = 0.;
  return 0;

}

inline double
DUMBBELL::DUMBBELL_VARIABLE::
record_virial_into_energy_array(MATRIX& energy)
{
  double time_st = dsecnd();
  BROWNIAN::BROWNIAN_VARIABLE::
    record_virial_into_energy_array(energy);

  // note that the index number in here is temporary used
  // it might occure confusion between connector virial component compared with stochastic_HEUR
  // however, to keep homogeneous data structure, monotonic increasing function for index number is used.
  
  energy(18) = RF_connector_xx/(2.*volume_PBC_box);
  energy(19) = RF_connector_yy/(2.*volume_PBC_box);
  energy(20) = RF_connector_zz/(2.*volume_PBC_box);

  energy(21) = RF_connector_xy/(2.*volume_PBC_box);
  energy(22) = RF_connector_xz/(2.*volume_PBC_box);
  energy(23) = RF_connector_yz/(2.*volume_PBC_box);

  energy(3) = energy_elastic_potential;
  
  return dsecnd() - time_st;
}

inline double
DUMBBELL::
sum_virial_components(MATRIX& energy)
{
  double time_st = dsecnd();
  BROWNIAN::sum_virial_components(energy);
  for(MKL_LONG i=6; i<12; i++)
    {
      energy(i) += energy(i + 6*2);
    }

  energy(1) = energy(2) + energy(3); // related with total energy of the system
  
  return dsecnd() - time_st;
}


#endif
