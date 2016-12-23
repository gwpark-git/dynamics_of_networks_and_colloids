
#ifndef STOCHASTIC_HEUR_FLOWERS_H
#define STOCHASTIC_HEUR_FLOWERS_H

#include "../lib/matrix.h"
#include "../lib/trajectory.h"
#include "../lib/time_evolution.h"
#include "../lib/association.h"
#include "../lib/handle_association.h"
#include "../lib/potential.h"
#include "../lib/parallel.h"
#include "../lib/geometry.h"

#include "repulsive_brownian.h"
#include "omp.h"

namespace HEUR
{
  struct
  TEMPORAL_VARIABLE_HEUR;

  MKL_LONG
  stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ,
				     POTENTIAL_SET& POTs,
				     ASSOCIATION& CONNECT,
				     CHAIN_HANDLE& CHAIN,
				     RECORD_DATA& DATA,
				     COND& given_condition);

  double
  record_simulation_data(RECORD_DATA& DATA,
			 TRAJECTORY& TRAJ,
			 ASSOCIATION& CONNECT,
			 CHAIN_HANDLE& CHAIN,
			 const MKL_LONG index_t_now);

  double
  record_RDIST_BRIDGE(ofstream &file_DIST,
                      RDIST& R_boost, ASSOCIATION& CONNECT);

  double
  record_RDIST_PARTICLE(ofstream &file_DIST,
                        RDIST& R_boost);

    
  double
  report_simulation_info(TRAJECTORY& TRAJ,
			 MATRIX& energy,
			 TEMPORAL_VARIABLE_HEUR& VAR);

  double
  OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next,
			   ASSOCIATION& CONNECT,
			   POTENTIAL_SET& POTs, MATRIX* force_repulsion, MATRIX* force_random, MATRIX* force_spring,
			   RDIST& R_boost, MATRIX* vec_boost_Nd_parallel,
			   RNG_BOOST& RNG,
			   const MKL_LONG N_THREADS_BD,
			   COND& given_condition, TEMPORAL_VARIABLE_HEUR& VAR);

  double
  write_MC_LOG_if_TRUE(bool flag_MC_LOG,
		       RECORD_DATA& DATA,
		       ASSOCIATION& CONNECT,
		       const INDEX_MC& IDX,
		       const MKL_LONG cnt, const MKL_LONG* cnt_arr,
		       const double rolling_dCDF, const double rolling_dCDF_U);

  double
  transition_single_chain_end(ASSOCIATION& CONNECT,
			      POTENTIAL_SET& POTs,
			      CHAIN_HANDLE& CHAIN,
			      RDIST& R_boost,
			      INDEX_MC& IDX, const MKL_LONG IDENTIFIER_ACTION);

  double
  check_dissociation_probability(ASSOCIATION& CONNECT,
				 POTENTIAL_SET& POTs,
				 RDIST& R_boost,
				 INDEX_MC& IDX,
				 gsl_rng* RNG_BOOST_SS_IT,
				 MKL_LONG& IDENTIFIER_ACTION);


  double
  LOCKING_PARALLEL(LOCK& LOCKER,
		   TEMPORAL_VARIABLE_HEUR& VAR,
		   const INDEX_MC& IDX,
		   MKL_LONG& IDENTIFIER_ACTION, MKL_LONG& IDENTIFIER_LOCKING);
  double
  release_LOCKING(LOCK& LOCKER, INDEX_MC& IDX);

  double
  micelle_selection(ASSOCIATION& CONNECT,
		    gsl_rng* RNG_BOOST_SS_IT,
		    INDEX_MC& IDX,
		    RDIST& R_boost,
		    TEMPORAL_VARIABLE_HEUR& VAR,
		    double& rolling_dCDF, double& rolling_dCDF_U, MKL_LONG& time_step);
  bool
  get_index_bridge_distance_larger_than_cutoff_range(MKL_LONG& index_particle, MKL_LONG& index_hash, MKL_LONG& index_target,
                                                     ASSOCIATION& CONNECT,
                                                     RDIST& R_boost,
                                                     MKL_LONG cutoff_range);

  double
  OMP_SS_update_topology(ASSOCIATION& CONNECT,
			 POTENTIAL_SET& POTs,
			 RDIST& R_boost,
			 CHAIN_HANDLE& CHAIN,
			 RNG_BOOST& RNG,
			 RECORD_DATA& DATA,
			 INDEX_MC* IDX_ARR,
			 LOCK& LOCKER,
			 TEMPORAL_VARIABLE_HEUR& VAR);

  double
  OMP_SS_update_STATISTICS(ASSOCIATION& CONNECT,
			   POTENTIAL_SET& POTs,
			   RDIST& R_boost,
			   TEMPORAL_VARIABLE_HEUR& VAR,
			   RECORD_DATA& DATA);

  double
  OMP_SS_topological_time_evolution(const MKL_LONG time_step_LV,
				    ASSOCIATION& CONNECT,
				    CHAIN_HANDLE& CHAIN,
				    POTENTIAL_SET& POTs,
				    RDIST& R_boost,
				    RNG_BOOST& RNG,
				    INDEX_MC* IDX_ARR,
				    RECORD_DATA& DATA,
				    COND& given_condition,
				    LOCK& LOCKER,
				    TEMPORAL_VARIABLE_HEUR& VAR);

  // inline functions
  double
  sum_virial_components(MATRIX& energy);

}

struct
HEUR::
  TEMPORAL_VARIABLE_HEUR
  : REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE
{
  MKL_LONG N_steps_block;
  MKL_LONG N_THREADS_SS;
  double time_SS, time_SS_CORE;
  double time_SS_index, time_SS_LOCK, time_SS_check, time_SS_transition, time_SS_update_info;
  double time_SS_update, time_SS_update_ASSOCIATION_MAP, time_SS_update_CHAIN_SUGGESTION;

  double time_LV_force_repulsion, time_LV_force_random, time_LV_force_connector;
  MKL_LONG cnt_arr[6];
  MKL_LONG cnt_SS; // counting number of steps for stochastic simulations
  MKL_LONG N_associations;
  
  MATRIX *force_spring;

  MKL_LONG N_tot_associable_chain;
  MKL_LONG N_diff_associations;
  
  // checking input files
  bool MC_renewal;
  bool MC_LOG;

  // for debugging or analysis of instance of time
  bool MC_ASSOCIATION;

  // for cutoff_scheme
  bool CUTOFF_BRIDGES;
  bool SIGN_DESTROY;
  double CUTOFF_RANGE;
  
  // related with virials
  double RF_connector_xx, RF_connector_yy, RF_connector_zz;
  double RF_connector_xy, RF_connector_xz, RF_connector_yz;

  double energy_elastic_potential;
  
  // member functions for virials
  MKL_LONG
  virial_initial();
  double
  record_virial_into_energy_array(MATRIX& energy);
  
  TEMPORAL_VARIABLE_HEUR()
  {};
  TEMPORAL_VARIABLE_HEUR(COND& given_condition, MKL_LONG given_N_basic);
  
  ~TEMPORAL_VARIABLE_HEUR()
  {
    if(INITIALIZATION)
      {
	delete[] force_spring;
      }
  }
};

inline MKL_LONG
HEUR::TEMPORAL_VARIABLE_HEUR::
virial_initial()
{
  // note that this function is not inherite from previous function, directly
  // rather than, it explicitly call the function defined in the function in parents class
  REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::
    virial_initial();
  
  RF_connector_xx = 0.; RF_connector_yy = 0.; RF_connector_zz = 0.;
  RF_connector_xy = 0.; RF_connector_xz = 0.; RF_connector_yz = 0.;

  energy_elastic_potential = 0.;
  return 0;
}

inline double
HEUR::TEMPORAL_VARIABLE_HEUR::
record_virial_into_energy_array(MATRIX& energy)
{
  double time_st = dsecnd();

  REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE::
    record_virial_into_energy_array(energy);
  
  energy(24) = RF_connector_xx;
  energy(25) = RF_connector_yy;
  energy(26) = RF_connector_zz;

  energy(27) = RF_connector_xy;
  energy(28) = RF_connector_xz;
  energy(29) = RF_connector_yz;

  // energy(3) = energy_repulsive_potential + energy_elastic_potential;
  // energy(2) = energy_repulsive_potential; // it is already in the REPULSIVE_BROWNIAN
  energy(3) = energy_elastic_potential;
  return time_st - dsecnd();
}
				
inline double
HEUR::
sum_virial_components(MATRIX& energy)
{
  double time_st = dsecnd();
  REPULSIVE_BROWNIAN::sum_virial_components(energy);
  
  for(MKL_LONG i=6; i<12; i++)
    {
      energy(i) += energy(i + 6*3);
    }

  energy(1) = energy(2) + energy(3);
  
  return dsecnd() - time_st;
}


#endif
