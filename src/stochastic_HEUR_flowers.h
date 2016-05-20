
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

struct TEMPORAL_VARIABLE_HEUR : REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE
{
  MKL_LONG N_steps_block;
  MKL_LONG N_THREADS_SS;
  double time_SS, time_SS_CORE;
  double time_SS_index, time_SS_LOCK, time_SS_check, time_SS_transition, time_SS_update;
  MKL_LONG cnt_arr[6];
  MKL_LONG N_associations;
  
  MATRIX *force_spring;

  MKL_LONG N_tot_associable_chain;

  // checking input files
  bool MC_renewal;
  bool MC_LOG;
  
  TEMPORAL_VARIABLE_HEUR()
    {};
 TEMPORAL_VARIABLE_HEUR(COND& given_condition, MKL_LONG given_N_basic) : REPULSIVE_BROWNIAN::TEMPORAL_VARIABLE(given_condition, given_N_basic)
    {
      MKL_LONG Np = atoi(given_condition("Np").c_str());
      MKL_LONG N_dimension = atoi(given_condition("N_dimension").c_str());
      N_steps_block = atol(given_condition("N_steps_block").c_str());
      N_THREADS_SS = atol(given_condition("N_THREADS_SS").c_str());
      time_SS = 0.;
      time_SS_CORE = 0.;
      time_SS_update = 0.;
      time_SS_index = 0.;
      time_SS_LOCK = 0.;
      time_SS_check = 0.;
      time_SS_transition = 0.;
      
      for(MKL_LONG i=0; i<6; i++)
        cnt_arr[i] = 0;
      
      force_spring = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
      
      for(MKL_LONG i=0; i<Np; i++)
        force_spring[i].initial(N_dimension, 1, 0.);

      N_tot_associable_chain = Np*atoi(given_condition("N_chains_per_particle").c_str());

      MC_renewal = FALSE;
      if(given_condition("MC_renewal")=="TRUE")
        MC_renewal = TRUE;

      MC_LOG = FALSE;
      if(given_condition("MC_LOG")=="TRUE")
        MC_LOG = TRUE;
    }
 ~TEMPORAL_VARIABLE_HEUR()
    {
      if(INITIALIZATION)
        {
          mkl_free(force_spring);
          
        }
    }
};


MKL_LONG stochastic_simulation_HEUR_flowers(TRAJECTORY& TRAJ, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, COND& given_condition);

double record_simulation_data(RECORD_DATA& DATA, TRAJECTORY& TRAJ, CHAIN_HANDLE& CHAIN, MATRIX& energy, const MKL_LONG index_t_now) ;

double report_simulation_info(TRAJECTORY& TRAJ, MATRIX& energy, TEMPORAL_VARIABLE_HEUR& VAR);
double OMP_time_evolution_Euler(TRAJECTORY& TRAJ, const MKL_LONG index_t_now, const MKL_LONG index_t_next, POTENTIAL_SET& POTs, RDIST& R_boost, MATRIX* vec_boost_Nd_parallel, MATRIX* force_repulsion, MATRIX* force_random, MATRIX* force_spring, RNG_BOOST& RNG, const MKL_LONG N_THREADS_BD, COND& given_condition);

double write_MC_LOG_if_TRUE(bool flag_MC_LOG, RECORD_DATA& DATA, ASSOCIATION& CONNECT, const MKL_LONG cnt, const MKL_LONG index_itself, const double rolling_dCDF, const MKL_LONG index_attached_bead, const MKL_LONG index_new_attached_bead, const double rolling_dCDF_U, MKL_LONG* cnt_arr);

double transition_single_chain_end(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, INDEX_MC& IDX);

double check_dissociation_probability(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, INDEX_MC& IDX, RNG_BOOST& RNG);


double LOCKING_PARALLEL(LOCK& LOCKER, TEMPORAL_VARIABLE_HEUR& VAR, const MKL_LONG index_thread);
double release_LOCKING(LOCK& LOCKER, INDEX_MC& IDX);

double micelle_selection(ASSOCIATION& CONNECT, RNG_BOOST& RNG, MKL_LONG& index_itself, MKL_LONG& index_hash_attached_bead, MKL_LONG& index_attached_bead, MKL_LONG& index_new_attached_bead, const MKL_LONG index_thread);

double OMP_SS_update_topology(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, CHAIN_HANDLE& CHAIN, RECORD_DATA& DATA, INDEX_MC* IDX_ARR, TEMPORAL_VARIABLE_HEUR& VAR);

double OMP_SS_update_STATISTICS(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, RDIST& R_boost, TEMPORAL_VARIABLE_HEUR& VAR);

double OMP_SS_topological_time_evolution(ASSOCIATION& CONNECT, CHAIN_HANDLE& CHAIN, POTENTIAL_SET& POTs, RECORD_DATA& DATA, COND& given_condition, TEMPORAL_VARIABLE_HEUR& VAR);








#endif
