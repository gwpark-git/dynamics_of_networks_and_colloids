
#ifndef LIB_ASSOCIATION_H
#define LIB_ASSOCIATION_H

#define TRUE 1
#define FALSE 0

#include <iostream>
#include <omp.h>
#include <mkl.h>

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}
#include "matrix_ed.h"
#include "lib_connectivity.h"
#include "lib_potential.h"
#include "lib_traj.h"
#include "read_file_condition.h"


class ASSOCIATION : public CONNECTIVITY
{
 public:
  MKL_LONG Nc; // number of chains per micelle
  MKL_LONG Tec; // tolerance of end chains
  MKL_LONG N_min; // 2Nc - Tec
  MKL_LONG N_max; // 2Nc + Tec
  MKL_LONG N_ASSOCIATION;
  bool MULTIPLE_CONNECTIONS;
  MATRIX *CASE;
  MATRIX *dPDF;
  MATRIX *dCDF;
  double *Z;

  MATRIX *weight;

  // member variables
  
  /* MKL_LONG flag_other, flag_new, flag_itself, flag_hash_other; */
  static const MKL_LONG flag_other = 0;
  static const MKL_LONG flag_new = 1;
  static const MKL_LONG flag_itself = 2;
  static const MKL_LONG flag_hash_other = 3;
  bool (*ADD_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CANCEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_ADD_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_MOV_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);

  
  MKL_LONG N_TOTAL_ASSOCIATION();
  MKL_LONG N_CONNECTED_ENDS(MKL_LONG given_index);
  MKL_LONG N_TOTAL_CONNECTED_ENDS();
  bool CHECK_EXIST(MKL_LONG index_A, MKL_LONG index_B);
  bool CHECK_EXIST_1D(MKL_LONG index_A, MKL_LONG index_B);
  MKL_LONG add_association(MKL_LONG index_particle, MKL_LONG index_target);
  MKL_LONG add_association_INFO(POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG index_target, double distance);
  MKL_LONG del_association(MKL_LONG index_particle, MKL_LONG index_target);
  MKL_LONG del_association_hash(MKL_LONG index_particle, MKL_LONG hash_index_target);
  MKL_LONG del_association_IK(MKL_LONG index_I, MKL_LONG hash_index_K);
  MKL_LONG del_association_grab_IK(MKL_LONG index_I, MKL_LONG hash_index_K);
  
  MKL_LONG GET_INDEX_HASH_FROM_ROLL(MKL_LONG index_particle, double rolled_p);
  MKL_LONG GET_HASH_FROM_ROLL(MKL_LONG index_particle, double rolled_p);
  MKL_LONG FIND_HASH_INDEX(MKL_LONG index_A, MKL_LONG index_B);


  // At the moment, there is no namespace dependence between KINETICS::NORMALIZED and KINETICS::METROPOLIS.
  // The dependency is set with the w_function and transtion_probability that is developed with lib_potential.
  // For better performance, it would be better to tune the following with specific dependencies (future works).
  double update_CASE_particle_hash_target(POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG hash_index_target, double distance);

  double update_CASE_particle_target(POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG index_target, double distance);
  /* MKL_LONG CONNECTIVITY_update_CASE_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, double distance); */
  double update_Z_particle(MKL_LONG index_particle);
  double update_dPDF_particle(MKL_LONG index_particle);
  double update_dCDF_particle(MKL_LONG index_particle);


  
 ASSOCIATION() : CONNECTIVITY(){}

  ASSOCIATION(TRAJECTORY& TRAJ, COND& given_condition);
  ASSOCIATION(MKL_LONG number_of_particles, MKL_LONG number_of_chains_per_particles, MKL_LONG tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS);
  MKL_LONG set_initial_condition();
  MKL_LONG initial();
  virtual ~ASSOCIATION()
    {
      mkl_free(CASE);
      mkl_free(dPDF);
      mkl_free(dCDF);
      mkl_free(Z);
      mkl_free(weight);
    }

  MKL_LONG read_exist_weight(const char* fn_weight);
  
};


namespace TRUTH_MAP
{
  namespace MULTIPLE
  {
    bool CHECK_N_ADD_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_DEL_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_MOV_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  }
  namespace SINGLE
  {
    bool CHECK_N_ADD_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_MOV_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_DEL_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  }
  bool CANCEL_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool ADD_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool DEL_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  
  bool DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  
  bool NEW_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  
  bool MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool MOV_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
}



namespace KINETICS
{
  double CONNECTIVITY_update_dCDF_particle(ASSOCIATION* CONNECT, MKL_LONG index_particle);
  double CONNECTIVITY_update_CASE_particle_hash_target(ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG hash_index_target, double distance);
  double CONNECTIVITY_update_CASE_particle_target(ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG index_target, double distance);

  double CONNECTIVITY_update_Z_particle(ASSOCIATION* CONNECT, MKL_LONG index_particle);
  double CONNECTIVITY_update_dPDF_particle(ASSOCIATION* CONNECT, MKL_LONG index_particle);
}

  

#endif
