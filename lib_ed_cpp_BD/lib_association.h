
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
  long Nc; // number of chains per micelle
  long Tec; // tolerance of end chains
  long N_min; // 2Nc - Tec
  long N_max; // 2Nc + Tec
  long N_ASSOCIATION;
  bool MULTIPLE_CONNECTIONS;
  MATRIX *CASE;
  MATRIX *dPDF;
  MATRIX *dCDF;
  double *Z;

  MATRIX *weight;
  /* MATRIX dist_map; */
  /* MATRIX *dist_map; */
  /* MATRIX dPDF_U; // this is for removing potential overhead for computing (ceiling) PDF for potential */
  

  // member variables
  
  /* bool (*ADD_ASSOCIATION)(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new); */
  /* bool (*MOV_ASSOCIATION)(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new); */
  /* bool (*DEL_ASSOCIATION)(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new); */
  long flag_other, flag_new, flag_itself, flag_hash_other;
  bool (*ADD_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  /* bool (*MOV_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]); */
  bool (*DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CANCEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_ADD_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_MOV_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);

  
  long N_TOTAL_ASSOCIATION();
  long N_CONNECTED_ENDS(long given_index);
  long N_TOTAL_CONNECTED_ENDS();
  bool CHECK_EXIST(long index_A, long index_B);
  bool CHECK_EXIST_1D(long index_A, long index_B);
  long add_association(long index_particle, long index_target);
  long add_association_INFO(POTENTIAL_SET& POTs, long index_particle, long index_target, double distance);
  long del_association(long index_particle, long index_target);
  long del_association_hash(long index_particle, long hash_index_target);
  long del_association_IK(long index_I, long hash_index_K);
  long del_association_grab_IK(long index_I, long hash_index_K);
  
  long GET_INDEX_HASH_FROM_ROLL(long index_particle, double rolled_p);
  long GET_HASH_FROM_ROLL(long index_particle, double rolled_p);
  long FIND_HASH_INDEX(long index_A, long index_B);
 ASSOCIATION() : CONNECTIVITY(){}

  ASSOCIATION(TRAJECTORY& TRAJ, COND& given_condition);
  ASSOCIATION(long number_of_particles, long number_of_chains_per_particles, long tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS);
  long set_initial_condition();
  long initial();
  virtual ~ASSOCIATION()
    {
      mkl_free(CASE);
      mkl_free(dPDF);
      mkl_free(dCDF);
      mkl_free(Z);
      mkl_free(weight);
    }

  long read_exist_weight(const char* fn_weight);
  
};


namespace TRUTH_MAP
{
  long IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);  
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
  
  bool DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  
  bool NEW_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  bool NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  
  bool MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  bool MOV_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
}


double CONNECTIVITY_update_CASE_particle_hash_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, long index_particle, long hash_index_target, double distance);

double CONNECTIVITY_update_CASE_particle_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, long index_particle, long index_target, double distance);
/* long CONNECTIVITY_update_CASE_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, long index_particle, double distance); */
double CONNECTIVITY_update_Z_particle(ASSOCIATION& CONNECT, long index_particle);
double CONNECTIVITY_update_dPDF_particle(ASSOCIATION& CONNECT, long index_particle);
double CONNECTIVITY_update_dCDF_particle(ASSOCIATION& CONNECT, long index_particle);

long backsearch(MATRIX& given_arr, double p);
long bisection_search(MATRIX& given_arr, double p);
  

#endif
