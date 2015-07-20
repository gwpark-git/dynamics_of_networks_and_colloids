
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
#include "matrix_long_ed.h"
#include "lib_connectivity.h"
#include "lib_potential.h"

class ASSOCIATION : public CONNECTIVITY
{
 public:
  MKL_LONG Nc; // number of chains per micelle
  MKL_LONG Tec; // tolerance of end chains
  MKL_LONG N_min; // 2Nc - Tec
  MKL_LONG N_max; // 2Nc + Tec
  MKL_LONG N_ASSOCIATION;
  /* double *CASE_function_variables; */
  /* double (*CASE_function)(double distance, double* given_variables); */
  MATRIX CASE;
  MATRIX dPDF;
  MATRIX dCDF;
  MATRIX Z;

  MATRIX_LONG weight;

  MATRIX dist_map;

  /* MATRIX dPDF_U; // this is for removing potential overhead for computing (ceiling) PDF for potential */
  

  // member variables
  
  /* double& N_CHAIN_ENDS(MKL_LONG index_particle); */

  bool (*NEW_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool (*MOV_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool (*DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);

  MKL_LONG N_TOTAL_ASSOCIATION();
  MKL_LONG N_CONNECTED_ENDS(MKL_LONG given_index);
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
  /* bool CHECK_NC_ADD_BTA_1D(MKL_LONG index_A, MKL_LONG index_B); */
  /* bool CHECK_NC_ADD_BTA(MKL_LONG index_A, MKL_LONG index_B); */
  /* bool CHECK_NC_DEL_hBFA(MKL_LONG index_A, MKL_LONG index_hash_B); */
  /* bool CHECK_NC_DEL_BFA(MKL_LONG index_A, MKL_LONG index_B); */
  /* bool CHECK_NC_DEL_BFA_1D(MKL_LONG index_A, MKL_LONG index_B); */
 ASSOCIATION() : CONNECTIVITY(){}
  
  ASSOCIATION(MKL_LONG number_of_particles, MKL_LONG number_of_chains_per_particles, MKL_LONG tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS);
  MKL_LONG initial();
  virtual ~ASSOCIATION()
    {
    }

  
};


namespace TRUTH_MAP
{

  bool DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  
  bool NEW_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  
  bool MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool MOV_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
}

double CONNECTIVITY_update_CASE_particle_hash_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG hash_index_target, double distance);

double CONNECTIVITY_update_CASE_particle_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG index_target, double distance);
/* MKL_LONG CONNECTIVITY_update_CASE_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, double distance); */
double CONNECTIVITY_update_Z_particle(ASSOCIATION& CONNECT, MKL_LONG index_particle);
double CONNECTIVITY_update_dPDF_particle(ASSOCIATION& CONNECT, MKL_LONG index_particle);
double CONNECTIVITY_update_dCDF_particle(ASSOCIATION& CONNECT, MKL_LONG index_particle);
/* MKL_LONG CONNECTIVITY_update_information_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle); */
/* MKL_LONG CONNECTIVITY_update_information(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs); */


#endif
