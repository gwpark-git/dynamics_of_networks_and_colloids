
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
#include "lib_traj.h"
#include "read_file_condition.h"

/* class MC_IDENTIFIER */
/* { */
/*  public: */
/*   double IDENTIFIER; */
/*   double FACTOR_PRE, FACTOR_NOW; */
/*   double TOLERANCE; */
/*   MC_IDENTIFIER() */
/*     { */
/*       BASIC_INIT(); */
/*       TOLERANCE = 0.01; */
/*     } */
/*   MC_IDENTIFIER(double given_tolerance) */
/*     { */
/*       BASIC_INIT(); */
/*       TOLERANCE = given_tolerance; */
/*     } */
/*   long BASIC_INIT() */
/*     { */
/*       IDENTIFIER=1.; */
/*       FACTOR_PRE=0.; */
/*       FACTOR_NOW=0.; */
/*       return 0; */
/*     } */
/*   ~MC_IDENTIFIER(){} */
   
/* }; */


class ASSOCIATION : public CONNECTIVITY
{
 public:
  long Nc; // number of chains per micelle
  long Tec; // tolerance of end chains
  long N_min; // 2Nc - Tec
  long N_max; // 2Nc + Tec
  long N_ASSOCIATION;
  /* MATRIX CASE; */
  /* MATRIX dPDF; */
  /* MATRIX dCDF; */
  /* MATRIX Z; */

  MATRIX *CASE;
  MATRIX *dPDF;
  MATRIX *dCDF;
  double *Z;

  /* MATRIX_LONG weight; */
  /* MATRIX_LONG * weight; */
  MATRIX *weight;
  /* MATRIX dist_map; */
  /* MATRIX *dist_map; */
  /* MATRIX dPDF_U; // this is for removing potential overhead for computing (ceiling) PDF for potential */
  

  // member variables
  
  /* double& N_CHAIN_ENDS(long index_particle); */

  bool (*NEW_ASSOCIATION)(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  bool (*MOV_ASSOCIATION)(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);
  bool (*DEL_ASSOCIATION)(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new);

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
  /* bool CHECK_NC_ADD_BTA_1D(long index_A, long index_B); */
  /* bool CHECK_NC_ADD_BTA(long index_A, long index_B); */
  /* bool CHECK_NC_DEL_hBFA(long index_A, long index_hash_B); */
  /* bool CHECK_NC_DEL_BFA(long index_A, long index_B); */
  /* bool CHECK_NC_DEL_BFA_1D(long index_A, long index_B); */
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
/* long CONNECTIVITY_update_information_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, long index_particle); */
/* long CONNECTIVITY_update_information(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs); */


long backsearch(MATRIX& given_arr, double p);
long bisection_search(MATRIX& given_arr, double p);
  

#endif
