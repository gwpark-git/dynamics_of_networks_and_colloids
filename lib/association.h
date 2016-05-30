
#ifndef ASSOCIATION_H
#define ASSOCIATION_H

#define TRUE 1
#define FALSE 0

#include <iostream>
#include <omp.h>
#include <mkl.h>

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_const_mksa.h>
}
#include "matrix.h"
#include "connectivity.h"
#include "potential.h"
/* #include "trajectory.h" */
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "geometry.h"

class ASSOCIATION : public CONNECTIVITY
{
 public:
  MKL_LONG Nc; // number of chains per micelle
  MKL_LONG Tec; // tolerance of end chains
  MKL_LONG N_min; // 2Nc - Tec
  MKL_LONG N_max; // 2Nc + Tec
  MKL_LONG N_ASSOCIATION;
  bool MULTIPLE_CONNECTIONS;

  MATRIX *weight;
  // related with suggestion probability (selecting a degenerated pair)
  MATRIX *CASE_SUGGESTION;
  MATRIX *dPDF_SUGGESTION;
  MATRIX *dCDF_SUGGESTION;
  double *Z_SUGGESTION;

  // related with association probability
  MATRIX *dCDF_ASSOCIATION;
  /* size_t **INDEX_ASSOCIATION; */
  MATRIX *INDEX_ASSOCIATION;
  MKL_LONG *TOKEN_ASSOCIATION;

  // member variables
  
  /* MKL_LONG flag_other, flag_new, flag_itself, flag_hash_other; */
  static const MKL_LONG flag_other = 0;
  static const MKL_LONG flag_new = 1;
  static const MKL_LONG flag_itself = 2;
  static const MKL_LONG flag_hash_other = 3;
  static const MKL_LONG flag_hash_backtrace = 4;
  
  bool (*ADD_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*OPP_DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CANCEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_ADD_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_MOV_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool (*CHECK_N_OPP_DEL_ASSOCIATION)(ASSOCIATION& CONNECT, MKL_LONG index_set[]);

  
  MKL_LONG N_TOTAL_ASSOCIATION();
  MKL_LONG N_CONNECTED_ENDS(MKL_LONG given_index);
  MKL_LONG N_TOTAL_CONNECTED_ENDS();
  bool CHECK_EXIST(MKL_LONG index_A, MKL_LONG index_B);
  bool CHECK_EXIST_1D(MKL_LONG index_A, MKL_LONG index_B);
  MKL_LONG add_association(const MKL_LONG index_particle, MKL_LONG index_target);
  MKL_LONG add_association_INFO(POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG index_target, double distance);
  MKL_LONG opp_del_association(const MKL_LONG index_particle, MKL_LONG index_target);
  MKL_LONG opp_del_association_hash(const MKL_LONG index_particle, MKL_LONG hash_index_target);
  MKL_LONG opp_del_association_IK(MKL_LONG index_I, MKL_LONG hash_index_K);
  MKL_LONG opp_del_association_grab_IK(MKL_LONG index_I, MKL_LONG hash_index_K);
  
  MKL_LONG GET_INDEX_HASH_FROM_ROLL(const MKL_LONG index_particle, double rolled_p);
  MKL_LONG GET_HASH_FROM_ROLL(const MKL_LONG index_particle, double rolled_p);
  MKL_LONG FIND_HASH_INDEX(MKL_LONG index_A, MKL_LONG index_B);


  // At the moment, there is no namespace dependence between KINETICS::NORMALIZED and KINETICS::METROPOLIS.
  // The dependency is set with the w_function and transtion_probability that is developed with lib_potential.
  // For better performance, it would be better to tune the following with specific dependencies (future works).
  double update_CASE_SUGGESTION_particle_hash_target(POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG hash_index_target, double distance);

  double update_CASE_SUGGESTION_particle_target(POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG index_target, double distance);
  /* MKL_LONG CONNECTIVITY_update_CASE_SUGGESTION_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, const MKL_LONG index_particle, double distance); */
  double update_Z_SUGGESTION_particle(const MKL_LONG index_particle);
  double update_dPDF_SUGGESTION_particle(const MKL_LONG index_particle);
  double update_dCDF_SUGGESTION_particle(const MKL_LONG index_particle);
  double update_CHAIN_SUGGESTION_MAP_particle(const MKL_LONG index_particle, POTENTIAL_SET& POTs, RDIST& R_boost);


  double update_ASSOCIATION_MAP_particle(const MKL_LONG index_particle, POTENTIAL_SET& POTs, RDIST& R_boost);

  MKL_LONG return_multiplicity(const MKL_LONG index_particle, const MKL_LONG target_particle)
  // this additional function is useful for tracking individual chain
  {
    for(MKL_LONG j=0; j<TOKEN[index_particle]; j++)
      {
        if(HASH[index_particle](j) == target_particle)
          return weight[index_particle](j);
      }
    return 0;
  }
  
 ASSOCIATION() : CONNECTIVITY(){}

  ASSOCIATION(COND& given_condition);
  
  ASSOCIATION(MKL_LONG number_of_particles, MKL_LONG number_of_chains_per_particles, MKL_LONG tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS);
  MKL_LONG dynamic_alloc();
  MKL_LONG set_initial_condition();
  MKL_LONG initial();
  MKL_LONG initial_inheritance();
  virtual ~ASSOCIATION()
    {
      mkl_free(weight);

      // related with suggestion probability
      mkl_free(CASE_SUGGESTION);
      mkl_free(dPDF_SUGGESTION);
      mkl_free(dCDF_SUGGESTION);
      mkl_free(Z_SUGGESTION);

      // related with association probability
      mkl_free(dCDF_ASSOCIATION);
      mkl_free(INDEX_ASSOCIATION);
      /* for(MKL_LONG i=0; i<Np; i++) */
      /*   mkl_free(INDEX_ASSOCIATION[i]); */
      /* mkl_free(INDEX_ASSOCIATION); */
      
      mkl_free(TOKEN_ASSOCIATION);
    }

  MKL_LONG read_exist_weight(const char* fn_weight);
  
};

MKL_LONG index_set_val(size_t* index, MKL_LONG given_val, MKL_LONG size_of_arr);


namespace TRUTH_MAP
{
  namespace MULTIPLE
  {
    bool CHECK_N_ADD_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_OPP_DEL_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_MOV_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  }
  namespace SINGLE
  {
    bool CHECK_N_ADD_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_MOV_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
    bool CHECK_N_OPP_DEL_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  }
  bool CANCEL_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool ADD_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  bool OPP_DEL_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[]);
  
  bool OPP_DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  
  bool NEW_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  
  bool MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
  bool MOV_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new);
}



namespace KINETICS
{
  double CONNECTIVITY_update_dCDF_SUGGESTION_particle(ASSOCIATION* CONNECT, const MKL_LONG index_particle);
  double CONNECTIVITY_update_CASE_SUGGESTION_particle_hash_target(ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG hash_index_target, double distance);
  double CONNECTIVITY_update_CASE_SUGGESTION_particle_target(ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG index_target, double distance);

  double CONNECTIVITY_update_Z_SUGGESTION_particle(ASSOCIATION* CONNECT, const MKL_LONG index_particle);
  double CONNECTIVITY_update_dPDF_SUGGESTION_particle(ASSOCIATION* CONNECT, const MKL_LONG index_particle);
}


class CHAIN_INFORMATION
{
  /*
    The CHAIN class is using CHAIN_NODE without inheritance.
    It contains dynamically allocated CHAIN_NODE.
  */
 public:

  CHAIN_NODE *CHAIN;
  MKL_LONG N_chains;
  MKL_LONG N_particles;

  bool INITIALIZATION;
  MKL_LONG& HEAD(MKL_LONG given_chain_index)
    {
      return CHAIN[given_chain_index].HEAD;
    }
  MKL_LONG& TAIL(MKL_LONG given_chain_index)
    {
      return CHAIN[given_chain_index].TAIL;
    }

  /*
    On the same reason, the ATTACHED will return pointer variable rather than reference variable. This break the existing interface.
  */
  MKL_LONG& ATTACHED(MKL_LONG given_chain_index, MKL_LONG flag_HEAD_TAIL)
    {
      return CHAIN[given_chain_index].index(flag_HEAD_TAIL);
    }
  
  
  MKL_LONG mov_attachment(MKL_LONG target_particle, MKL_LONG given_chain_index, MKL_LONG flag_HEAD_TAIL)
  {
    /*
      flag_HEAD_TAIL is 0 for HEAD while 1 for TAIL.
    */
    ATTACHED(given_chain_index, flag_HEAD_TAIL) = target_particle;

    return 0;
  }

  MKL_LONG add_attachment(const MKL_LONG index_particle, MKL_LONG given_chain_index, MKL_LONG flag_HEAD_TAIL)
  {
    // there is no distingushable between add and mov attachment for chain node
    mov_attachment(index_particle, given_chain_index, flag_HEAD_TAIL);
    return 0;
  }

  MKL_LONG del_attachment(MKL_LONG given_chain_index, MKL_LONG flag_HEAD_TAIL)
  {
    // it make detachment for subjected chain ends
    // then attach to particles which other chain ends are attached
    mov_attachment(ATTACHED(given_chain_index, (flag_HEAD_TAIL + 1)%2), given_chain_index, flag_HEAD_TAIL);
    return 0;
  }
    
  MKL_LONG initial(MKL_LONG number_of_chains, MKL_LONG number_of_particles)
  {
    printf("ST:initial:CHAIN_INFORMATION\n");
    N_chains = number_of_chains;
    N_particles = number_of_particles;
    CHAIN = (CHAIN_NODE*)mkl_malloc(N_chains*sizeof(CHAIN_NODE), BIT);
    printf("end_dynamic_allocate\n");
    for(MKL_LONG i=0; i<N_chains; i++)
      {
        /*
          it allocate the equivalent number of chain ends per particles
        */
        HEAD(i) = i%N_particles;
        TAIL(i) = i%N_particles;
      }
    printf("end_particle_allocation\n");
    INITIALIZATION = TRUE;
    return INITIALIZATION;
  }

  MKL_LONG initial(ASSOCIATION& CONNECT)
  {
    printf("ST:initial:CHAIN_INFORMATION\n");
    N_chains = CONNECT.Nc*CONNECT.Np;
    N_particles = CONNECT.Np;
    CHAIN = (CHAIN_NODE*)mkl_malloc(N_chains*sizeof(CHAIN_NODE), BIT);
    printf("end_dynamic_allocate\n");
    MKL_LONG cnt = 0;
    for(MKL_LONG i=0; i<N_particles; i++)
      {
        for(MKL_LONG j=i; j<N_particles; j++) // it includes i=j which is loop and prevents j<i means all the pair is once accounted. 
          {
            MKL_LONG wij = CONNECT.return_multiplicity(i, j);
            if(i==j)
              {
                // note that weight count chain end rather than chain
                // when i != j, the duplication is removed because we only count the pair with condition j>i
                // when i == j, the duplication happens in the wij.
                wij /=2;
              }
            for(MKL_LONG w=0; w<wij; w++)
              {
                HEAD(cnt) = i;
                TAIL(cnt) = j;
                cnt++;
              }
          }
      }
    printf("end_particle_allocation\n");
    INITIALIZATION = TRUE;
    return INITIALIZATION;
  }

  MKL_LONG inheritance_chain(COND& given_condition, ASSOCIATION& CONNECT)
  {
    N_chains = CONNECT.Nc*CONNECT.Np;
    N_particles = CONNECT.Np;
    CHAIN = (CHAIN_NODE*)mkl_malloc(N_chains*sizeof(CHAIN_NODE), BIT);
    MKL_LONG cnt = 0;
    ifstream GIVEN_CHAIN_FILE;
    GIVEN_CHAIN_FILE.open(given_condition("CONTINUATION_CHAIN_FN").c_str());
    string line;
    while(getline(GIVEN_CHAIN_FILE, line))
      cnt ++;
    GIVEN_CHAIN_FILE.clear();
    GIVEN_CHAIN_FILE.seekg(0);
    for(MKL_LONG i=0; i<cnt-1; i++)
      getline(GIVEN_CHAIN_FILE, line);

    for(MKL_LONG k=0; k<N_chains; k++)
      GIVEN_CHAIN_FILE >> HEAD(k); // index from 0 to N_chains-1 : HEAD
    
    for(MKL_LONG k=0; k<N_chains; k++)
      GIVEN_CHAIN_FILE >> TAIL(k); // index from N_chains to 2*N_chains-1 : TAIL
    GIVEN_CHAIN_FILE.close();
    INITIALIZATION = TRUE;
    return INITIALIZATION;
  }
  
  CHAIN_INFORMATION() 
    {
      std::cout << "ERR: CHAIN Class must have initialization argument\n";
    }
  CHAIN_INFORMATION(MKL_LONG number_of_chains, MKL_LONG number_of_particles) 
    {
      initial(number_of_chains, number_of_particles);
    }
  CHAIN_INFORMATION(COND& given_condition)
    {
      if(given_condition("tracking_individual_chain") == "TRUE")
        {
          if(given_condition("CONTINUATION_CHAIN") == "TRUE")
            {
              printf("ERR: Without CONTINUATION TOPOLOGICAL INFORMATION (.hash, .weight), the tracking individual chain cannot inherit existing .chain file information\n");
              /* inheritance_chain(atoi(given_condition("N_chains_per_particle").c_str())*atoi(given_condition("Np").c_str()), atoi(given_condition("Np").c_str())); */
            }
          else
            {
              printf("ST:Constructor:CHAIN_INFORMATION\n");
              initial(atoi(given_condition("N_chains_per_particle").c_str())*atoi(given_condition("Np").c_str()), atoi(given_condition("Np").c_str()));
              printf("END:Constructor:CHAIN_INFORMATION\n");
            }
        }
    }
  CHAIN_INFORMATION(COND& given_condition, ASSOCIATION& CONNECT)
    // it initialize based on the CONNECT information 
    {
      if(given_condition("tracking_individual_chain") == "TRUE")
        {
          if(given_condition("CONTINUATION_CHAIN") == "TRUE")
            {
              inheritance_chain(given_condition, CONNECT);
            }
          else
            {
              printf("ST:Constructor:CHAIN_INFORMATION\n");
              initial(CONNECT);
              printf("END:Constructor:CHAIN_INFORMATION\n");
            }
        }
    }
  virtual ~CHAIN_INFORMATION()
    {
      if(INITIALIZATION)
        mkl_free(CHAIN);
    }
};


#endif
