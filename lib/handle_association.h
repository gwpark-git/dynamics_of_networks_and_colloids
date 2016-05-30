
#ifndef HANDLE_ASSOCIATION_H
#define HANDLE_ASSOCIATION_H

#define TRUE 1
#define FALSE 0

#include "association.h"
#include "geometry.h"
#include "potential.h"
/* #include "trajectory.h" */ // the trajectory dependecy is removed since the R_boost is applied
#include "random.h"

// this library is designed to handle association in easiler way.
// the reason to seperate it from lib_association is for removing dependencies of library

class INDEX_MC
{
  // index set for MC steps
 public:
  // INDEX_MC::beads information include the all related bead indice.
  // The index value is related with the flag in the lib_association, which is not syncronized on this sourcecode at the moment.
  // Shortly, index_itself refers beads[2] this is because boolean optimization and handling indce easiler ways, and beads[0] is refer index_other and beads[1] is index_new. beads[3] is added in order to handle hash for index_other.
  // For better understanding, the reference variables are tagged and initialized in the constructor. For details, see the constructor in the sourcecode.

  MKL_LONG beads[5];
  /* MKL_LONG &bead_selected_chain_end, &bead_opp_selected_chain_end, &bead_new_opp_selected_chain_end, &hash_opp_selected_chain_end, &hash_backtrace; // this will be initialize when constructor have been called */
  
  /* // Commented for copy constructor validity. the references will be defined outside class object
     MKL_LONG &itself, &attached_bead, &new_attached_bead, &hash_attached_bead, &hash_backtrace;
  */
  // action_arry setting for boosting up the boolean identifier
  // ACTION_ARR[IDX.CANCEL] act CANCEL. The IDX.CANCEL can be changed by IDX.ADD, IDX.OPP_DEL, IDX.MOV, respectively.
  MKL_LONG (*ACTION_ARR[4])(POTENTIAL_SET&, ASSOCIATION&, MKL_LONG[], MATRIX*);

  // Note that the last MATRIX* is related with R_minimum_distance_boost that has whole information
  
  /*
    The static const wokring on this way:
    const: will work as natural const
    static: will be shared with all the object from this class
  */
  static const MKL_LONG CANCEL = 0;
  static const MKL_LONG ADD = 1;
  static const MKL_LONG OPP_DEL = 2;
  static const MKL_LONG MOV = 3;
  static const MKL_LONG LOCK = 4;
  static const MKL_LONG N_BOOST_COUNT[]; // defined on sourcefile
  INDEX_MC();
  ~INDEX_MC(){}

  /* 
     This copy constructor is of importance to prevent the potential problem.
     Remind rule of three.
  */ 
  MKL_LONG copy_INDEX(const INDEX_MC& given_IDX);
  INDEX_MC& operator=(const INDEX_MC& given_IDX);
  INDEX_MC(const INDEX_MC& given_IDX); // copy constructore
  
  MKL_LONG set_initial_variables();
  MKL_LONG initial();
  
};

namespace SEARCHING
{
  //the basic backward searching method. Note that the probability map is given by sorted manner (this is due to the fact that the proability map is computed once a each time steps.) Therefore, if the cut-off scheme for FENE connector is properly working, there is no problem for the backtrace.
  MKL_LONG backtrace(MATRIX& given_arr, double p);
  MKL_LONG backtrace_cell_list(MATRIX& given_arr, MKL_LONG TOKEN, double p, MKL_LONG index_particle, RDIST& R_boost);
  
  // For general purpose, the bisection method is good aspect for approaching. Since the domain of this bisection is index space, that is defined on Integer space, the steps to reach solution is log_2(Np). In this scheme, however, some additional step is needed in order to compensate the loss of odd number division. Therefore, the overall step is determined log_2(Np) + #compensation, which approximately 11 in the Np640 cases (log2(640) is 8).
  MKL_LONG bisection(MATRIX& given_arr, double p);
}

namespace ACTION
{
  MKL_LONG ACT(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX* R_minimum_distance_boost, MKL_LONG const IDENTIFIER_ACTION);
  
  MKL_LONG IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX);  
  double UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION);

  MKL_LONG CANCEL(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG MOV(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG OPP_DEL(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG ADD(POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  
}

class CHAIN_HANDLE : public CHAIN_INFORMATION
{
  /*
    This class is designed to handle the CHAIN_INFORMATION with hash table.
    The role of this hash table is of importance to remove possible overhead during simulation, where the overhead introduced for searching the chain end attached to subjected particle.
    Hence, the hash table has information for direct map between particles and chain ends index.

    Notice that even if the class is using CHAIN_NODE which basically has both of chain head and tail, the hash table design to use chain end index (2*N_chains). Which means 0~N_chains-1 means the HEAD chain ends while N_chains ~ 2N_chains -1 means TAIL chain ends.

    For performance issue, almost of all member functions are defined inside class declaration in order to make inline function. These functions are not so long, which means the function call takes relatively high time compared with the overall time to excute the function.
  */
 public:

  MKL_LONG **PARTICLE;
  MKL_LONG *P_TOKEN;
  gsl_rng *r_degeneracy_check;
  MKL_LONG *degeneracy_index_array;
  MKL_LONG& CE_ATTACHED_REF(MKL_LONG given_chain_end_index)
    {
      return ATTACHED(given_chain_end_index%N_chains, (MKL_LONG)(given_chain_end_index/N_chains));
    }
  
  MKL_LONG opp_chain_end_index(MKL_LONG given_chain_end_index)
  {
    if (given_chain_end_index < N_chains)
      return given_chain_end_index + N_chains;
    return given_chain_end_index - N_chains;
  }
  
  MKL_LONG get_hash_index(MKL_LONG index_particle, MKL_LONG index_chain_end)
  {
    for(MKL_LONG i=0; i<P_TOKEN[index_particle]; i++)
      {
        if(PARTICLE[index_particle][i] == index_chain_end)
          return i;
      }
    printf("ERR: get_hash_index in CHAIN_HANDLE should be used for existing one. index for particle and chain end: (%ld, %ld)\n", index_particle, index_chain_end);
    return -1;
  }

  MKL_LONG mov_attachment(MKL_LONG target_particle, MKL_LONG given_chain_end_index);
  MKL_LONG get_index_degeneracy(MKL_LONG particle_subject, MKL_LONG particle_target);
  MKL_LONG TRACKING_ACTION(ASSOCIATION& CONNECT, MKL_LONG flag_ACTION, INDEX_MC& IDX);


  MKL_LONG write(std::ofstream& file);
  MKL_LONG allocate_array();
  MKL_LONG hash_initial(MKL_LONG seed);
 CHAIN_HANDLE() : CHAIN_INFORMATION(){}
 CHAIN_HANDLE(MKL_LONG number_of_chains, MKL_LONG number_of_particles) : CHAIN_INFORMATION(number_of_chains, number_of_particles)
    {
      /* printf("CHECK\n"); */
      hash_initial(random());
    }
 CHAIN_HANDLE(COND& given_condition) : CHAIN_INFORMATION(given_condition)
  {
    /* printf("CHECK\n"); */
    if(given_condition("tracking_individual_chain") == "TRUE")
      {
        hash_initial(atoi(given_condition("basic_random_seed_SS").c_str()));
      }
    else
      {
        INITIALIZATION = FALSE;
      }
  }
 CHAIN_HANDLE(COND& given_condition, ASSOCIATION& CONNECT) : CHAIN_INFORMATION(given_condition, CONNECT)
  {
    if(given_condition("tracking_individual_chain") == "TRUE")
      {
        /* if(given_condition("CONTINUATION_CHAIN") == "TRUE") */
        /*   { */
        /*     hash_initial_ */
        /*   } */
        /* else */
        /*   { */
            hash_initial(atoi(given_condition("basic_random_seed_SS").c_str()));
          /* } */
      }
  }
  
  ~CHAIN_HANDLE()
    {
      if(INITIALIZATION)
        {
          for(MKL_LONG i=0; i<N_particles; i++)
            {
              mkl_free(PARTICLE[i]);
            }
          mkl_free(PARTICLE);
          mkl_free(P_TOKEN);
          gsl_rng_free(r_degeneracy_check);
          mkl_free(degeneracy_index_array);
        }
    }
  MKL_LONG allocate_existing_bridges(ASSOCIATION& CONNECT);

};

// inline member function definition



#endif
