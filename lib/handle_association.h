
#ifndef HANDLE_ASSOCIATION_H
#define HANDLE_ASSOCIATION_H

#include "association.h"
#include "geometry.h"
#include "potential.h"
#include "trajectory.h"

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

  MKL_LONG beads[4];
  /* // Commented for copy constructor validity. the references will be defined outside class object
  MKL_LONG &itself, &attached_bead, &new_attached_bead, &hash_attached_bead;
  */
  // action_arry setting for boosting up the boolean identifier
  // ACTION_ARR[IDX.CANCEL] act CANCEL. The IDX.CANCEL can be changed by IDX.ADD, IDX.DEL, IDX.MOV, respectively.
  MKL_LONG (*ACTION_ARR[4])(TRAJECTORY&, MKL_LONG, POTENTIAL_SET&, ASSOCIATION&, MKL_LONG[], MATRIX*);
  // Note that the last MATRIX* is related with R_minimum_distance_boost that has whole information
  
  /*
    The static const wokring on this way:
    const: will work as natural const
    static: will be shared with all the object from this class
   */
  static const MKL_LONG CANCEL = 0;
  static const MKL_LONG ADD = 1;
  static const MKL_LONG DEL = 2;
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
  MKL_LONG ACT(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX* R_minimum_distance_boost, MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX);  
  MKL_LONG UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
}

class CHAIN_HANDLE : public CHAIN_INFORMATION
{
  /*
    This class is designed to handle the CHAIN_INFORMATION with hash table.
    The role of this hash table is of importance to remove possible overhead during simulation, where the overhead introduced for searching the chain end attached to subjected particle.
    Hence, the hash table has information for direct map between particles and chain ends index.

    Notice that even if the class is using CHAIN_NODE which basically has both of chain head and tail, the hash table design to use chain end index (2*N_chains). Which means 0~N_chains-1 means the HEAD chain ends while N_chains ~ 2N_chains -1 means TAIL chain ends.
   */
 public:

  MKL_LONG **PARTICLE;
  MKL_LONG *P_TOKEN;
  MKL_LONG get_hash_index(MKL_LONG index_particle, MKL_LONG index_chain_end)
  {
    for(MKL_LONG i=0; i<P_TOKEN[index_particle]; i++)
      {
        if(PARTICLE[index_particle][i] == index_chain_end)
          return i;
      }
    std::cout<< "ERR: get_hash_index in CHAIN_HANDLE should be used for existing one\n";
    return -1;
  }

  MKL_LONG add_attachment(MKL_LONG index_particle, MKL_LONG index_chain_end) : add_attachment(index_particle, index_chain_end%N_chains, (MKL_LONG)(index_chain_end/N_chains))
    {
      PARTICLE[index_particle][P_TOKEN[index_particle]++] = index_chain_end;
      return 0;
    }

  MKL_LONG del_attachment(MKL_LONG index_chain_end) : del_attachment(index_chain_end%N_chains, (MKL_LONG)(index_chain_end/N_chains))
    {
      for(MKL_LONG i=get_hash_index(index_particle); i<P_TOKEN[index_particle]-1; i++)
        {
          PARTICLE[index_particle][token] = PARTICLE[index_particle][token+1];
        }
      PARTICLE[index_particle][P_TOKEN[index_particle]] = -1;
      P_TOKEN[index_particle]--;
      return 0;
    }
  
  /* MKL_LONG particle_add_CE(MKL_LONG index_particle, MKL_LONG index_chain_end) */
  /* { */
  /*   PARTICLE[index_particle][P_TOKEN[index_particle]++] = index_chain_end; */
  /* } */
  
  /* MKL_LONG particle_del_CE(MKL_LONG index_particle, MKL_LONG index_chain_end) */
  /* { */
  /*   /\* MKL_LONG token = 0; *\/ */
  /*   /\* while(PARTICLE[index_particle][token++]!=index_chain_end){} *\/ */
  /*   for(MKL_LONG i=get_hash_index(index_particle); i<P_TOKEN[index_particle]-1; i++) */
  /*     { */
  /*       PARTICLE[index_particle][token] = PARTICLE[index_particle][token+1]; */
  /*     } */
  /*   PARTICLE[index_particle][P_TOKEN[index_particle]] = -1; */
  /*   P_TOKEN[index_particle]--; */
  /*   return 0; */
  /* } */

  MKL_LONG& CE_ATTAHCED(MKL_LONG given_chain_end_index)
    {
      return ATTACHED(given_chain_end_index%N_chains, (MKL_LONG)(given_chain_end_index/N_chains));
    }
  
  MKL_LONG mov_attachment(MKL_LONG index_particle, MKL_LONG given_chain_end_index) : mov_attachment(index_particle, given_chain_end_index%N_chains, (MKL_LONG)(given_chain_end_index/N_chains))
    {
      del_attachment(given_chain_end_index);
      add_attachment(index_particle, given_chain_end_index);
      /* // deleting existing bridge */
      /* particle_del_CE(CE_ATTACHED(given_chain_end_index), given_chain_end_index); */
      /* // adding new bridge */
      /* particle_add_CE(index_particle, given_chain_end); */
      return 0;
    }

  MKL_LONG get_index_degeneracy(MKL_LONG index_subject_chain_end, MKL_LONG index_opposite_end)
  {
    // working in progress
  }
    
  TRACKING_ACTION(ASSOCIATION& CONNECT, MKL_LONG flag_ACTION, MKL_LONG *index_set)
    {
      /*
        It is of importance to using ij indices by identify subjected chain end.
        The basic mechanism in association, handle_association library detaches opponents of subjected chain ends because of convenience. 
        Therefore, the index_subjected_chain_end is using flag_itself and flag_other but the define should be inside sequence of if-phrase.
        It might be changed for future convenience, then it would be changed in the Boltzmann distribution in stochastic_simulation.cpp, ACTION::ADD and ACTION::MOV in handle_association.cpp, del_association and mov_association in association.cpp.
        At this moment, keep the conventional rule even if it is in-convinience with our logic.
        Keep in mind that there is no distingushable for the transition probability of both chain ends belong one chain.
       */
      MKL_LONG index_subject_chain_end = -1;
      if(flag_ACTION == CONNECT.CANCEL)
        return 0;
      else if(flag_ACTION == CONNECT.ADD)
        {
          index_subject_chain_end = get_index_degeneracy(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_other]);
          add_attachment(index_set[CONNECT.flag_new], index_subject_chain_end)
        }
      else if(flag_ACTION == CONNECT.DEL)
        {
          index_subject_chain_end = get_index_degeneracy(index_set[CONNECT.flag_other], index_set[CONNECT.flag_itself]);
          del_attachment(index_subject_chain_end)
        }
      else if(flag_ACTION == CONNECT.MOV)
        {
          /* mov_attachment(index_set[CONNECT.flag_new], index_subject_chain_end); */
          index_subject_chain_end = get_index_degeneracy(index_set[CONNECT.flag_other], index_set[CONNECT.flag_itself]);
          del_attachment(index_subject_chain_end);
          add_attachment(index_set[CONNECT.flag_new], index_subject_chain_end);
        }
      else
        {
          std::cout << "ERR: No left option\n";
        }
      return 0;
    }
  
  MKL_LONG hash_initial()
  {
    PARTICLE = (MKL_LONG**) mkl_malloc(N_particles*sizeof(MKL_LONG*), BIT);
    P_TOKEN = (MKL_LONG*) mkl_malloc(N_particles*sizeof(MKL_LONG), BIT);
    for(MKL_LONG i=0; i<N_particles; i++)
      {
        /* 
           The PARTICLE is N_particles by 2*N_chains matrix.
           For default test (1000 particles with 10 chains per particles), the array might have 160 Mb, which is quite large but not too much.
         */
        PARTICLE[i] = (MKL_LONG*) mkl_malloc(2*N_chains*sizeof(MKL_LONG), BIT);
        for(MKL_LONG j=0; j<2*N_chains; j++)
          {
            PARTICLE[i][j] = -1; // this means the hash direct no chain ends (index for chain ends is started with 0)
          }
        P_TOKEN[i] = 0;
      }
    for(MKL_LONG i=0; i<2*N_chains; i++)
      {
        /* CE_ATTACHED(i) */
        particle_add_CE(CE_ATTACHED(i), i);
      }
  }
 CHAIN_HANDLE(MKL_LONG number_of_chains, MKL_LONG number_of_particles) : CHAIN_INFORMATION(number_of_chains, number_of_particles)
    {
      hash_initial();
    }
 CHAIN_HANDLE(COND& given_condition) : CHAIN_INFORMATION(given_condition)
  {
    hash_initial();
  }
  ~CHAIN_HANDLE()
    {
      if(INITIALIZATION)
        {
          for(MKL_LONG i=0; i<N_particles; i++)
            {
              mkl_free(PARTICLES[i]);
            }
          mkl_free(PARTICLES);
          mkl_free(P_TOKEN);
        }
    }
}

#endif
