
#ifndef HANDLE_ASSOCIATION_H
#define HANDLE_ASSOCIATION_H

#include "association.h"
#include "geometry.h"
#include "potential.h"
#include "trajectory.h"
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

  MKL_LONG beads[4];
  /* // Commented for copy constructor validity. the references will be defined outside class object
  MKL_LONG &itself, &attached_bead, &new_attached_bead, &hash_attached_bead;
  */
  // action_arry setting for boosting up the boolean identifier
  // ACTION_ARR[IDX.CANCEL] act CANCEL. The IDX.CANCEL can be changed by IDX.ADD, IDX.OPP_DEL, IDX.MOV, respectively.
  MKL_LONG (*ACTION_ARR[4])(TRAJECTORY&, MKL_LONG, POTENTIAL_SET&, ASSOCIATION&, MKL_LONG[], MATRIX*);
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
  MKL_LONG ACT(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX* R_minimum_distance_boost, MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX);  
  MKL_LONG UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION);
  MKL_LONG CANCEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG MOV(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG OPP_DEL(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
  MKL_LONG ADD(TRAJECTORY& TRAJ, MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost);
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
  /* MKL_LONG CE_ATTACHED(MKL_LONG given_chain_end_index) */
  /* { */
    
  /* } */
  
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

  MKL_LONG add_attachment(MKL_LONG index_particle, MKL_LONG index_chain_end)
    {
      CHAIN_INFORMATION::add_attachment(index_particle, index_chain_end%N_chains, (MKL_LONG)(index_chain_end/N_chains)); // it calls existing function in CHAIN_INFORMATION
      PARTICLE[index_particle][P_TOKEN[index_particle]++] = index_chain_end;
      return 0;
    }

  MKL_LONG del_attachment(MKL_LONG index_chain_end)
    {
      MKL_LONG index_particle = (MKL_LONG)CE_ATTACHED_REF(index_chain_end);
      CHAIN_INFORMATION::del_attachment(index_chain_end%N_chains, (MKL_LONG)(index_chain_end/N_chains)); // it calls existing function in CHAIN_INFORMATION
      for(MKL_LONG i=get_hash_index(index_particle, index_chain_end); i<P_TOKEN[index_particle]-1; i++)
        {
          PARTICLE[index_particle][i] = PARTICLE[index_particle][i+1];
        }
      PARTICLE[index_particle][P_TOKEN[index_particle]] = -1;
      P_TOKEN[index_particle]--;
      return 0;
    }

  MKL_LONG mov_attachment(MKL_LONG index_particle, MKL_LONG given_chain_end_index)
    {
      CHAIN_INFORMATION::mov_attachment(index_particle, given_chain_end_index%N_chains, (MKL_LONG)(given_chain_end_index/N_chains)); // it calls existing function in CHAIN_INFORMATION
      del_attachment(given_chain_end_index);
      add_attachment(index_particle, given_chain_end_index);
      /* // deleting existing bridge */
      /* particle_del_CE(CE_ATTACHED_REF(given_chain_end_index), given_chain_end_index); */
      /* // adding new bridge */
      /* particle_add_CE(index_particle, given_chain_end); */
      return 0;
    }
  
  MKL_LONG particle_add_CE(MKL_LONG index_particle, MKL_LONG index_chain_end)
  {
    PARTICLE[index_particle][P_TOKEN[index_particle]++] = index_chain_end;
    return 0;
  }
  
  MKL_LONG particle_del_CE(MKL_LONG index_particle, MKL_LONG index_chain_end)
  {
    /* MKL_LONG token = 0; */
    /* while(PARTICLE[index_particle][token++]!=index_chain_end){} */
    for(MKL_LONG i=get_hash_index(index_particle, index_chain_end); i<P_TOKEN[index_particle]-1; i++)
      {
        PARTICLE[index_particle][i] = PARTICLE[index_particle][i+1];
      }
    PARTICLE[index_particle][P_TOKEN[index_particle]] = -1;
    P_TOKEN[index_particle]--;
    return 0;
  }


  MKL_LONG get_index_degeneracy(MKL_LONG particle_subject, MKL_LONG particle_target)
  {
    /*
      This index function is the core for the tracking algorithm.
      Basically, the tracking individual chain related with the data structure.
      The used adjacence list in our algorithm, it is not distingushable between different chain with the same pair of attachment, but counting multiple connection individually.
      That is one of the reason to reduce computational time dramatically from the original code since we can extracting out the computational overhead for degeneracy.
      To track individual chain, the degeneracy is matter which means we need one more selecting procedure in the degenerated chains.
      This algorithm is selecting one of the degenerated chain using generated random number index.
    */
    MKL_LONG degeneracy=0;
    // check the chain: HEAD - particle_subject, TAIL - particle_target
    for(MKL_LONG i=0; i<P_TOKEN[particle_subject]; i++)
      {
        if(PARTICLE[particle_subject][i] == particle_target)
          degeneracy_index_array[degeneracy++] = PARTICLE[particle_subject][i];
      }
    MKL_LONG count_head_degeneracy = degeneracy;
    // check the chain: HEAD - particle_target, TAIL - particle_subject
    for(MKL_LONG i=0; i<P_TOKEN[particle_target]; i++)
      {
        if(PARTICLE[particle_target][i] == particle_subject)
          degeneracy_index_array[degeneracy++] = PARTICLE[particle_target][i];
      }
    MKL_LONG selecting_index = gsl_rng_uniform_int(r_degeneracy_check, degeneracy);
    MKL_LONG subjected_chain_end = degeneracy_index_array[selecting_index];
    if (selecting_index >= count_head_degeneracy) // when the subjected chain end attached on subjected particle in TAIL
      subjected_chain_end = subjected_chain_end%N_chains;
    return subjected_chain_end;
  }
    
  MKL_LONG TRACKING_ACTION(ASSOCIATION& CONNECT, MKL_LONG flag_ACTION, INDEX_MC& IDX)
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
      if(flag_ACTION == INDEX_MC::CANCEL)
        return 0;
      else if(flag_ACTION == INDEX_MC::ADD)
        {
          // head = subjected chain end, head = other chain end
          index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_other]);
          add_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end);
        }
      else if(flag_ACTION == INDEX_MC::OPP_DEL)
        {
          // head = other chain end, tail = subjected chain end
          index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_itself]);
          // note that it calls the typical 'del' function rather than 'opp_del'
          // since the given index already transposed.
          del_attachment(index_subject_chain_end);
        }
      else if(flag_ACTION == INDEX_MC::MOV)
        {
          /* mov_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end); */
          // the following delete scheme is exactly the same for OPP_DEL functionality.
          index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_itself]);
          del_attachment(index_subject_chain_end);
          // Because of the previous detachment (way back to the origin), both of chain ends are attached to the same particle.
          // Therefore, selection HEAD or TAIL is NOT importance to make new bridge.
          // Note that HEAD and TAIL are not related with the degeneracy, which means the checking degeneracy function 'get_index_degeneracy' checked for both pairs for HEAD-TAIL and TAIL-HEAD.
          add_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end);
        }
      else
        {
          std::cout << "ERR: No left option in the TRACKING_ACTION function in handle_association.h \n";
        }
      return 0;
    }
  
  MKL_LONG hash_initial();
 CHAIN_HANDLE() : CHAIN_INFORMATION(){}
 CHAIN_HANDLE(MKL_LONG number_of_chains, MKL_LONG number_of_particles) : CHAIN_INFORMATION(number_of_chains, number_of_particles)
    {
      /* printf("CHECK\n"); */
      hash_initial();
    }
 CHAIN_HANDLE(COND& given_condition) : CHAIN_INFORMATION(given_condition)
  {
    /* printf("CHECK\n"); */
    hash_initial();
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


#endif
