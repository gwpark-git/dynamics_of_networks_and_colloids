
#ifndef HANDLE_ASSOCIATION_H
#define HANDLE_ASSOCIATION_H

#define TRUE 1
#define FALSE 0

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

  MKL_LONG mov_attachment(MKL_LONG target_particle, MKL_LONG given_chain_end_index)
    {
      /* del_attachment(given_chain_end_index); */
      /* add_attachment(target_particle, given_chain_end_index); */
      MKL_LONG index_particle = (MKL_LONG)CE_ATTACHED_REF(given_chain_end_index);
      /* MKL_LONG index_hash_opp = get_hash_index(index_particle, given_chain_end_index); */
      /* MKL_LONG opp_particle = PARTICLE[index_particle][index_hash_opp]; */

      // delete existing information
      MKL_LONG hash_check = get_hash_index(index_particle, given_chain_end_index);
      if(hash_check == -1)
	{
	  printf("given chain end index: %ld: attached to HEAD and TAIL: (%ld, %ld) and index particle = %ld\n", given_chain_end_index, HEAD(given_chain_end_index%N_chains), TAIL(given_chain_end_index%N_chains), index_particle);
	  for(MKL_LONG i=0; i<N_particles; i++)
	    {
	      for(MKL_LONG j=0; j<P_TOKEN[i]; j++)
		{
		  if(PARTICLE[i][j] == given_chain_end_index)
		    printf("found subjected chain end in PARTICLE[%ld][%ld] = %ld\n", i, j, given_chain_end_index);
		}
	    }
	  printf("attached particle: %ld\n", CE_ATTACHED_REF(given_chain_end_index));
	  for(MKL_LONG i=0; i<P_TOKEN[index_particle] + 10; i++)
	    {
	      printf("[%ld, %ld] = %ld\t [%ld, %ld] = %ld\n", HEAD(given_chain_end_index%N_chains), i, PARTICLE[HEAD(given_chain_end_index%N_chains)][i], TAIL(given_chain_end_index%N_chains), i, PARTICLE[TAIL(given_chain_end_index%N_chains)][i]);
	    }
	}
      for(MKL_LONG i=get_hash_index(index_particle, given_chain_end_index); i<P_TOKEN[index_particle]-1; i++)
	{
	  PARTICLE[index_particle][i] = PARTICLE[index_particle][i+1];
	}
      PARTICLE[index_particle][P_TOKEN[index_particle]--] = -1;

      // add new information
      /* if(given_chain_end_index == 1691) */
      /* 	{ */
      /* 	  printf("HISTORY for CE = 1691: ATTACHED to %ld\n", target_particle); */
      /* 	} */
      PARTICLE[target_particle][P_TOKEN[target_particle]++] = given_chain_end_index;
      /* MKL_LONG opp_chain_end = opp_chain_end_index(given_chain_end_index); */
      /* PARTICLE[opp_particle][get_hash_index(opp_particle, opp_chain_end_index(given_chain_end_index))] = target_particle; */
      // it just change from the existing one to the new one. Therefore, the P_TOKEN will not be changed
      CHAIN_INFORMATION::mov_attachment(target_particle, given_chain_end_index%N_chains, (MKL_LONG)(given_chain_end_index/N_chains)); // it calls existing function in CHAIN_INFORMATION

      return 0;
    }
  /* MKL_LONG add_attachment(MKL_LONG target_particle, MKL_LONG given_chain_end_index) */
  /*   { */
  /*     /\* // from loop to bridge *\/ */
  /*     /\* CHAIN_INFORMATION::add_attachment(index_particle, index_chain_end%N_chains, (MKL_LONG)(index_chain_end/N_chains)); // it calls existing function in CHAIN_INFORMATION *\/ */
  /*     /\* PARTICLE[index_particle][P_TOKEN[index_particle]++] = index_chain_end; *\/ */
  /*     /\* return 0; *\/ */
  /*     mov_attachment(target_particle, given_chain_end_index); */
  /*     return 0; */
  /*   } */

  /* MKL_LONG del_attachment(MKL_LONG index_chain_end) */
  /*   { */
  /*     // from bridge to loop */
  /*     MKL_LONG index_particle = (MKL_LONG)CE_ATTACHED_REF(index_chain_end); */
  /*     MKL_LONG index_hash_opp = get_hash_index(index_particle, index_chain_end); */
  /*     MKL_LONG opp_particle = PARTICLE[index_particle][index_hash_opp]; */
  /*     CHAIN_INFORMATION::del_attachment(index_chain_end%N_chains, (MKL_LONG)(index_chain_end/N_chains)); // it calls existing function in CHAIN_INFORMATION */
  /*     // following deletes particle information from the index_particle */
  /*     for(MKL_LONG i=index_hash_opp; i<P_TOKEN[index_particle]-1; i++) */
  /*       { */
  /*         PARTICLE[index_particle][i] = PARTICLE[index_particle][i+1]; */
  /*       } */
  /*     PARTICLE[index_particle][P_TOKEN[index_particle]] = -1; */
  /*     P_TOKEN[index_particle]--; */

  /*     // following adding new information for forming bridge chain */
  /*     PARTICLE[opp_particle][P_TOKEN[opp_particle]++] = index_chain_end; */
  /*     return 0; */
  /*   } */

  
  
  /* MKL_LONG particle_add_CE(MKL_LONG index_particle, MKL_LONG index_chain_end) */
  /* { */
  /*   PARTICLE[index_particle][P_TOKEN[index_particle]++] = index_chain_end; */
  /*   return 0; */
  /* } */
  
  /* MKL_LONG particle_del_CE(MKL_LONG index_particle, MKL_LONG index_chain_end) */
  /* { */
  /*   /\* MKL_LONG token = 0; *\/ */
  /*   /\* while(PARTICLE[index_particle][token++]!=index_chain_end){} *\/ */
  /*   for(MKL_LONG i=get_hash_index(index_particle, index_chain_end); i<P_TOKEN[index_particle]-1; i++) */
  /*     { */
  /*       PARTICLE[index_particle][i] = PARTICLE[index_particle][i+1]; */
  /*     } */
  /*   PARTICLE[index_particle][P_TOKEN[index_particle]] = -1; */
  /*   P_TOKEN[index_particle]--; */
  /*   return 0; */
  /* } */


  MKL_LONG get_index_degeneracy(MKL_LONG particle_subject, MKL_LONG particle_target)
  {
    /*
      This index function is the core for the tracking algorithm.
      Basically, the tracking individual chain related with the data structure.
      The used adjacence list in our algorithm, it is not distingushable between different chain with the same pair of attachment, but counting multiple connection individually.
      That is one of the reason to reduce computational time dramatically from the original code since we can extracting out the computational overhead for degeneracy.
      To track individual chain, the degeneracy is matter which means we need one more selecting procedure in the degenerated chains.
      This algorithm is selecting one of the degenerated chain using generated random number index.

      It is of importance that the returning chain end index should be related with particle_subjection because the flag controls in the TRACKING_ACTION function is ordered function.
    */
    MKL_LONG degeneracy=0;
    // check the chain: HEAD - particle_subject, TAIL - particle_target
    // note that it include all the chain end attached to particle_subject
    for(MKL_LONG i=0; i<P_TOKEN[particle_subject]; i++)
      {
        if((MKL_LONG)CE_ATTACHED_REF(opp_chain_end_index(PARTICLE[particle_subject][i])) == particle_target)
	  // it check the bridges. Therefore, opp_chain_end_index is of importance since PARTICLE[particle_subject][i] is the chain end attached to particle_subject
          degeneracy_index_array[degeneracy++] = PARTICLE[particle_subject][i];
      }
    /*
      following reverse map is not necessary
     */
    /*     MKL_LONG count_head_degeneracy = degeneracy;
    /* // check the chain: HEAD - particle_target, TAIL - particle_subject */
    /* for(MKL_LONG i=0; i<P_TOKEN[particle_target]; i++) */
    /*   { */
    /*     if((MKL_LONG)CE_ATTACHED_REF(PARTICLE[particle_target][i]) == particle_subject) */
    /*       degeneracy_index_array[degeneracy++] = PARTICLE[particle_target][i]; */
    /*   } */
    /* printf("N_DEGENERACY = %ld\n", degeneracy); */
    if (degeneracy < 1)
      {
	printf("ERR: degeneracy is %ld, which is not valid for get_index_degeneracy\n", degeneracy);
	return -1;
      }
    MKL_LONG selecting_index = gsl_rng_uniform_int(r_degeneracy_check, degeneracy);
    MKL_LONG subjected_chain_end = degeneracy_index_array[selecting_index];
    /* if (selecting_index >= count_head_degeneracy) // when the subjected chain end attached on subjected particle in TAIL */
    /*   subjected_chain_end = subjected_chain_end%N_chains; */
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
      if(flag_ACTION == INDEX_MC::CANCEL)
        return 0;
      else
	{
	  MKL_LONG index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_itself]);
	  if(index_subject_chain_end == -1)
	    {
	      printf("particle info: (itself, other, new) = (%ld, %ld, %ld)\n", IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_new]);
	      for(MKL_LONG i=0; i<P_TOKEN[IDX.beads[CONNECT.flag_other]]; i++)
		{
		  printf("(C, OC)([%ld, %ld] = %ld) = (%ld, %ld)\n", IDX.beads[CONNECT.flag_other], i, PARTICLE[IDX.beads[CONNECT.flag_other]][i], (MKL_LONG)CE_ATTACHED_REF(PARTICLE[IDX.beads[CONNECT.flag_other]][i]), (MKL_LONG)CE_ATTACHED_REF(opp_chain_end_index(PARTICLE[IDX.beads[CONNECT.flag_other]][i])));
		}
	    }
	  
	  if(flag_ACTION == INDEX_MC::ADD)
	    {
	      // for loop chains
	      // head = subjected chain end, head = other chain end
	      /* index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_other]); */
	      /* add_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end); */
	      mov_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end);
	    }
	  else if(flag_ACTION == INDEX_MC::OPP_DEL)
	    {
	      // for bridge chains
	      // head = other chain end, tail = subjected chain end
	      /* index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_other]); */
	      /* index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_itself]); */
	      // note that it calls the typical 'del' function rather than 'opp_del'
	      // the subject of dell is the first argument, which means the first argument will be deleted
	  
	      mov_attachment(IDX.beads[CONNECT.flag_itself], index_subject_chain_end);
	    }
	  else if(flag_ACTION == INDEX_MC::MOV)
	    {
	      // for bridge chains
	      // the following delete scheme is exactly the same for OPP_DEL functionality.
	      /* index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_itself], IDX.beads[CONNECT.flag_other]); */
	      /* index_subject_chain_end = get_index_degeneracy(IDX.beads[CONNECT.flag_other], IDX.beads[CONNECT.flag_itself]); */
	      mov_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end);
	      /* del_attachment(index_subject_chain_end); */
	      // Because of the previous detachment (way back to the origin), both of chain ends are attached to the same particle.
	      // Therefore, selection HEAD or TAIL is NOT importance to make new bridge.
	      // Note that HEAD and TAIL are not related with the degeneracy, which means the checking degeneracy function 'get_index_degeneracy' checked for both pairs for HEAD-TAIL and TAIL-HEAD.
	      /* add_attachment(IDX.beads[CONNECT.flag_new], index_subject_chain_end); */
	    }
	  else
	    {
	      std::cout << "ERR: No left option in the TRACKING_ACTION function in handle_association.h \n";
	    }
	}
      return 0;
    }

  MKL_LONG write(const char *fn)
  {
    std::ofstream FILE1;
    FILE1.open(fn, std::ios_base::app);
    if(!FILE1)
      {
	std::cout << "ERROR OCCURS DURING WRITING CHAIN INFORMATION\n";
	return -1;
      }
    MKL_LONG cnt = 0;
    for(MKL_LONG i=0; i<2*N_chains; i++)
      {
	FILE1 << CE_ATTACHED_REF(i) << '\t';
      }
    FILE1 << std::endl;
    FILE1.close();
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
    if(given_condition("tracking_individual_chain") == "TRUE")
      {
        hash_initial();
      }
    else
      {
        INITIALIZATION = FALSE;
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


#endif
