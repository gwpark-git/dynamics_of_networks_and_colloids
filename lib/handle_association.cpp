

#include "handle_association.h"

MKL_LONG CHAIN_HDNALE::write(std::ofstream& file)
{
  MKL_LONG cnt = 0;
  for(MKL_LONG i=0; i<2*N_chains; i++)
    {
      file << CE_ATTACHED_REF(i) << '\t';
    }
  file << std::endl;
  /* FILE1.close(); */
  return 0;
}


MKL_LONG CHAIN_HANDLE::mov_attachment(MKL_LONG target_particle, MKL_LONG given_chain_end_index)
{
  /* del_attachment(given_chain_end_index); */
  /* add_attachment(target_particle, given_chain_end_index); */
  MKL_LONG index_particle = (MKL_LONG)CE_ATTACHED_REF(given_chain_end_index);
  /* MKL_LONG index_hash_opp = get_hash_index(index_particle, given_chain_end_index); */
  /* MKL_LONG opp_particle = PARTICLE[index_particle][index_hash_opp]; */

  // delete existing information
  MKL_LONG hash_check = get_hash_index(index_particle, given_chain_end_index);
  /* if(hash_check == -1) */
  /*   { */
  /*     printf("given chain end index: %ld: attached to HEAD and TAIL: (%ld, %ld) and index particle = %ld\n", given_chain_end_index, HEAD(given_chain_end_index%N_chains), TAIL(given_chain_end_index%N_chains), index_particle); */
  /*     for(MKL_LONG i=0; i<N_particles; i++) */
  /*       { */
  /*         for(MKL_LONG j=0; j<P_TOKEN[i]; j++) */
  /*           { */
  /*             if(PARTICLE[i][j] == given_chain_end_index) */
  /*               printf("found subjected chain end in PARTICLE[%ld][%ld] = %ld\n", i, j, given_chain_end_index); */
  /*           } */
  /*       } */
  /*     printf("attached particle: %ld\n", CE_ATTACHED_REF(given_chain_end_index)); */
  /*     for(MKL_LONG i=0; i<P_TOKEN[index_particle] + 10; i++) */
  /*       { */
  /*         printf("[%ld, %ld] = %ld\t [%ld, %ld] = %ld\n", HEAD(given_chain_end_index%N_chains), i, PARTICLE[HEAD(given_chain_end_index%N_chains)][i], TAIL(given_chain_end_index%N_chains), i, PARTICLE[TAIL(given_chain_end_index%N_chains)][i]); */
  /*       } */
  /*   } */
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


MKL_LONG CHAIN_HANDLE::TRACKING_ACTION(ASSOCIATION& CONNECT, MKL_LONG flag_ACTION, INDEX_MC& IDX)
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

MKL_LONG CHAIN_HANDLE::get_index_degeneracy(MKL_LONG particle_subject, MKL_LONG particle_target)
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

const MKL_LONG INDEX_MC::N_BOOST_COUNT[] = {0, 2, 2, 3};
/*
  This is defined outside of class definition in the header file.
  The initialization for MKL_LONG CACEL and so on are fine with static constant definition.
  However, in the case of array, it should be defined on here.
  N_BOOST_COUNT[IDX.CANCEL] = 0;
  N_BOOST_COUNT[IDX.ADD] = 2;
  N_BOOST_COUNT[IDX.OPP_DEL] = 2;
  N_BOOST_COUNT[IDX.MOV] = 3;
 */

/*
  the reason for inheritance for reference variable is that the reference variables cannot be along without given variables. In this way, the reference variables automatically allocated with exiting variables with beads array.
*/
INDEX_MC::INDEX_MC()
{
  initial();
  set_initial_variables();
}

INDEX_MC::INDEX_MC(const INDEX_MC& given_IDX)
{
  copy_INDEX(given_IDX);
}

MKL_LONG INDEX_MC::initial()
{
  ACTION_ARR[CANCEL] = ACTION::CANCEL ;
  ACTION_ARR[ADD]    = ACTION::ADD    ;
  ACTION_ARR[OPP_DEL]    = ACTION::OPP_DEL;
  ACTION_ARR[MOV]    = ACTION::MOV    ;
  return 0;
}

MKL_LONG INDEX_MC::copy_INDEX(const INDEX_MC& given_IDX)
{
  initial();
  set_initial_variables();
  for(MKL_LONG i=0; i<4; i++)
    beads[i] = given_IDX.beads[i];
  return 0;
}

INDEX_MC& INDEX_MC::operator=(const INDEX_MC& given_IDX)
{
  copy_INDEX(given_IDX);
  return *this;
}

MKL_LONG INDEX_MC::set_initial_variables()
{
  // Handling the -1 value as initialization is good part for re-usability from other library.
  for(MKL_LONG i=0; i<4; i++)
    beads[i] = -1; // -1 flag means nothing allocated
  return 0;
}


MKL_LONG ACTION::ACT(MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, INDEX_MC& IDX, MATRIX* R_minimum_distance_boost, MKL_LONG const IDENTIFIER_ACTION)
{
  IDX.ACTION_ARR[IDENTIFIER_ACTION](index_t_now, POTs, CONNECT, IDX.beads, R_minimum_distance_boost);
  return 0;
}


MKL_LONG ACTION::UPDATE_INFORMATION(ASSOCIATION& CONNECT, INDEX_MC& IDX, MKL_LONG cnt_arr[], MKL_LONG const IDENTIFIER_ACTION)
{
  for(MKL_LONG i=0; i<IDX.N_BOOST_COUNT[IDENTIFIER_ACTION]; i++)
    {
      CONNECT.update_Z_SUGGESTION_particle(IDX.beads[i]);
      CONNECT.update_dPDF_SUGGESTION_particle(IDX.beads[i]);
      CONNECT.update_dCDF_SUGGESTION_particle(IDX.beads[i]);
    }
  return 0;
}

MKL_LONG ACTION::CANCEL(MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost)
{
  return 0;
}

MKL_LONG ACTION::MOV(MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost)
{
  /*
    It is of importance to nocie that the opp_del_association (i and j) are detaching j chain ends rather than i chain end.
    This is the reason to use 'opp_del' for it's name rather than 'del'.
    In consequence, opp_del_association break the opposite chain ends, then grap this chain end to the subjected particle.
    Some details are described in the opp_del_association_hash file which is the original form in the ASSOCIATION class definition.
    Note that this is related with the Boltzmann distribution generated during stochastic_simulation code is based on the itself particle rather than other particle.
    If this scheme is changed to opposite way, CHAIN_HANDLE class is also affected.
   */
  CONNECT.opp_del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], R_minimum_distance_boost[index_set[CONNECT.flag_other]](index_set[CONNECT.flag_new]));
  
  return 0;
}

MKL_LONG ACTION::OPP_DEL(MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost)
{
  CONNECT.opp_del_association_hash(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_hash_other]);
  return 0;
}

MKL_LONG ACTION::ADD(MKL_LONG index_t_now, POTENTIAL_SET& POTs, ASSOCIATION& CONNECT, MKL_LONG *index_set, MATRIX* R_minimum_distance_boost)
{
  CONNECT.add_association_INFO(POTs, index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new], R_minimum_distance_boost[index_set[CONNECT.flag_itself]](index_set[CONNECT.flag_new]));
  return 0;
}


MKL_LONG ACTION::IDENTIFIER_ACTION_BOOLEAN_BOOST(ASSOCIATION& CONNECT, INDEX_MC& IDX)
{
  // MKL_LONG &index_itself = index_set[2], &index_other = index_set[0], &index_new = index_set[1];
  // check index_other == index_new
  if (CONNECT.CANCEL_ASSOCIATION(CONNECT, IDX.beads))
    {
      // this is the most probable case for canceling
      return IDX.CANCEL;
    }
  // check index_itself == index_new
  else if (CONNECT.ADD_ASSOCIATION(CONNECT, IDX.beads))
    {
      if (CONNECT.CHECK_N_ADD_ASSOCIATION(CONNECT, IDX.beads))
	return IDX.ADD;
      else
        return IDX.CANCEL;
    }
  // check index_itself == index_other
  else if (CONNECT.OPP_DEL_ASSOCIATION(CONNECT, IDX.beads))
    {
      if (CONNECT.CHECK_N_OPP_DEL_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.OPP_DEL;
      else return IDX.CANCEL;
    }
  // with nested case, this is only the case of MOV
  else
    {
      if (CONNECT.CHECK_N_MOV_ASSOCIATION(CONNECT, IDX.beads))
        return IDX.MOV;
      return IDX.CANCEL;
    }
  printf("ERR:TRUTH_MAP::NO SPECIFIED ACTION CASE\n");
  return IDX.CANCEL;
}

MKL_LONG SEARCHING::bisection(MATRIX& given_arr, double p)
{
  MKL_LONG N = given_arr.size;
  MKL_LONG k = N/2;
  // double dk = N/2.;
  MKL_LONG dk = N/2;

  do
    {
      // dk /= 2;
      // the conditional phrase (dk%2 != 0) will return 0 or 1 when the modulo is zero or not, respectively
      // this compensate the loss of approaching because of given domain is in integer rather than real number
      // because of it, the convergence rate is slower than the previous one (log2(Np))
      dk = dk/2 + (dk%2 != 0);
      if (given_arr(k) < p)
        k += dk;
      else if (given_arr(k) > p)
        k -= dk;
      else
        return k;
    } while (dk > 1);
    
  if (given_arr(k) < p)
    k += 1;
  return k;
}

MKL_LONG SEARCHING::backtrace(MATRIX& given_arr, double p)
{
  for(MKL_LONG k= given_arr.size - 1; k >=0; k--)
    {
      if(given_arr(k) < p)
        {
          return k+1;
        }
    }
  return 0;
}

MKL_LONG SEARCHING::backtrace_cell_list(MATRIX& given_arr, MKL_LONG TOKEN, double p, MKL_LONG index_particle, RDIST& R_boost)
{
  // MKL_LONG TOKEN = R_boost.TOKEN[R_boost.cell_index[index_particle]];
  // printf("====\n");

  for(MKL_LONG k= given_arr.size - 1; k >= given_arr.size - TOKEN; k--)
    {
      // printf("dCDF_SUGGESTION[%ld][%ld (min = %ld)] = %4.1e\n", index_particle, k, given_arr.size - TOKEN, given_arr(k));
      if(given_arr(k) < p)
        {
          return k+1;
        }
    }
  return given_arr.size - TOKEN;
}



MKL_LONG CHAIN_HANDLE::hash_initial()
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

  for(MKL_LONG i=0; i<N_chains; i++)
    {
      // note that the initial in association.h already allocate HEAD(i) and TAIL(i) as i%N_particles (all connections are loop)
      // hence, this initialization is only allocate the PARTICLE HASH table for further approaches
      PARTICLE[HEAD(i)][P_TOKEN[HEAD(i)]++] = i;
      PARTICLE[TAIL(i)][P_TOKEN[TAIL(i)]++] = i + N_chains;
    }
  // for(MKL_LONG i=0; i<2*N_chains; i++)
  //   {
  //     /* CE_ATTACHED(i) */
  //     particle_add_CE(CE_ATTACHED_REF(i), i);
  //   }

  // this is for generating random number stream which will be necessary for selecting chain in degenerated state
  const gsl_rng_type *T;
  /* gsl_rng *r; */
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r_degeneracy_check = gsl_rng_alloc(T);
  gsl_rng_set(r_degeneracy_check, random());
  degeneracy_index_array = (MKL_LONG*)mkl_malloc(N_chains*sizeof(MKL_LONG), BIT);
  return 0;
}

// MKL_LONG CHAIN_HANDLE::initialize_without_connection()
// {
//   // zeroing PARTICLE array 
//   for(MKL_LONG i=0; i<N_particles; i++)
//     {
//       P_TOKEN[i] = 0;
//       for(MKL_LONG j=0; j<2*N_chains; j++)
// 	{
// 	  PARTICLE[i][j] = 0;
// 	}
//     }
//   for(MKL_LONG i=0; i<N_chains; i++)
//     {
//       PARTICLE[(MKL_LONG)CE_ATTACHED_REF(i)
//     }
//   return 0;
// }

MKL_LONG CHAIN_HANDLE::allocate_existing_bridges(ASSOCIATION& CONNECT)
{
  /*
    This function is of importance when simulation inherit the existing connectivity information.
    At this moment, however, the function is not fully functional. There are bugs inside (logically).
   */
  // zeroing PARTICLE array 
  for(MKL_LONG i=0; i<N_particles; i++)
    {
      P_TOKEN[i] = 0;
      for(MKL_LONG j=0; j<2*N_chains; j++)
	{
	  PARTICLE[i][j] = 0;
	}
    }
  // zeroing CHAIN information
  for(MKL_LONG i=0; i<2*N_chains; i++)
    CE_ATTACHED_REF(i) = -1;

  MKL_LONG **PAIR_CHECK = (MKL_LONG**)mkl_malloc(N_particles*sizeof(MKL_LONG*), BIT);
  for(MKL_LONG i=0; i<N_particles; i++)
    {
      PAIR_CHECK[i] = (MKL_LONG*)mkl_malloc(N_particles*sizeof(MKL_LONG), BIT);
      for(MKL_LONG j=0; j<N_particles; j++)
	{
	  PAIR_CHECK[i][j] = -1;
	}
    }
  
  MKL_LONG chain_count = 0;
  for(MKL_LONG i=0; i<N_particles; i++)
    {
      for(MKL_LONG k=0; k<CONNECT.TOKEN[i]; k++)
        {
	  MKL_LONG connected_particle = (MKL_LONG)CONNECT.HASH[i](k);
	  if(PAIR_CHECK[i][connected_particle] == -1 && PAIR_CHECK[connected_particle][i] == -1)
	    {
	      MKL_LONG degeneracy = CONNECT.weight[i](k);
	      for(MKL_LONG p=0; p<degeneracy; p++)
		{
		  // CE_ATTACHED_REF(chain_count) = (MKL_LONG)CONNECT.HASH[i](k); // individual chain information
		  CE_ATTACHED_REF(chain_count) = i;
		  CE_ATTACHED_REF(chain_count + N_chains) = connected_particle;
		  // PARTICLE[i][P_TOKEN[i]++] = (MKL_LONG)CONNECT.HASH[i](k); // hash table information
		  PARTICLE[i][P_TOKEN[i]++] = chain_count;
		  PARTICLE[connected_particle][P_TOKEN[connected_particle]++] = chain_count + N_chains;
		  chain_count++;
		  // CHAIN[chain_count].HEAD = i;
		  // CHAIN[chain_count].TAIL = j;
		}
	      PAIR_CHECK[i][connected_particle] = degeneracy;
	      PAIR_CHECK[connected_particle][i] = degeneracy;
	    }
        }
    }
  for(MKL_LONG i=0; i<N_particles; i++)
    mkl_free(PAIR_CHECK[i]);
  mkl_free(PAIR_CHECK);

  if (chain_count != N_chains)
    printf("ERR: allocating_existing_brdiges have unbalanced number between chain_count (%ld) and N_chains (%ld)\n", chain_count, N_chains);
  return 0;
}
