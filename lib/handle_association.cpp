

#include "handle_association.h"

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
      CONNECT.update_Z_particle(IDX.beads[i]);
      CONNECT.update_dPDF_particle(IDX.beads[i]);
      CONNECT.update_dCDF_particle(IDX.beads[i]);
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
      // printf("dCDF[%ld][%ld (min = %ld)] = %4.1e\n", index_particle, k, given_arr.size - TOKEN, given_arr(k));
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
