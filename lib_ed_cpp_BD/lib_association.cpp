
#include "lib_association.h"

MKL_LONG bisection_search_log2Np(MATRIX& given_arr, double p, MKL_LONG log2_Np)
{
  MKL_LONG N = given_arr.size;
  MKL_LONG k = N/2;
  double dk = N/2.;
  for(MKL_LONG i=0; i<=log2_Np; i++)
    {
      dk /= 2.;
      if (given_arr(k) < p)
        k += (MKL_LONG)dk;
      else if (given_arr(k) > p)
        k -= (MKL_LONG)dk;
      else
        return k;
    }
  if (given_arr(k) < p)
    k+= 1;
  return k;
}

MKL_LONG bisection_search(MATRIX& given_arr, double p)
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

MKL_LONG backsearch(MATRIX& given_arr, double p)
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


MKL_LONG ASSOCIATION::read_exist_weight(const char* fn_weight)
{
  ifstream GIVEN_WEIGHT;
  GIVEN_WEIGHT.open(fn_weight);
  long cnt = 0;
  string line;
  while(getline(GIVEN_WEIGHT, line))
    {
      cnt++;
    }
  GIVEN_WEIGHT.clear();
  GIVEN_WEIGHT.seekg(0);
  for(long i=0; i<cnt-Np; i++)
    {
      getline(GIVEN_WEIGHT, line);
    }
  for(long i=0; i<Np; i++)
    {
      MKL_LONG weight_k = 0;
      for(long k=0; k<TOKEN(i); k++)
        {
          GIVEN_WEIGHT >> weight_k;
          weight(k) = weight_k;
        }
    }
  GIVEN_WEIGHT.close();
  return 0;
}

MKL_LONG ASSOCIATION::initial()
{
  CASE.initial(Np, N_max, 0.);
  dPDF.initial(Np, N_max, 0.);
  dCDF.initial(Np, N_max, 0.);
  Z.initial(Np, 1, 0.);
  HASH.initial(Np, N_max, -1);
  weight.initial(Np, N_max, 0);
  dist_map.initial(Np, Np, 0);
  N_ASSOCIATION = 0;

  for(MKL_LONG i=0; i<Np; i++)
    {
      dPDF(i, 0) = 1.0;
      dCDF(i, 0) = 1.0;
      Z(i, 0) = 2*Nc;
      HASH(i, 0) = i;
      CASE(i, 0) = 1.0;
      weight(i, 0) = 2*Nc;
    }
  return 0;
}

// ASSOCIATION::ASSOCIATION(TRAJECTORY& TRAJ, COND& given_condition) : CONNECTIVITY(TRAJ.Np, 2*atol(given_condition("N_chains_per_particle").c_str()) + atol(given_condition("tolerance_allowing_connections").c_str()))
ASSOCIATION::ASSOCIATION(TRAJECTORY& TRAJ, COND& given_condition) : CONNECTIVITY(given_condition)
{
  Nc = atol(given_condition("N_chains_per_particle").c_str());
  Tec = atol(given_condition("tolerance_allowing_connections").c_str());
  N_min = 2*Nc - Tec;
  N_max = 2*Nc + Tec;
  initial();

  if (given_condition("CONTINUATION_CONNECTION")=="TRUE")
    {
      read_exist_weight(given_condition("CONTINUATION_WEIGHT_FN").c_str());
    }
  if (given_condition("allowing_multiple_connections") == "TRUE")
    {
      NEW_ASSOCIATION = TRUTH_MAP::NEW_ASSOCIATION_BASIC;
      MOV_ASSOCIATION = TRUTH_MAP::MOV_ASSOCIATION_BASIC;
    }
  else
    {
      NEW_ASSOCIATION = TRUTH_MAP::NEW_ASSOCIATION_SINGLE;
      MOV_ASSOCIATION = TRUTH_MAP::MOV_ASSOCIATION_SINGLE;
    }
  DEL_ASSOCIATION = TRUTH_MAP::DEL_ASSOCIATION_BASIC;
}



ASSOCIATION::ASSOCIATION(MKL_LONG number_of_particles, MKL_LONG number_of_chains_per_particles, MKL_LONG tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS) : CONNECTIVITY(number_of_particles, 2*number_of_chains_per_particles + tolerance_connection)
{
  Nc = number_of_chains_per_particles;
  Tec = tolerance_connection;
  N_min = 2*Nc - Tec;
  N_max = 2*Nc + Tec;
  initial();
  if(ALLOWING_MULTIPLE_CONNECTIONS)
    {
      NEW_ASSOCIATION = TRUTH_MAP::NEW_ASSOCIATION_BASIC;
      MOV_ASSOCIATION = TRUTH_MAP::MOV_ASSOCIATION_BASIC;
    }
  else
    {
      NEW_ASSOCIATION = TRUTH_MAP::NEW_ASSOCIATION_SINGLE;
      MOV_ASSOCIATION = TRUTH_MAP::MOV_ASSOCIATION_SINGLE;
    }
  DEL_ASSOCIATION = TRUTH_MAP::DEL_ASSOCIATION_BASIC;
}

// MKL_LONG ASSOCIATION::TEST_VALID()
// {
//   for(MKL_LONG i=0; i<TOKEN(index_A); i++)
//     {

//     }
// }

bool TRUTH_MAP::DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  if (index_new == index_itself && index_target != index_itself && CONNECT.N_CONNECTED_ENDS(index_itself) < CONNECT.N_max && CONNECT.N_CONNECTED_ENDS(index_target) > CONNECT.N_min)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::NEW_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  if (index_new != index_itself && index_target == index_itself && CONNECT.N_CONNECTED_ENDS(index_itself) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_new) < CONNECT.N_max)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  if (NEW_ASSOCIATION_BASIC(CONNECT, index_itself, index_target, index_new) && !CONNECT.CHECK_EXIST(index_itself, index_new))
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  // importance update:: index_new != index_itself => index_target
  if(index_new != index_target && index_target != index_itself && CONNECT.N_CONNECTED_ENDS(index_target) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_new) < CONNECT.N_max)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::MOV_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  if (TRUTH_MAP::MOV_ASSOCIATION_BASIC(CONNECT, index_itself, index_target, index_new) && !CONNECT.CHECK_EXIST(index_itself, index_new))
    return TRUE;
  return FALSE;
}

MKL_LONG ASSOCIATION::N_TOTAL_ASSOCIATION()
{
  MKL_LONG total_association = 0;
// #pragma omp parallel for default(none) shared(total_association)
  for(MKL_LONG i=0; i<Np; i++)
    {
      for(MKL_LONG j=1; j<TOKEN(i); j++)
        {
          total_association += weight(i,j);
        }
      // total_association += TOKEN(i) - 1;
    }
  return total_association;
}

MKL_LONG ASSOCIATION::N_CONNECTED_ENDS(MKL_LONG given_index)
{
  MKL_LONG cnt = 0;
  for(MKL_LONG j=0; j<TOKEN(given_index); j++)
    cnt += weight(given_index, j);
  return cnt;
}

MKL_LONG ASSOCIATION::N_TOTAL_CONNECTED_ENDS()
{
  MKL_LONG cnt = 0;
  for(MKL_LONG i=0; i<Np; i++)
    {
      cnt += N_CONNECTED_ENDS(i);
    }
  return cnt;
}

bool ASSOCIATION::CHECK_EXIST_1D(MKL_LONG index_A, MKL_LONG index_B)
{
  if (FIND_HASH_INDEX(index_A, index_B) == TOKEN(index_A))
    return FALSE;
  return TRUE;
}

bool ASSOCIATION::CHECK_EXIST(MKL_LONG index_A, MKL_LONG index_B)
{
  if (FIND_HASH_INDEX(index_A, index_B) == TOKEN(index_A) || FIND_HASH_INDEX(index_B, index_A) == TOKEN(index_B))
    return FALSE;
  return TRUE;
}

MKL_LONG ASSOCIATION::FIND_HASH_INDEX(MKL_LONG index_A, MKL_LONG index_B)
{
  MKL_LONG i=1;
  for(i=1; i<TOKEN(index_A); i++)
    {
      if (HASH(index_A, i) == index_B)
        return i;
    }
  // printf("CANNOT FOUND THE HASH\n");
  return TOKEN(index_A);
}

// bool ASSOCIATION::CHECK_NC_ADD_BTA(MKL_LONG index_A, MKL_LONG index_B)
// {
//   if (weight(index_A, 0) + 1< N_max)
//     return TRUE;
//   return FALSE;
// }

// // bool ASSOCIATION::CHECK_NC_ADD_BTA(MKL_LONG index_A, MKL_LONG index_B)
// // {
// //   if( CHECK_NC_ADD_BTA_1D(index_A, index_B) && CHECK_NC_ADD_BTA_1D(index_B, index_A))
// //     return TRUE;
// //   return FALSE;
// // }

// bool ASSOCIATION::CHECK_NC_DEL_hBFA(MKL_LONG index_A, MKL_LONG index_hash_B)
// {
//   if (weight(index_A, 0) - 1 > N_min)
//     return TRUE;
//   return FALSE;
// }

// bool ASSOCIATION::CHECK_NC_DEL_BFA(MKL_LONG index_A, MKL_LONG index_B)
// {
//   return CHECK_NC_DEL_hBFA(index_A, FIND_HASH_INDEX(index_A, index_B));
// }

// bool ASSOCIATION::CHECK_NC_DEL_BFA(MKL_LONG index_A, MKL_LONG index_B)
// {
//   if(CHECK_NC_DEL_BFA_1D(index_A, index_B) && CHECK_NC_DEL_BFA_1D(index_B, index_A))
//     return TRUE;
//   return FALSE;
// }


MKL_LONG ASSOCIATION::GET_INDEX_HASH_FROM_ROLL(MKL_LONG index_particle, double rolled_p)
{
  MKL_LONG i=0;
  // bool IDENTIFIER = TRUE;
  for(i=0; i<TOKEN(index_particle); i++)
    {
      if (dCDF(index_particle, i) > rolled_p)
        return i;
    }
  // if(HASH(index_particle, i) == -1)
  //   printf("HASH index overflow: (i_P, i_H, i_T, dCDF, rolled_P) = (%ld, %ld, %ld, %6.3f, %6.3f)\n", index_particle, i, HASH(index_particle, i), dCDF(index_particle, i), rolled_p);
  return -1;
}

MKL_LONG ASSOCIATION::GET_HASH_FROM_ROLL(MKL_LONG index_particle, double rolled_p)
{
  return HASH(index_particle, GET_INDEX_HASH_FROM_ROLL(index_particle, rolled_p));
}


// double& ASSOCIATION::N_CHAIN_ENDS(MKL_LONG index_particle)
// {
//   return CASE(index_particle, 0);
// }

MKL_LONG ASSOCIATION::add_association(MKL_LONG index_particle, MKL_LONG index_target)
{
  // HASH(index_particle, TOKEN(index_particle)) = index_target;
  MKL_LONG hash_index_target = FIND_HASH_INDEX(index_particle, index_target); // if there is no connection, it will return TOKEN(index_particle)

  // weight(index_particle, 0) --; 
  // the chain ends attached to itself decreases with number 1
  // but, one weight should transfer from the itself index, 0, to the other has index
  // for this reason, the weight(index_particle, 0) is decreases with number 2
  // basically, it preserver total number of chain ends on the system
  weight(index_particle, 0) -= 2; 
  weight(index_particle, hash_index_target) ++; 

  if (hash_index_target == TOKEN(index_particle))
    {
      HASH(index_particle, hash_index_target) = index_target;
      TOKEN(index_particle) += 1;
    }

  MKL_LONG hash_index_particle = FIND_HASH_INDEX(index_target, index_particle);
  weight(index_target, hash_index_particle) ++ ;
  if (hash_index_particle == TOKEN(index_target))
    {
      HASH(index_target, hash_index_particle) = index_particle;
      TOKEN(index_target) += 1;
    }
  return 0;
}

MKL_LONG ASSOCIATION::add_association_INFO(POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG index_target, double distance)
{
  double val_CASE = POTs.w_function(distance, POTs.f_connector(distance, POTs.force_variables), POTs.force_variables);
  add_association(index_particle, index_target);
  CASE(index_particle, FIND_HASH_INDEX(index_particle, index_target)) = val_CASE;
  CASE(index_target, FIND_HASH_INDEX(index_target, index_particle)) = val_CASE;
  // MKL_LONG hash_index_reverse_target = add_association(index_target, index_particle);
  // CASE(index_target, hash_index_reverse_target) = val_CASE;
  return 0;
  // CASE(index_particle, TOKEN(index_particle)) = val_CASE;
  // CASE(index_target, TOKEN(index_target)) = val_CASE;
  // return add_association(index_particle, index_target);
}
// note that the following code is using hash table.
// however, this case have potential overhead compared with linked-list
// for further implementation, the infrastructure for the interface should be refined.
MKL_LONG ASSOCIATION::del_association_IK(MKL_LONG index_I, MKL_LONG hash_index_K)
{
  // if(weight(index_I, hash_index_K) > 1)
  //   {
  //     weight(index_I, hash_index_K) --;
  //   }
  // else
  //   {
  if(--weight(index_I, hash_index_K) == 0)
    {
      for(MKL_LONG j=hash_index_K; j<TOKEN(index_I) - 1; j++)
        {
          // draw for deleted hash position
          HASH(index_I, j) = HASH(index_I, j+1);
          CASE(index_I, j) = CASE(index_I, j+1);
          weight(index_I, j) = weight(index_I, j+1);
        }
      // removing end-tail of the hash table
      HASH(index_I, TOKEN(index_I) - 1) = -1;
      weight(index_I, TOKEN(index_I) - 1) = 0;
      CASE(index_I, TOKEN(index_I) - 1) = 0.;
      TOKEN(index_I) -= 1;
    }
  
  // N_CHAIN_ENDS(index_I) ++;
  // CASE(index_I, 0) += 1.0;
  return 0; 
}

MKL_LONG ASSOCIATION::del_association_grab_IK(MKL_LONG index_I, MKL_LONG hash_index_K)
{
  del_association_IK(index_I, hash_index_K);
  // CASE(index_I, 0) += 1.0;
  // weight(index_I, 0) ++;
  // The same reason with add_association, the weight on here is increases by number of 2
  weight(index_I, 0) += 2;

  return 0;
}

MKL_LONG ASSOCIATION::del_association_hash(MKL_LONG index_particle, MKL_LONG hash_index_target)
{
  MKL_LONG index_target = HASH(index_particle, hash_index_target);
  del_association_grab_IK(index_particle, hash_index_target);
  del_association_IK(index_target, FIND_HASH_INDEX(index_target, index_particle));
  return 0;
}

MKL_LONG ASSOCIATION::del_association(MKL_LONG index_particle, MKL_LONG index_target)
{
  MKL_LONG hash_index_target = FIND_HASH_INDEX(index_particle, index_target);
  del_association_hash(index_particle, hash_index_target);
  // del_association_grab_IK(index_particle, hash_index_target);
  // del_association_IK(index_target, FIND_HASH_INDEX(index_target, index_particle));
  // for(MKL_LONG i=1; i<TOKEN(index_particle); i++)
  //   {
  //     if (HASH(index_particle, i) == index_target)
  //       {
  //         del_association_grab_IK(index_particle, i);
  //         del_association_IK(HASH(index_particle, i), FIND_HASH_INDEX(index_particle)); // reversible deleting
  //         break; 
  //       }
  //   }
  return 0;
}


double CONNECTIVITY_update_CASE_particle_hash_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG hash_index_target, double distance)
{
  CONNECT.CASE(index_particle, hash_index_target) = POTs.w_function(distance, POTs.f_connector(distance, POTs.force_variables), POTs.force_variables);
  return CONNECT.CASE(index_particle, hash_index_target);
}

double CONNECTIVITY_update_CASE_particle_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, MKL_LONG index_target, double distance)
{
  // double distance = GEOMETRY::return_minimum_distance(TRAJ, index_t, index_particle, CONNECT.HASH(index_particle, k));
  MKL_LONG hash_index_target = CONNECT.FIND_HASH_INDEX(index_particle, index_target);
  return CONNECTIVITY_update_CASE_particle_hash_target(CONNECT, POTs, index_particle, hash_index_target, distance);
}


// MKL_LONG CONNECTIVITY_update_CASE_particle(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, MKL_LONG index_particle, double distance)
// {
//   // MKL_LONG Ndc = CONNECT.TOKEN(index_particle) - 1;
//   // CONNECT.CASE(index_particle, 0) = 2*CONNECT.Nc - Ndc;
//   for(MKL_LONG k=0; k<CONNECT.TOKEN(index_particle); k++)
//     {
//       CONNECTIVITY_update_CASE_particle_target(CONNECT, POTs, index_particle, k, distance);
//     }
//   return 0;
// }

double CONNECTIVITY_update_Z_particle(ASSOCIATION& CONNECT, MKL_LONG index_particle)
{
  CONNECT.Z(index_particle) = 0.;
  for(MKL_LONG k=0; k<CONNECT.TOKEN(index_particle); k++)
    {
      CONNECT.Z(index_particle) += (double)CONNECT.weight(index_particle, k)*CONNECT.CASE(index_particle, k);
    }
  return CONNECT.Z(index_particle);
}



double CONNECTIVITY_update_dPDF_particle(ASSOCIATION& CONNECT, MKL_LONG index_particle)
{
  // double sum = 0.;
  // CONNECT.dPDF(index_particle, 0) = 0.;

// #pragma omp parallel for default(none) shared(CONNECT, index_particle)
  for(MKL_LONG k=0; k<CONNECT.TOKEN(index_particle); k++)
    {
      CONNECT.dPDF(index_particle, k) = (double)CONNECT.weight(index_particle, k)*CONNECT.CASE(index_particle, k)/CONNECT.Z(index_particle);
      // sum += CONNECT.P(index_particle, k);
    }
  // return sum;
  return 0;
}

double CONNECTIVITY_update_dCDF_particle(ASSOCIATION& CONNECT, MKL_LONG index_particle)
{
  CONNECT.dCDF(index_particle, 0) = CONNECT.dPDF(index_particle, 0);
  for(MKL_LONG k=1; k<CONNECT.TOKEN(index_particle); k++)
    {
      CONNECT.dCDF(index_particle, k) = CONNECT.dCDF(index_particle, k-1) + CONNECT.dPDF(index_particle, k);
    }
  return CONNECT.dCDF(index_particle, CONNECT.TOKEN(index_particle));
}


