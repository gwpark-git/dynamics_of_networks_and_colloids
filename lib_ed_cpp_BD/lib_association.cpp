
#include "lib_association.h"

long bisection_search(MATRIX& given_arr, double p)
{
  long N = given_arr.size;
  long k = N/2;
  // double dk = N/2.;
  long dk = N/2;

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

long backsearch(MATRIX& given_arr, double p)
{
  for(long k= given_arr.size - 1; k >=0; k--)
    {
      if(given_arr(k) < p)
        {
          return k+1;
        }
    }
  return 0;
}


long ASSOCIATION::read_exist_weight(const char* fn_weight)
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
      long weight_k = 0;
      for(long k=0; k<TOKEN[i]; k++)
        {
          GIVEN_WEIGHT >> weight_k;
          weight[i](k) = weight_k;
        }
    }
  GIVEN_WEIGHT.close();
  return 0;
}
long ASSOCIATION::set_initial_condition()
{
  for(long i=0; i<Np; i++)
    {
      HASH[i].initial(N_max, 1, -1); // this is temporal initialization
      CASE[i].initial(N_max, 1, 0.);
      dPDF[i].initial(N_max, 1, 0.);
      dCDF[i].initial(N_max, 1, 0.);
      Z[i] = 0.;
      weight[i].initial(N_max, 1, 0);
    }
  N_ASSOCIATION = 0;

  for(long i=0; i<Np; i++)
    {
      dPDF[i](0) = 1.0;
      dCDF[i](0) = 1.0;
      Z[i] = 2*Nc;
      HASH[i](0) = i;
      CASE[i](0) = 1.0;
      weight[i](0) = 2*Nc;
      TOKEN[i] = 1;
    }
  return 0;
}

long ASSOCIATION::initial() // it should not be called by outside
{

  CASE = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  dPDF = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  dCDF = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  Z = (double*) mkl_malloc(Np*sizeof(double), BIT);

  weight = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  set_initial_condition();
  return 0;
}

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



ASSOCIATION::ASSOCIATION(long number_of_particles, long number_of_chains_per_particles, long tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS) : CONNECTIVITY(number_of_particles, 2*number_of_chains_per_particles + tolerance_connection)
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


bool TRUTH_MAP::DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new)
{
  if (index_new == index_itself && index_target != index_itself && CONNECT.N_CONNECTED_ENDS(index_itself) < CONNECT.N_max && CONNECT.N_CONNECTED_ENDS(index_target) > CONNECT.N_min)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::NEW_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new)
{
  if (index_new != index_itself && index_target == index_itself && CONNECT.N_CONNECTED_ENDS(index_itself) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_new) < CONNECT.N_max)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new)
{
  // importance update:: index_new != index_itself => index_target
  if(index_new != index_target && index_target != index_itself && CONNECT.N_CONNECTED_ENDS(index_target) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_new) < CONNECT.N_max)
    return TRUE;
  return FALSE;
}


bool TRUTH_MAP::NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new)
{
  if (NEW_ASSOCIATION_BASIC(CONNECT, index_itself, index_target, index_new) && !CONNECT.CHECK_EXIST(index_itself, index_new))
    return TRUE;
  return FALSE;
}


bool TRUTH_MAP::MOV_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, long index_itself, long index_target, long index_new)
{
  if (TRUTH_MAP::MOV_ASSOCIATION_BASIC(CONNECT, index_itself, index_target, index_new) && !CONNECT.CHECK_EXIST(index_itself, index_new))
    return TRUE;
  return FALSE;
}

long ASSOCIATION::N_TOTAL_ASSOCIATION()
{
  long total_association = 0;
  for(long i=0; i<Np; i++)
    {
      for(long j=1; j<TOKEN[i]; j++)
        {
          total_association += (long)weight[i](j);
        }
    }
  return total_association;
}

long ASSOCIATION::N_CONNECTED_ENDS(long given_index)
{
  long cnt = 0;
  for(long j=0; j<TOKEN[given_index]; j++)
    cnt += (long)weight[given_index](j);
  return cnt;
}

long ASSOCIATION::N_TOTAL_CONNECTED_ENDS()
{
  long cnt = 0;
  for(long i=0; i<Np; i++)
    {
      cnt += N_CONNECTED_ENDS(i);
    }
  return cnt;
}

bool ASSOCIATION::CHECK_EXIST_1D(long index_A, long index_B)
{
  if (FIND_HASH_INDEX(index_A, index_B) == TOKEN[index_A])
    return FALSE;
  return TRUE;
}

bool ASSOCIATION::CHECK_EXIST(long index_A, long index_B)
{
  if (FIND_HASH_INDEX(index_A, index_B) == TOKEN[index_A] || FIND_HASH_INDEX(index_B, index_A) == TOKEN[index_B])
    return FALSE;
  return TRUE;
}

long ASSOCIATION::FIND_HASH_INDEX(long index_A, long index_B)
{
  long i=1;
  for(i=1; i<TOKEN[index_A]; i++)
    {
      if ((long)HASH[index_A](i) == index_B)
        return i;
    }
  return TOKEN[index_A];
}


long ASSOCIATION::GET_INDEX_HASH_FROM_ROLL(long index_particle, double rolled_p)
{
  long i=0;
  for(i=0; i<TOKEN[index_particle]; i++)
    {
      if (dCDF[index_particle](i) > rolled_p)
        return i;
    }
  return -1;
}

long ASSOCIATION::GET_HASH_FROM_ROLL(long index_particle, double rolled_p)
{
  return (long)HASH[index_particle](GET_INDEX_HASH_FROM_ROLL(index_particle, rolled_p));
}



long ASSOCIATION::add_association(long index_particle, long index_target)
{
  long hash_index_target = FIND_HASH_INDEX(index_particle, index_target); // if there is no connection, it will return TOKEN(index_particle)

  // weight(index_particle, 0) --; 
  // the chain ends attached to itself decreases with number 1
  // but, one weight should transfer from the itself index, 0, to the other has index
  // for this reason, the weight(index_particle, 0) is decreases with number 2
  // basically, it preserve total number of chain ends on the system
  weight[index_particle](0) -= 2; 
  weight[index_particle](hash_index_target) += 1; 

  if (hash_index_target == TOKEN[index_particle])
    {
      HASH[index_particle](hash_index_target) = index_target;
      TOKEN[index_particle] += 1;
    }

  long hash_index_particle = FIND_HASH_INDEX(index_target, index_particle);
  weight[index_target](hash_index_particle) += 1 ;
  if (hash_index_particle == TOKEN[index_target])
    {
      HASH[index_target](hash_index_particle) = index_particle;
      TOKEN[index_target] += 1;
    }
  // following is test function for the hash table
  // if ((long)HASH[index_particle](TOKEN[index_particle]) != -1) 
  //   {
  //     printf("(ADD) index_target=%ld, hash_index_K=%ld, TOKEN[index_particle]=%ld, HASH[%ld](%ld)=%ld, HASH[%ld](%ld)=%ld\n", index_target, hash_index_target, TOKEN[index_particle], index_particle, TOKEN[index_particle]-1, (long)HASH[index_particle](TOKEN[index_particle]-1), index_particle, TOKEN[index_particle], (long)HASH[index_particle](TOKEN[index_particle]));
  //     for(long q=0; q<HASH[index_particle].size; q++)
  //       {
  //         printf("%ld\t", (long)HASH[index_particle](q));
  //       }
  //     printf("\n");
  //   }

  return 0;
}

long ASSOCIATION::add_association_INFO(POTENTIAL_SET& POTs, long index_particle, long index_target, double distance)
{
  double val_CASE = POTs.w_function(distance, POTs.f_connector(distance, POTs.force_variables), POTs.force_variables);
  add_association(index_particle, index_target);
  CASE[index_particle](FIND_HASH_INDEX(index_particle, index_target)) = val_CASE;
  CASE[index_target](FIND_HASH_INDEX(index_target, index_particle)) = val_CASE;
  return 0;
}
// note that the following code is using hash table.
// however, this case have potential overhead compared with linked-list
// for further implementation, the infrastructure for the interface should be refined.
long ASSOCIATION::del_association_IK(long index_I, long hash_index_K)
{
  weight[index_I](hash_index_K) -= 1;
  if((long)weight[index_I](hash_index_K) == 0)
    {
      for(long j=hash_index_K; j<TOKEN[index_I] - 1; j++)
        {
          // draw for deleted hash position
          HASH[index_I](j) = HASH[index_I](j+1);
          CASE[index_I](j) = CASE[index_I](j+1);
          weight[index_I](j) = weight[index_I](j+1);
        }
      // removing end-tail of the hash table
      HASH[index_I](TOKEN[index_I] - 1) = -1;
      weight[index_I](TOKEN[index_I] - 1) = 0;
      CASE[index_I](TOKEN[index_I] - 1) = 0.;
      TOKEN[index_I] -= 1;
    }
  // following are testing code for hash table
  // for (long j=TOKEN[index_I]; j<HASH[index_I].size; j++)
  //   {
  //     if ((long)HASH[index_I](j) != -1)
  //       {
  //         printf("(DEL) index_I=%ld, hash_index_K=%ld, TOKEN[index_I]=%ld, HASH[index_I](%ld)=%ld, HASH[index_I](%ld)=%ld\n", index_I, hash_index_K, TOKEN[index_I], j-1, (long)HASH[index_I](j-1), j, (long)HASH[index_I](j));
  //         for(long q=0; q<HASH[index_I].size; q++)
  //           {
  //             printf("%ld\t", (long)HASH[index_I](q));
  //           }
  //         printf("\n");
  //       }
  //   }
  return 0; 
}

long ASSOCIATION::del_association_grab_IK(long index_I, long hash_index_K)
{
  del_association_IK(index_I, hash_index_K);
  // CASE(index_I, 0) += 1.0;
  // weight(index_I, 0) ++;
  // The same reason with add_association, the weight on here is increases by number of 2
  weight[index_I](0) += 2;

  return 0;
}

long ASSOCIATION::del_association_hash(long index_particle, long hash_index_target)
{
  long index_target = (long)HASH[index_particle](hash_index_target);
  del_association_grab_IK(index_particle, hash_index_target);
  del_association_IK(index_target, FIND_HASH_INDEX(index_target, index_particle));
  return 0;
}

long ASSOCIATION::del_association(long index_particle, long index_target)
{
  long hash_index_target = FIND_HASH_INDEX(index_particle, index_target);
  del_association_hash(index_particle, hash_index_target);
  return 0;
}


double CONNECTIVITY_update_CASE_particle_hash_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, long index_particle, long hash_index_target, double distance)
{
  CONNECT.CASE[index_particle](hash_index_target) = POTs.w_function(distance, POTs.f_connector(distance, POTs.force_variables), POTs.force_variables);
  return CONNECT.CASE[index_particle](hash_index_target);
}

double CONNECTIVITY_update_CASE_particle_target(ASSOCIATION& CONNECT, POTENTIAL_SET& POTs, long index_particle, long index_target, double distance)
{
  long hash_index_target = CONNECT.FIND_HASH_INDEX(index_particle, index_target);
  return CONNECTIVITY_update_CASE_particle_hash_target(CONNECT, POTs, index_particle, hash_index_target, distance);
}

double CONNECTIVITY_update_Z_particle(ASSOCIATION& CONNECT, long index_particle)
{
  
  // result = dot(x, y)
  CONNECT.Z[index_particle] = cblas_ddot(CONNECT.TOKEN[index_particle], // N
                                         CONNECT.weight[index_particle].data, // x
                                         1,
                                         CONNECT.CASE[index_particle].data, // y
                                         1);
  // note that the previous code is already checked with the below code,
  // CONNECT.Z[index_particle] = 0.;
  // double result = 0.;
  // for(long k=0; k<CONNECT.HASH[index_particle].size; k++)
  //   {
  //     // CONNECT.Z[index_particle] += CONNECT.weight[index_particle](k)*CONNECT.CASE[index_particle](k);
  //     result += CONNECT.weight[index_particle](k)*CONNECT.CASE[index_particle](k);
  //   }
  // if (fabs(CONNECT.Z[index_particle] - result) > 0.0001)
  //   printf("Z err %lf, %lf\n", CONNECT.Z[index_particle], result);
  return CONNECT.Z[index_particle];
}



double CONNECTIVITY_update_dPDF_particle(ASSOCIATION& CONNECT, long index_particle)
{
  // Note that the element-wise multiplication can be done to make diagonalization. BLAS has good aspect for this functionality.
  // DSBMV which gave us y = alpha*A*x + beta*y,
  // which is designed for symmetric banded matrix for A
  // In this case, we can set the third argument, k, set to 0

  // y = alpha*dot(A, x) + beta*y
  cblas_dsbmv(CblasRowMajor,
              CblasUpper,
              CONNECT.TOKEN[index_particle], // number of accounting
              0, // k is related with band. 0 for diagonal (or super-diagonal)
              1./CONNECT.Z[index_particle], // alpha
              CONNECT.weight[index_particle].data, // A (here is array for diagonal components)
              1,  // lda this is leading dimension which inc A. here is set as 1 because the A is only contained diagonal components
              CONNECT.CASE[index_particle].data, // x
              1, // incx
              0., // beta. 0 is removing the previous history because 0*y
              CONNECT.dPDF[index_particle].data, // y
              1); // incy
  // Note that the previous computation has been tested with the below code (the original one)
  // for(long k=0; k<CONNECT.TOKEN[index_particle]; k++)
  //   {
  //     double tmp = (double)CONNECT.weight[index_particle](k)*CONNECT.CASE[index_particle](k)/CONNECT.Z[index_particle];
  //     // this is testing fcn
  //     if (fabs(tmp - CONNECT.dPDF[index_particle](k)) > 0.00001 )
  //       printf("UPDATE: k=%ld, tmp=%lf, dPDF[index_particle](k)=%lf\n", k, tmp, CONNECT.dPDF[index_particle](k));
  //   }
  return 0;
}

double CONNECTIVITY_update_dCDF_particle(ASSOCIATION& CONNECT, long index_particle)
{
  CONNECT.dCDF[index_particle](0) = CONNECT.dPDF[index_particle](0);
  for(long k=1; k<CONNECT.TOKEN[index_particle]; k++)
    {
      CONNECT.dCDF[index_particle](k) = CONNECT.dCDF[index_particle](k-1) + CONNECT.dPDF[index_particle](k);
    }
  return CONNECT.dCDF[index_particle](CONNECT.TOKEN[index_particle]);
}


