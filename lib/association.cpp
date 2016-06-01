
#include "association.h"

MKL_LONG index_set_val(size_t* index, MKL_LONG given_val, MKL_LONG size_of_arr)
{
  for(MKL_LONG i=0; i<size_of_arr; i++)
    {
      index[i] = given_val;
    }
  return 0;
}

MKL_LONG ASSOCIATION::read_exist_weight(const char* fn_weight)
{
  ifstream GIVEN_WEIGHT;
  GIVEN_WEIGHT.open(fn_weight);
  if (!GIVEN_WEIGHT)
    {
      printf("ERR to read %s", fn_weight);
      return 0;
    }
  MKL_LONG cnt = 0;
  string line;
  while(getline(GIVEN_WEIGHT, line))
    {
      cnt++;
    }
  GIVEN_WEIGHT.clear();
  GIVEN_WEIGHT.seekg(0);
  for(MKL_LONG i=0; i<cnt-Np; i++)
    {
      getline(GIVEN_WEIGHT, line);
    }

  for(MKL_LONG i=0; i<Np; i++)
    {
      MKL_LONG weight_k = 0;
      for(MKL_LONG k=0; k<TOKEN[i]; k++)
        {
          GIVEN_WEIGHT >> weight_k;
          weight[i](k) = weight_k;
        }
      // for(MKL_LONG k=0; k<weight[i].size; k++)
      //   {
      //     GIVEN_WEIGHT >> weight_k;
      //     weight[i](k) = weight_k; // note that weight without association is zero.
      //   }
      
    }
  GIVEN_WEIGHT.close();
  return 0;
}

MKL_LONG ASSOCIATION::set_initial_condition()
{
  for(MKL_LONG i=0; i<Np; i++)
    {
      HASH[i].initial(N_max, 1, -1); // this is temporal initialization
      weight[i].initial(N_max, 1, 0);

      // related with suggestion probability
      // note that it is related with N_max ~ 10*Nc + Tec << Np
      CASE_SUGGESTION[i].initial(N_max, 1, 0.);
      dPDF_SUGGESTION[i].initial(N_max, 1, 0.);
      dCDF_SUGGESTION[i].initial(N_max, 1, 0.);
      Z_SUGGESTION[i] = 0.;

      // related with association probability
      // note that the size of association is related with number of particles
      TOKEN_ASSOCIATION[i] = 0;
      dCDF_ASSOCIATION[i].initial(Np, 1, 0.);
      INDEX_ASSOCIATION[i].initial(Np, 1, -1);
      // index_set_val(INDEX_ASSOCIATION[i], -1, Np);
    }
  N_ASSOCIATION = 0;

  for(MKL_LONG i=0; i<Np; i++)
    {
      HASH[i](0) = i;
      weight[i](0) = 2*Nc;
      TOKEN[i] = 1;

      // related with suggestion probability
      dPDF_SUGGESTION[i](0) = 1.0;
      dCDF_SUGGESTION[i](0) = 1.0;
      Z_SUGGESTION[i] = 2*Nc;
      CASE_SUGGESTION[i](0) = 1.0;
    }
  return 0;
}

MKL_LONG ASSOCIATION::dynamic_alloc()
{
  // intrinsic association member functions
  // weight = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  weight = new MATRIX [Np];
  // related with suggestion probability
  // CASE_SUGGESTION = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  // dPDF_SUGGESTION = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  // dCDF_SUGGESTION = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  CASE_SUGGESTION = new MATRIX [Np];
  dPDF_SUGGESTION = new MATRIX [Np];
  dCDF_SUGGESTION = new MATRIX [Np];
  // Z_SUGGESTION = (double*) mkl_malloc(Np*sizeof(double), BIT);
  Z_SUGGESTION = new double [Np];
  
  // related with association probability
  // dCDF_ASSOCIATION = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  // INDEX_ASSOCIATION = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  dCDF_ASSOCIATION = new MATRIX [Np];
  INDEX_ASSOCIATION = new MATRIX [Np];
  // INDEX_ASSOCIATION = (size_t**) mkl_malloc(Np*sizeof(size_t*), BIT);
  // for(MKL_LONG i=0; i<Np; i++)
  //   INDEX_ASSOCIATION[i] = (size_t*) mkl_malloc(Np*sizeof(size_t), BIT);
  
  // TOKEN_ASSOCIATION = (MKL_LONG*) mkl_malloc(Np*sizeof(MKL_LONG), BIT);
  TOKEN_ASSOCIATION = new MKL_LONG [Np];
    
  if(MULTIPLE_CONNECTIONS)
    {
      CHECK_N_ADD_ASSOCIATION = TRUTH_MAP::MULTIPLE::CHECK_N_ADD_BOOST;
      CHECK_N_MOV_ASSOCIATION = TRUTH_MAP::MULTIPLE::CHECK_N_MOV_BOOST;
      CHECK_N_OPP_DEL_ASSOCIATION = TRUTH_MAP::MULTIPLE::CHECK_N_OPP_DEL_BOOST;
    }
  else
    {
      CHECK_N_ADD_ASSOCIATION = TRUTH_MAP::SINGLE::CHECK_N_ADD_BOOST;
      CHECK_N_MOV_ASSOCIATION = TRUTH_MAP::SINGLE::CHECK_N_MOV_BOOST;
      CHECK_N_OPP_DEL_ASSOCIATION = TRUTH_MAP::SINGLE::CHECK_N_OPP_DEL_BOOST;
    }
  ADD_ASSOCIATION = TRUTH_MAP::ADD_ASSOCIATION_BOOST;
  OPP_DEL_ASSOCIATION = TRUTH_MAP::OPP_DEL_ASSOCIATION_BOOST;
  CANCEL_ASSOCIATION = TRUTH_MAP::CANCEL_ASSOCIATION_BOOST;

  return 0;
}

MKL_LONG ASSOCIATION::initial() // it should not be called by outside
{
  dynamic_alloc();
  // set_initial_condition reset all the variables inclduing variable belong to CONNECTIVITY class
  // this violate the inheritance of association information from the previous computation
  // hence, even if it is quite duplicate, the set_initial_condition() is replaced by the real sequence

  set_initial_condition();
  

  return 0;
}

MKL_LONG ASSOCIATION::initial_inheritance() // it should not be called by outside
{

  dynamic_alloc();
  //   // set_initial_condition reset all the variables inclduing variable belong to CONNECTIVITY class
  //   // this violate the inheritance of association information from the previous computation
  //   // hence, even if it is quite duplicate, the set_initial_condition() is replaced by the real sequence
  // set_initial_condition();
  for(MKL_LONG i=0; i<Np; i++)
    {
      weight[i].initial(N_max, 1, 0);

      // related with suggestion probability
      CASE_SUGGESTION[i].initial(N_max, 1, 0.);
      dPDF_SUGGESTION[i].initial(N_max, 1, 0.);
      dCDF_SUGGESTION[i].initial(N_max, 1, 0.);
      Z_SUGGESTION[i] = 0.;

      // related with association probability
      TOKEN_ASSOCIATION[i] = 0;
      dCDF_ASSOCIATION[i].initial(Np, 1, 0.);
      INDEX_ASSOCIATION[i].initial(Np, 1, -1);
      // index_set_val(INDEX_ASSOCIATION[i], -1, Np);
    }
  
  for(MKL_LONG i=0; i<Np; i++)
    {
      dPDF_SUGGESTION[i](0) = 1.0;
      dCDF_SUGGESTION[i](0) = 1.0;
      CASE_SUGGESTION[i](0) = 1.0;
    }

  return 0;
}

ASSOCIATION::ASSOCIATION(COND& given_condition) : CONNECTIVITY(given_condition)
{
  Nc = atol(given_condition("N_chains_per_particle").c_str());
  Tec = atol(given_condition("tolerance_allowing_connections").c_str());
  N_min = 2*Nc - Tec;
  N_max = 2*Nc + Tec;
  if(given_condition("allowing_multiple_connections") == "TRUE")
    MULTIPLE_CONNECTIONS = TRUE;

  if (given_condition("CONTINUATION_CONNECTION")=="TRUE")
    {
      initial_inheritance();      
      read_exist_weight(given_condition("CONTINUATION_WEIGHT_FN").c_str());
      N_ASSOCIATION = N_TOTAL_ASSOCIATION();
    }
  else
    {
      initial();
    }
}

ASSOCIATION::ASSOCIATION(MKL_LONG number_of_particles, MKL_LONG number_of_chains_per_particles, MKL_LONG tolerance_connection, bool ALLOWING_MULTIPLE_CONNECTIONS) : CONNECTIVITY(number_of_particles, 2*number_of_chains_per_particles + tolerance_connection)
{
  Nc = number_of_chains_per_particles;
  Tec = tolerance_connection;
  N_min = 2*Nc - Tec;
  N_max = 2*Nc + Tec;
  MULTIPLE_CONNECTIONS = ALLOWING_MULTIPLE_CONNECTIONS;
  initial();
}



bool TRUTH_MAP::MULTIPLE::CHECK_N_ADD_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if (CONNECT.N_CONNECTED_ENDS(index_set[CONNECT.flag_itself]) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_set[CONNECT.flag_new]) < CONNECT.N_max )
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::MULTIPLE::CHECK_N_OPP_DEL_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if (CONNECT.N_CONNECTED_ENDS(index_set[CONNECT.flag_itself]) < CONNECT.N_max && CONNECT.N_CONNECTED_ENDS(index_set[CONNECT.flag_other]) > CONNECT.N_min)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::MULTIPLE::CHECK_N_MOV_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if (CONNECT.N_CONNECTED_ENDS(index_set[CONNECT.flag_other]) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_set[CONNECT.flag_new]) < CONNECT.N_max)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::SINGLE::CHECK_N_ADD_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if(!CONNECT.CHECK_EXIST(index_set[CONNECT.flag_itself], index_set[CONNECT.flag_new]))
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::SINGLE::CHECK_N_MOV_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  return TRUTH_MAP::SINGLE::CHECK_N_ADD_BOOST(CONNECT, index_set); // it has same condition
}

bool TRUTH_MAP::SINGLE::CHECK_N_OPP_DEL_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  return TRUE; // this is because the single connection case
}

bool TRUTH_MAP::CANCEL_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if (index_set[CONNECT.flag_other] == index_set[CONNECT.flag_new])
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::ADD_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if (index_set[CONNECT.flag_itself] == index_set[CONNECT.flag_other])
    return TRUE;
  return FALSE;
}


bool TRUTH_MAP::OPP_DEL_ASSOCIATION_BOOST(ASSOCIATION& CONNECT, MKL_LONG index_set[])
{
  if (index_set[CONNECT.flag_itself] == index_set[CONNECT.flag_new])
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::OPP_DEL_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
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

bool TRUTH_MAP::MOV_ASSOCIATION_BASIC(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  // importance update:: index_new != index_itself => index_target
  if(index_new != index_target && index_target != index_itself && CONNECT.N_CONNECTED_ENDS(index_target) > CONNECT.N_min && CONNECT.N_CONNECTED_ENDS(index_new) < CONNECT.N_max)
    return TRUE;
  return FALSE;
}

bool TRUTH_MAP::NEW_ASSOCIATION_SINGLE(ASSOCIATION& CONNECT, MKL_LONG index_itself, MKL_LONG index_target, MKL_LONG index_new)
{
  if (NEW_ASSOCIATION_BASIC(CONNECT, index_itself, index_target, index_new) && !CONNECT.CHECK_EXIST(index_itself, index_new))
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
  for(MKL_LONG i=0; i<Np; i++)
    {
      for(MKL_LONG j=1; j<TOKEN[i]; j++)
        {
          total_association += (MKL_LONG)weight[i](j);
        }
    }
  return total_association/2.;
}

MKL_LONG ASSOCIATION::N_CONNECTED_ENDS(MKL_LONG given_index)
{
  MKL_LONG cnt = 0;
  for(MKL_LONG j=0; j<TOKEN[given_index]; j++)
    cnt += (MKL_LONG)weight[given_index](j);
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
  if (FIND_HASH_INDEX(index_A, index_B) == TOKEN[index_A])
    return FALSE;
  return TRUE;
}

bool ASSOCIATION::CHECK_EXIST(MKL_LONG index_A, MKL_LONG index_B)
{
  if (FIND_HASH_INDEX(index_A, index_B) == TOKEN[index_A] || FIND_HASH_INDEX(index_B, index_A) == TOKEN[index_B])
    return FALSE;
  return TRUE;
}

MKL_LONG ASSOCIATION::FIND_HASH_INDEX(MKL_LONG index_A, MKL_LONG index_B)
{
  MKL_LONG i=1;
  for(i=1; i<TOKEN[index_A]; i++)
    {
      if ((MKL_LONG)HASH[index_A](i) == index_B)
        return i;
    }
  return TOKEN[index_A];
}


MKL_LONG ASSOCIATION::GET_INDEX_HASH_FROM_ROLL(const MKL_LONG index_particle, double rolled_p)
{
  MKL_LONG i=0;
  for(i=0; i<TOKEN[index_particle]; i++)
    {
      if (dCDF_SUGGESTION[index_particle](i) > rolled_p)
        return i;
    }
  return -1;
}

MKL_LONG ASSOCIATION::GET_HASH_FROM_ROLL(const MKL_LONG index_particle, double rolled_p)
{
  return (MKL_LONG)HASH[index_particle](GET_INDEX_HASH_FROM_ROLL(index_particle, rolled_p));
}



MKL_LONG ASSOCIATION::add_association(const MKL_LONG index_particle, MKL_LONG index_target)
{
  MKL_LONG hash_index_target = FIND_HASH_INDEX(index_particle, index_target); // if there is no connection, it will return TOKEN(index_particle)

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

  MKL_LONG hash_index_particle = FIND_HASH_INDEX(index_target, index_particle);
  weight[index_target](hash_index_particle) += 1 ;
  if (hash_index_particle == TOKEN[index_target])
    {
      HASH[index_target](hash_index_particle) = index_particle;
      TOKEN[index_target] += 1;
    }

  return 0;
}

MKL_LONG ASSOCIATION::add_association_INFO(POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG index_target, double distance)
{
  double val_CASE_SUGGESTION = POTs.w_function(distance, POTs.f_connector(distance, POTs.force_variables), POTs.force_variables);
  add_association(index_particle, index_target);
  CASE_SUGGESTION[index_particle](FIND_HASH_INDEX(index_particle, index_target)) = val_CASE_SUGGESTION;
  CASE_SUGGESTION[index_target](FIND_HASH_INDEX(index_target, index_particle)) = val_CASE_SUGGESTION;
  return 0;
}

// note that the following code is using hash table.
// however, this case have potential overhead compared with linked-list
// for further implementation, the infrastructure for the interface should be refined.
MKL_LONG ASSOCIATION::opp_del_association_IK(MKL_LONG index_I, MKL_LONG hash_index_K)
{
  weight[index_I](hash_index_K) -= 1;
  if((MKL_LONG)weight[index_I](hash_index_K) == 0)
    {
      for(MKL_LONG j=hash_index_K; j<TOKEN[index_I] - 1; j++)
        {
          // draw for deleted hash position
          HASH[index_I](j) = HASH[index_I](j+1);
          CASE_SUGGESTION[index_I](j) = CASE_SUGGESTION[index_I](j+1);
          weight[index_I](j) = weight[index_I](j+1);
        }
      // removing end-tail of the hash table
      HASH[index_I](TOKEN[index_I] - 1) = -1;
      weight[index_I](TOKEN[index_I] - 1) = 0;
      CASE_SUGGESTION[index_I](TOKEN[index_I] - 1) = 0.;
      TOKEN[index_I] -= 1;
    }
  return 0; 
}

MKL_LONG ASSOCIATION::opp_del_association_grab_IK(MKL_LONG index_I, MKL_LONG hash_index_K)
{
  opp_del_association_IK(index_I, hash_index_K);
  // The same reason with add_association, the weight on here is increases by number of 2
  weight[index_I](0) += 2;

  return 0;
}

MKL_LONG ASSOCIATION::opp_del_association_hash(const MKL_LONG index_particle, MKL_LONG hash_index_target)
{
  /*
    This is the basic call for deleting the existing connection.
    The subjected chain ends is opposite chain ends that selected at this moment, because for the association distribution subjected by itself particle, which means the subjected chains. 
    The transition probability for chain ends in the same chain is the same with its opponents, which means we can tweak things by detachment opponent chain ends rather than subjected chain end.
    This is the reason the function is called 'opp_del' rather than 'del'
   */
  MKL_LONG index_target = (MKL_LONG)HASH[index_particle](hash_index_target);
  opp_del_association_grab_IK(index_particle, hash_index_target);
  opp_del_association_IK(index_target, FIND_HASH_INDEX(index_target, index_particle));
  return 0;
}

MKL_LONG ASSOCIATION::opp_del_association(const MKL_LONG index_particle, MKL_LONG index_target)
{
  MKL_LONG hash_index_target = FIND_HASH_INDEX(index_particle, index_target);
  opp_del_association_hash(index_particle, hash_index_target);
  return 0;
}


double KINETICS::CONNECTIVITY_update_CASE_SUGGESTION_particle_hash_target(ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG hash_index_target, double distance)
{
  CONNECT->CASE_SUGGESTION[index_particle](hash_index_target) = POTs.w_function(distance, POTs.f_connector(distance, POTs.force_variables), POTs.force_variables);
  return CONNECT->CASE_SUGGESTION[index_particle](hash_index_target);
}

double KINETICS::CONNECTIVITY_update_CASE_SUGGESTION_particle_target(ASSOCIATION* CONNECT, POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG index_target, double distance)
{
  MKL_LONG hash_index_target = CONNECT->FIND_HASH_INDEX(index_particle, index_target);
  return CONNECTIVITY_update_CASE_SUGGESTION_particle_hash_target(CONNECT, POTs, index_particle, hash_index_target, distance);
}

double KINETICS::CONNECTIVITY_update_Z_SUGGESTION_particle(ASSOCIATION* CONNECT, const MKL_LONG index_particle)
{
  
  // result = dot(x, y)
  CONNECT->Z_SUGGESTION[index_particle] = cblas_ddot(CONNECT->TOKEN[index_particle], // N
                                                     CONNECT->weight[index_particle].data, // x
                                                     1,
                                                     CONNECT->CASE_SUGGESTION[index_particle].data, // y
                                                     1);
  return CONNECT->Z_SUGGESTION[index_particle];
}



double KINETICS::CONNECTIVITY_update_dPDF_SUGGESTION_particle(ASSOCIATION* CONNECT, const MKL_LONG index_particle)
{
  // Note that the element-wise multiplication can be done to make diagonalization. BLAS has good aspect for this functionality.
  // DSBMV which gave us y = alpha*A*x + beta*y,
  // which is designed for symmetric banded matrix for A
  // In this case, we can set the third argument, k, set to 0

  // y = alpha*dot(A, x) + beta*y
  cblas_dsbmv(CblasRowMajor,
              CblasUpper,
              CONNECT->TOKEN[index_particle], // number of accounting
              0, // k is related with band. 0 for diagonal (or super-diagonal)
              1./CONNECT->Z_SUGGESTION[index_particle], // alpha
              CONNECT->weight[index_particle].data, // A (here is array for diagonal components)
              1,  // lda this is leading dimension which inc A. here is set as 1 because the A is only contained diagonal components
              CONNECT->CASE_SUGGESTION[index_particle].data, // x
              1, // incx
              0., // beta. 0 is removing the previous history because 0*y
              CONNECT->dPDF_SUGGESTION[index_particle].data, // y
              1); // incy
  return 0;
}

double KINETICS::CONNECTIVITY_update_dCDF_SUGGESTION_particle(ASSOCIATION* CONNECT, const MKL_LONG index_particle)
{
  CONNECT->dCDF_SUGGESTION[index_particle](0) = CONNECT->dPDF_SUGGESTION[index_particle](0);
  for(MKL_LONG k=1; k<CONNECT->TOKEN[index_particle]; k++)
    {
      CONNECT->dCDF_SUGGESTION[index_particle](k) = CONNECT->dCDF_SUGGESTION[index_particle](k-1) + CONNECT->dPDF_SUGGESTION[index_particle](k);
    }
  return CONNECT->dCDF_SUGGESTION[index_particle](CONNECT->TOKEN[index_particle]);
}


double ASSOCIATION::update_CASE_SUGGESTION_particle_hash_target(POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG hash_index_target, double distance)
{
  return KINETICS::CONNECTIVITY_update_CASE_SUGGESTION_particle_hash_target(this, POTs, index_particle, hash_index_target, distance);
}

double ASSOCIATION::update_CASE_SUGGESTION_particle_target(POTENTIAL_SET& POTs, const MKL_LONG index_particle, MKL_LONG index_target, double distance)
{
  return KINETICS::CONNECTIVITY_update_CASE_SUGGESTION_particle_target(this, POTs, index_particle, index_target, distance);
}
double ASSOCIATION::update_Z_SUGGESTION_particle(const MKL_LONG index_particle)
{
  return KINETICS::CONNECTIVITY_update_Z_SUGGESTION_particle(this, index_particle);
}

double ASSOCIATION::update_dPDF_SUGGESTION_particle(const MKL_LONG index_particle)
{
  return KINETICS::CONNECTIVITY_update_dPDF_SUGGESTION_particle(this, index_particle);
}

double ASSOCIATION::update_dCDF_SUGGESTION_particle(const MKL_LONG index_particle)
{
  return KINETICS::CONNECTIVITY_update_dCDF_SUGGESTION_particle(this, index_particle);
}

double ASSOCIATION::update_CHAIN_SUGGESTION_MAP_particle(const MKL_LONG index_particle, POTENTIAL_SET& POTs, RDIST& R_boost)
{
  double time_st = dsecnd();
  // 1. compute all the case (weight) for pair of particles
  for(MKL_LONG j=0; j<TOKEN[index_particle]; j++)
    {
      // CONNECT.HASH[i](j) gave us the index for target
      // which means we have to compute distance between i and k where k is given by CONNECT.HASH[i](j).
      // CONNECT.update_CASE_SUGGESTION_particle_hash_target(POTs, i, j, R_minimum_distance_boost[i](CONNECT.HASH[i](j))); // RDIS      
      update_CASE_SUGGESTION_particle_hash_target(POTs, index_particle, j, R_boost.Rsca[index_particle](this->HASH[index_particle](j)));
    }
  // 2. compute suggestion partition function
  update_Z_SUGGESTION_particle(index_particle);
  // 3. compute probability
  update_dPDF_SUGGESTION_particle(index_particle);
  // 4. cumulating probability
  update_dCDF_SUGGESTION_particle(index_particle);
  return dsecnd() - time_st;
}

double ASSOCIATION::update_ASSOCIATION_MAP_particle(const MKL_LONG index_particle, POTENTIAL_SET& POTs, RDIST& R_boost)
{
  double time_st = dsecnd();
  MKL_LONG cell_index_particle = R_boost.cell_index[index_particle];
  MKL_LONG count_CDF_TOKEN = 0;
  dCDF_ASSOCIATION[index_particle].set_value(0);
  INDEX_ASSOCIATION[index_particle].set_value(-1);
  // index_set_val(INDEX_ASSOCIATION[index_particle], -1, dCDF_ASSOCIATION[index_particle].size);
  TOKEN_ASSOCIATION[index_particle] = 0;
  for(MKL_LONG k=0; k<R_boost.N_neighbor_cells; k++)
    {
      MKL_LONG cell_index_neighbor = R_boost.NEIGHBOR_CELLS[cell_index_particle][k];
      for(MKL_LONG p=0; p<R_boost.TOKEN[cell_index_neighbor]; p++)
        {
          MKL_LONG index_target = R_boost(cell_index_neighbor, p);
          double distance = R_boost.Rsca[index_particle](index_target);
          INDEX_ASSOCIATION[index_particle](count_CDF_TOKEN) = index_target;
          // INDEX_ASSOCIATION[index_particle][count_CDF_TOKEN] = index_target;
          dCDF_ASSOCIATION[index_particle](count_CDF_TOKEN) = POTs.PDF_connector(distance, POTs.force_variables);
          if(dCDF_ASSOCIATION[index_particle](count_CDF_TOKEN) > 0.)
            {
              TOKEN_ASSOCIATION[index_particle] ++;
            }
          count_CDF_TOKEN++;
        }
    }
  dCDF_ASSOCIATION[index_particle].sort2(INDEX_ASSOCIATION[index_particle]);
  for(MKL_LONG k=Np - TOKEN_ASSOCIATION[index_particle] + 1; k<Np; k++)
    {
      dCDF_ASSOCIATION[index_particle](k) += dCDF_ASSOCIATION[index_particle](k-1);
    }
  for(MKL_LONG k=Np - TOKEN_ASSOCIATION[index_particle]; k<Np; k++)
    {
      dCDF_ASSOCIATION[index_particle](k) /= dCDF_ASSOCIATION[index_particle](Np - 1);
    }
  return dsecnd() - time_st;
}
