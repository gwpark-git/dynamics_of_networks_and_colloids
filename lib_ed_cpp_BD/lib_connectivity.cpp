
#include "lib_connectivity.h"

CONNECTIVITY::CONNECTIVITY(COND& given_condition)
{
  Np = atol(given_condition("Np").c_str());
  Mc = 2*atol(given_condition("N_chains_per_particle").c_str()) + atol(given_condition("tolerance_allowing_connections").c_str());
  // HASH.initial(Np, Mc, -1);
  // TOKEN.initial(Np, 1, 1);
  // HASH = (MATRIX_LONG*) mkl_malloc(Np*sizeof(MATRIX_LONG), BIT);
  HASH = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  for(long i=0; i<Np; i++)
    {
      HASH[i].initial(Mc, 1, -1);
    }
  TOKEN = (long*) mkl_malloc(Np*sizeof(long), BIT);
  for(long i=0; i<Np; i++)
    TOKEN[i] = 1;
  if (given_condition("CONTINUATION_CONNECTION")=="TRUE")
    {
      read_exist_hash(given_condition("CONTINUATION_HASH_FN").c_str());
    }
  else
    {
      for(long i=0; i<Np; i++)
        {
          HASH[i](0) = i;
        }
    }
}


CONNECTIVITY::CONNECTIVITY(long number_of_particles, long maximum_connections)
{
  Np = number_of_particles;
  Mc = maximum_connections;
  // HASH = (MATRIX_LONG*) mkl_malloc(Np*sizeof(MATRIX_LONG), BIT);
  HASH = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  TOKEN = (long*) mkl_malloc(Np*sizeof(long), BIT);
  for(long i=0; i<Np; i++)
    {
      HASH[i](0) = i;
      TOKEN[i] = 1;
    }

  // HASH.initial(Np, maximum_connections, -1);
  // for(long i=0; i<Np; i++)
  //   {
  //     HASH(i, 0) = i;
  //   }
  // TOKEN.initial(Np, 1, 1);
}
  // TOKEN = (long*) mkl_calloc(Np, sizeof(long), BIT);
  // HASH = (long**) mkl_calloc(Np, sizeof(long), BIT);

long CONNECTIVITY::read_exist_hash(const char* fn_hash)
{
  ifstream GIVEN_HASH;
  GIVEN_HASH.open(fn_hash);
  long cnt = 0;
  string line;
  while(getline(GIVEN_HASH, line))
    {
      cnt++;
    }
  GIVEN_HASH.clear();
  GIVEN_HASH.seekg(0);
  for(long i=0; i<cnt-Np; i++)
    {
      getline(GIVEN_HASH, line); 
    }
  // flag will set for the last time steps
  // From here, the Np lines will left as initial conditions
  for(long i=0; i<Np; i++)
    {
      long hash_k = 1;
      long count = 0;
      while(hash_k != -1) // note that the initialization condition for hash is -1
        {
          GIVEN_HASH >> hash_k;
          // GIVEN_HASH >> HASH(++count);//hash_k;
          HASH[i](++count) = hash_k;
        }
      // TOKEN(i) = count;
      TOKEN[i] = count;
    }
  GIVEN_HASH.close();
  return 0;
}

long CONNECTIVITY::check_valid(long index_particle, long index_target)
{
  if (index_target > TOKEN[index_particle])
    {
      std::cout << "ERR:: Class CONNECTIVITY:: index_target must larger then TOKEN[index_particle]\n";
    }
  return TRUE;
}

