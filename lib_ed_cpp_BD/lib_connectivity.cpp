
#include "lib_connectivity.h"

CONNECTIVITY::CONNECTIVITY(COND& given_condition)
{
  Np = atol(given_condition("Np").c_str());
  Mc = 2*atol(given_condition("N_chains_per_particle").c_str()) + atol(given_condition("tolerance_allowing_connections").c_str());
  HASH = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  for(MKL_LONG i=0; i<Np; i++)
    {
      HASH[i].initial(Mc, 1, -1);
    }
  TOKEN = (MKL_LONG*) mkl_malloc(Np*sizeof(MKL_LONG), BIT);
  for(MKL_LONG i=0; i<Np; i++)
    TOKEN[i] = 1;
  if (given_condition("CONTINUATION_CONNECTION")=="TRUE")
    {
      read_exist_hash(given_condition("CONTINUATION_HASH_FN").c_str());
    }
  else
    {
      for(MKL_LONG i=0; i<Np; i++)
        {
          HASH[i](0) = i;
        }
    }
}


CONNECTIVITY::CONNECTIVITY(MKL_LONG number_of_particles, MKL_LONG maximum_connections)
{
  Np = number_of_particles;
  Mc = maximum_connections;
  HASH = (MATRIX*) mkl_malloc(Np*sizeof(MATRIX), BIT);
  TOKEN = (MKL_LONG*) mkl_malloc(Np*sizeof(MKL_LONG), BIT);
  for(MKL_LONG i=0; i<Np; i++)
    {
      HASH[i](0) = i;
      TOKEN[i] = 1;
    }

}

MKL_LONG CONNECTIVITY::read_exist_hash(const char* fn_hash)
{
  ifstream GIVEN_HASH;
  GIVEN_HASH.open(fn_hash);
  MKL_LONG cnt = 0;
  string line;
  while(getline(GIVEN_HASH, line))
    {
      cnt++;
    }
  GIVEN_HASH.clear();
  GIVEN_HASH.seekg(0);
  for(MKL_LONG i=0; i<cnt-Np; i++)
    {
      getline(GIVEN_HASH, line); 
    }
  // flag will set for the last time steps
  // From here, the Np lines will left as initial conditions
  for(MKL_LONG i=0; i<Np; i++)
    {
      MKL_LONG hash_k = 1;
      MKL_LONG count = 0;
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

MKL_LONG CONNECTIVITY::check_valid(MKL_LONG index_particle, MKL_LONG index_target)
{
  if (index_target > TOKEN[index_particle])
    {
      std::cout << "ERR:: Class CONNECTIVITY:: index_target must larger then TOKEN[index_particle]\n";
    }
  return TRUE;
}

