
#include "lib_connectivity.h"

CONNECTIVITY::CONNECTIVITY(MKL_LONG number_of_particles, MKL_LONG maximum_connections)
{
  Np = number_of_particles;
  Mc = maximum_connections;
  HASH.initial(Np, maximum_connections, -1);
  for(MKL_LONG i=0; i<Np; i++)
    {
      HASH(i, 0) = i;
    }
  TOKEN.initial(Np, 1, 1);
}
  // TOKEN = (MKL_LONG*) mkl_calloc(Np, sizeof(MKL_LONG), BIT);
  // HASH = (MKL_LONG**) mkl_calloc(Np, sizeof(MKL_LONG), BIT);


MKL_LONG CONNECTIVITY::check_valid(MKL_LONG index_particle, MKL_LONG index_target)
{
  if (index_target > TOKEN(index_particle))
    {
      std::cout << "ERR:: Class CONNECTIVITY:: index_target must larger then TOKEN[index_particle]\n";
    }
  return TRUE;
}
