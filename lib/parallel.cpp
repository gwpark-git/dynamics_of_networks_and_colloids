
#include "parallel.h"

LOCK::LOCK()
{
  INITIALIZATION = FALSE;
}

LOCK::LOCK(MKL_LONG N_target)
{
  NL = N_target;
  // locker = (bool*) mkl_calloc(NL, sizeof(bool), BIT);
  locker = new bool [NL];
  for(MKL_LONG i=0; i<NL; i++)
    locker[i] = FALSE;
  INITIALIZATION = TRUE;
}

LOCK::~LOCK()
{
  if(INITIALIZATION)
    delete[] locker;
    // mkl_free(locker);
}

MKL_LONG LOCK::RESET()
{
  if(!INITIALIZATION)
    printf("ERR:: RESET LOCK class only possible when it is initialized\n");
  for(MKL_LONG i=0; i<NL; i++)
    locker[i] = FALSE;
  return 0;
}


// bool LOCK::CHECKER(MKL_LONG const index_target)
// {
//   if(locker[index_target])
//     return TRUE;
//   return FALSE;
// }

// bool& LOCK::operator()(MKL_LONG const index_target)
// {
//   return locker[index_target];
// }
