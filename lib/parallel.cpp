
#include "parallel.h"

LOCK::LOCK()
{
  INITIALIZATION = FALSE;
}

LOCK::LOCK(MKL_LONG N_target)
{
  NL = N_target;
  locker = new bool [NL];
  for(MKL_LONG i=0; i<NL; i++)
    locker[i] = FALSE;
  INITIALIZATION = TRUE;
}

LOCK::~LOCK()
{
  if(INITIALIZATION)
    delete[] locker;
}

MKL_LONG LOCK::RESET()
{
  if(!INITIALIZATION)
    printf("ERR:: RESET LOCK class only possible when it is initialized\n");
  for(MKL_LONG i=0; i<NL; i++)
    locker[i] = FALSE;
  return 0;
}

