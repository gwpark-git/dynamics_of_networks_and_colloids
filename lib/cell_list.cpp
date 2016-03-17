
#include "cell_list.h"

MKL_LONG CLIST::cell_index(MKL_LONG* index_vec)
{
  /*
   (basic) index_vec = [i, j, k] => i*N_div*N_div + j*N_div + k
   This is basically the same with
   index_vec[0]*pow(N_div, 2) + index_vec[1]*pow(N_div, 1) + index_vec[2]*pow(N_div, 0)

   In this approach, the neighbour cells can be checked by shift factor [-1, 0, +1] for each components for index_vec.
  */
  MKL_LONG re = 0;
  for(MKL_LONG n=0; n< N_dimension; n++)
    {
      re += index_vec[n]*pow(N_div, 2-n);
    }
  return re;
}

MKL_LONG CLIST::MAPPING(TRAJECTORY& TRAJ)
{
  // TODO parallizable mapping function
  // Quick in mind: the sequencial mapping function is not enough to improve the speed
  // OR: it is just called once when the initial trajectory is given
  //     On this regards, the update during simulation should not use this MAPPING function since it will overhead because parallelization is not implemented in this function.
  return 0;
}

