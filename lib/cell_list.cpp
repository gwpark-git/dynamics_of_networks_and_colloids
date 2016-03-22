#include "cell_list.h"

MKL_LONG CLIST::identify_cell_from_given_position(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG index_particle, MKL_LONG *index_vec_boost)
{
  MKL_LONG re = 0;  
  for(MKL_LONG k=0; k<TRAJ.dimension; k++)
    {
      index_vec_boost[k] = (MKL_LONG)(TRAJ(index_t_now, i, k)/cell_length);
    }
  index_vec2sca(index_vec_boost, re);
  return re;
}

MKL_LONG CLIST::allocate_cells_from_positions(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG *index_vec_boost)
{
  MKL_LONG index = 0;
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      index = identify_cell_from_given_position(TRAJ, index_t_now, i, index_vec_boost);
      CELL[index](TOKEN[index]++) = i;
    }
  return 0;
}
  // TRAJ(index_t_now, 
  // for(MKL_LONG i=0; i<TRAJ.Np; i++)
  //   {
  //     for(MKL_LONG k=0; k<TRAJ.dimension; k++)
  //       {
  //         index_vec_boost[i][k] = (MKL_LONG)(TRAJ(index_t_now, i, k)/cell_length);
  //       }

  //   }
//   return 0;
// }

MKL_LONG UTILITY::index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca, const MKL_LONG N_dimension, const MKL_LONG N_div)
{
  /*
    This is the original index mapping function for general N-dimensional case. 
   */
  MKL_LONG re = 0;
  for(MKL_LONG n=0; n<N_dimension; n++)
    {
      re += index_vec[n]*pow(N_div, N_dimension - (n+1));
    }
  index_sca = re;
  return re;
}

MKL_LONG UTILITY::index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec, const MKL_LONG N_dimension, const MKL_LONG N_div)
{
  /*
    This inverse map applied to 
   */
  for(MKL_LONG n=0; n<N_dimension; n++)
    {
      index_vec[n] = ((MKL_LONG)(index_sca/pow(N_div, N_dimension - (n+2))))%N_div;
    }
  return index_sca;
}

MKL_LONG CLIST::allocate_index_neighbor_cell_list()
{
  MKL_LONG* index_vec_boost = (MKL_LONG*)mkl_malloc(N_dimension*sizeof(MKL_LONG), BIT);
  MKL_LONG* sf_vec_boost = (MKL_LONG*)mkl_malloc(3*sizeof(MKL_LONG), BIT); // 3 means number of shift factors {-1, 0, +1}
  for(MKL_LONG i=0; i<N_neighbor_cells; i++)
    {
      get_neighbor_cell_list(i, NEIGHBOR_CELLS[i], index_vec_boost, sf_vec_boost);
    }
  mkl_free(index_vec_boost);
  mkl_free(sf_vec_boost);
}

MKL_LONG CLIST::get_neighbor_cell_list(const MKL_LONG& index_sca, MKL_LONG* index_neighbor_cells, MKL_LONG* self_index_vec_boost, MKL_LONG* sf_vec_boost)
{
  /*
    get_neighbor_cell_list will store the index of neighbor_list into the index_neighbor_cells based on the scalar values.
    However, to call this function during the simulation make big overhead to compute index sets.
    Hence, we need to generate the array for the list when it is initialized, then re-used the generated list during simulation.
    For non-equilibrium simulation such as simple shear flow, the re-usability will decrease which eventually adds some overhead to compute additional lists.
   */
  CLIST::index_sca2vec(index_sca, self_index_vec_boost);
  MKL_LONG shift_factors[3] = {-1, 0, 1};
  MKL_LONG N_sf = 3; // {-1, 0, +1}
  for(MKL_LONG nsf=0; nsf<pow(N_sf, N_dimension); nsf++)
    {
      UTILITY::index_sca2vec(nsf, sf_index, N_dimension, N_sf);
      for(MKL_LONG n=0; n<N_dimension; n++)
        {
          sf_index[n] += self_index_vec[n] - 1;
          // the following function related with the periodic boundary condition
          // whenever the neighbor beyond the PBC box, it re-mapped inside of the box
          if(sf_index[n] < 0)
            sf_index[n] += N_div;
          else if(sf_index[n] >= N_div)
            sf_index[n] -= N_div;
          
          // sf_index[n] has the values in [0, 1, 2] since the N_sf is set with 3.
          // sf_index[n] -1 becomes [-1, 0, +1] which is exactly the same with shift factor array
          // self_index_vec[n] + (sf_index[n] - 1) -> sf_index[n]
        }
      UTILITY::index_vec2sca(sf_index, index_neighbor_cells[nsf], N_dimension, N_sf);
    }
  return 0;
}

MKL_LONG CLIST::index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca)
{
  return UTILITY::index_vec2sca(index_vec, index_sca, N_dimension, N_div);
}

MKL_LONG CLIST::index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec)
{
  return UTILITY::index_sca2vec(index_sca, index_vec, N_dimension, N_div);
}


// MKL_LONG CLIST::index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca)
// {
//   /*
//    (basic) index_vec = [i, j, k] => i*N_div*N_div + j*N_div + k
//    This is basically the same with
//    index_vec[0]*pow(N_div, 2) + index_vec[1]*pow(N_div, 1) + index_vec[2]*pow(N_div, 0)

//    In this approach, the neighbour cells can be checked by shift factor [-1, 0, +1] for each components for index_vec.
//   */
//   // MKL_LONG re = 0;
//   // for(MKL_LONG n=0; n< N_dimension; n++)
//   //   {
//   //     re += index_vec[n]*pow(N_div, 2-n);
//   //   }
//   // MKL_LONG re = 0;
//   // for (MKL_LONG n=0; n<N_dimension; n++)
//   //   P{

      
//   return re;
// }

// MKL_LONG CLIST::index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec)
// {
//   /*
//     This function will recover index set for each axis (inverse map of the index_vec2sca).
//     Note that the com(index_vec2sca, index_sca2vec) is identity function.
//    */
//   index_vec[2] = (MKL_LONG)index_sca%N_div;
//   index_vec[1] = (MKL_LONG)index_sca/(N_div);
//   index_vec[0] = (MKL_LONG)index_sca/(N_div*N_div);
//   for(MKL_LONG n=0; n<N_dimension-1; n++)
//     {
//       index_vec[n] = (MKL_LONG)index_sca/pow(N_div, 2-n);
//     }
//   index_vec[N_dimension-1] = (MKL_LONG)index_sca%N_div;
// }



MKL_LONG CLIST::MAPPING(TRAJECTORY& TRAJ)
{
  // TODO parallizable mapping function
  // Quick in mind: the sequencial mapping function is not enough to improve the speed
  // OR: it is just called once when the initial trajectory is given
  //     On this regards, the update during simulation should not use this MAPPING function since it will overhead because parallelization is not implemented in this function.
  return 0;
}

