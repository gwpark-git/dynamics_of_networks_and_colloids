#include "cell_list.h"

CLIST::CLIST(COND& given_condition)
{
  /* cut_off_radius = given_cut_off_radius; */
  /* N_dimension = given_N_dimension; */
  /* box_dimension = given_box_dimension; */
  printf("\tCLIST initialization");
  box_dimension = atof(given_condition("box_dimension").c_str());
  N_dimension = atoi(given_condition("N_dimension").c_str());
  cut_off_radius = atof(given_condition("cutoff_connection").c_str());
  Np = atoi(given_condition("Np").c_str());
  MAX_IN_CELL = Np; // default
  if (given_condition("cell_list") == "TRUE")
    CELL_LIST_BOOST = TRUE;
  else
    CELL_LIST_BOOST = FALSE;
  printf("\tcheck cut-off scheme");
  
  if((cut_off_radius <= 0 || cut_off_radius > box_dimension/3.) || (!CELL_LIST_BOOST))
    {
      // it prevent to use wrong value of cut_off_radius
      // when it is set with 0 or -1, the cut_off will not be applied
      // when the cut_off sectionize the box as number 3, the cell_list is not benefit to use it
      // at least, the number of cells per axis should be 4, so we have benefit to compute
      // note that it is of importance to prevent the cell_list functionality is turned on/off.
      cut_off_radius = box_dimension;
      if(CELL_LIST_BOOST)
        {
          CELL_LIST_BOOST = FALSE; 
          printf("\nWARNING:cell_list=TRUE but the cut_off_radius is given by %4.3f while box_dimension is %4.3f\n", cut_off_radius, box_dimension);
          printf("\tcell list is turned off\n");
        }
    }
      
  if(CELL_LIST_BOOST)
    {
      N_div = (MKL_LONG) box_dimension/cut_off_radius; // floor refine the given divisor
      cell_length = (double) box_dimension/(double)N_div; // real length scale of cells in each axis
      N_cells = (MKL_LONG)pow(N_div, N_dimension); // N_cells = N_div^N_dimension
      N_neighbor_cells = pow(3, N_dimension); // 3 means number for shift factor {-1, 0, +1}
    }
  else // in the case for cell list is turned off
    {
      N_div = 1;
      cell_length = box_dimension;
      N_cells = 1;
      N_neighbor_cells = 1;
    }

  // the following structure will not be affected either the cell list is turned on or off.
  // this is of importance to use the same interface for the main function
  // N_cells = 1;
  // MAX_IN_CELL=400;
  // Np=400;
  // N_neighbor_cells=1;
  printf("\tdynamic allocating CLIST member variables");
  TOKEN = (MKL_LONG*)mkl_malloc(N_cells*sizeof(MKL_LONG), BIT);
  CELL = (MATRIX*)mkl_malloc(N_cells*sizeof(MATRIX), BIT);
  NEIGHBOR_CELLS = (MKL_LONG**)mkl_malloc(N_cells*sizeof(MKL_LONG*), BIT);
  for(MKL_LONG i=0; i<N_cells; i++)
    {
      CELL[i].initial(MAX_IN_CELL, 1, 0.);
      TOKEN[i] = -1; // -1 is default value when there is no index (note that the particle index is started with 0, which is the reason to avoid using 0 identifier
      NEIGHBOR_CELLS[i] = (MKL_LONG*)mkl_malloc(N_neighbor_cells*sizeof(MKL_LONG), BIT);
      /* get_neighbor_cell_list(i,  */
    }
  printf("\tallocating neighbor cell list information");
  allocate_index_neighbor_cell_list(); // this is safety for the flag, CELL_LIST_BOOST
  printf("\tDONE\n");
}



MKL_LONG CLIST::identify_cell_from_given_position(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG index_particle, MKL_LONG *index_vec_boost)
{
  MKL_LONG re = 0;  
  for(MKL_LONG k=0; k<TRAJ.dimension; k++)
    {
      index_vec_boost[k] = (MKL_LONG)(TRAJ(index_t_now, index_particle, k)/cell_length);
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
  if(CELL_LIST_BOOST)
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
  else
    {
      // it prevent to violate the interface
      NEIGHBOR_CELLS[0] = 0;
    }
  return 0;
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
      UTILITY::index_sca2vec(nsf, sf_vec_boost, N_dimension, N_sf);
      for(MKL_LONG n=0; n<N_dimension; n++)
        {
          sf_vec_boost[n] += self_index_vec_boost[n] - 1;
          // the following function related with the periodic boundary condition
          // whenever the neighbor beyond the PBC box, it re-mapped inside of the box
          if(sf_vec_boost[n] < 0)
            sf_vec_boost[n] += N_div;
          else if(sf_vec_boost[n] >= N_div)
            sf_vec_boost[n] -= N_div;
          
          // sf_vec_boost[n] has the values in [0, 1, 2] since the N_sf is set with 3.
          // sf_vec_boost[n] -1 becomes [-1, 0, +1] which is exactly the same with shift factor array
          // self_index_vec[n] + (sf_vec_boost[n] - 1) -> sf_vec_boost[n]
        }
      UTILITY::index_vec2sca(sf_vec_boost, index_neighbor_cells[nsf], N_dimension, N_sf);
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



// MKL_LONG CLIST::MAPPING(TRAJECTORY& TRAJ)
// {
//   // TODO parallizable mapping function
//   // Quick in mind: the sequencial mapping function is not enough to improve the speed
//   // OR: it is just called once when the initial trajectory is given
//   //     On this regards, the update during simulation should not use this MAPPING function since it will overhead because parallelization is not implemented in this function.
//   return 0;
// }

