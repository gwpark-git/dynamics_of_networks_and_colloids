#include "cell_list.h"


CLIST::CLIST(COND& given_condition)
{
  /* cut_off_radius = given_cut_off_radius; */
  /* N_dimension = given_N_dimension; */
  /* box_dimension = given_box_dimension; */
  printf("\tCLIST initialization");
  INITIALIZATION = TRUE;
  box_dimension = atof(given_condition("box_dimension").c_str());
  N_dimension = atol(given_condition("N_dimension").c_str());
  cut_off_radius = atof(given_condition("cutoff_connection").c_str());
  Np = atol(given_condition("Np").c_str());
  MAX_IN_CELL = Np; // default
  CELL_LIST_BOOST = FALSE;
  if (given_condition("cell_list") == "TRUE")
    CELL_LIST_BOOST = TRUE;
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

  // related with simple shear
  if(given_condition("SIMPLE_SHEAR") == "TRUE")
    {
      SIMPLE_SHEAR = TRUE;
      shear_axis = atoi(given_condition("shear_axis").c_str());
      shear_grad_axis = atoi(given_condition("shear_grad_axis").c_str());
      map_to_central_box_image = 0.; // started with zero (equilibrium PBC box)
    }

  
  N_div = 1;
  cell_length = box_dimension;
  N_cells = 1;
  N_neighbor_cells = 1;
  if(CELL_LIST_BOOST)
    {
      N_div = (MKL_LONG) box_dimension/cut_off_radius; // floor refine the given divisor
      cell_length = (double) box_dimension/(double)N_div; // real length scale of cells in each axis
      N_cells = (MKL_LONG)pow(N_div, N_dimension); // N_cells = N_div^N_dimension
      N_neighbor_cells = (MKL_LONG)pow(3, N_dimension); // 3 means number for shift factor {-1, 0, +1}
    }
  // printf("given cond: N_div = %ld, cell_length = %ld, N_cells = %ld, N_neighbor_cells = %ld, N_dimension = %ld\n", N_div, cell_length, N_cells, N_neighbor_cells, N_dimension);
  // the following structure will not be affected either the cell list is turned on or off.
  // this is of importance to use the same interface for the main function
  // N_cells = 1;
  // MAX_IN_CELL=400;
  // Np=400;
  // N_neighbor_cells=1;
  printf("\tdynamic allocating CLIST member variables");
  // cell_index = (MKL_LONG*)mkl_malloc(Np*sizeof(MKL_LONG), BIT);
  // TOKEN = (MKL_LONG*)mkl_malloc(N_cells*sizeof(MKL_LONG), BIT);
  // CELL = (MKL_LONG**)mkl_malloc(N_cells*sizeof(MKL_LONG*), BIT);
  // NEIGHBOR_CELLS = (MKL_LONG**)mkl_malloc(N_cells*sizeof(MKL_LONG*), BIT);
  // BEYOND_BOX = (MKL_LONG***)mkl_malloc(N_cells*sizeof(MKL_LONG**), BIT);

  cell_index = new MKL_LONG [Np];
  for(MKL_LONG i=0; i<Np; i++)
    cell_index[i] = -1;

  TOKEN = new MKL_LONG [N_cells];
  CELL = new MKL_LONG* [N_cells];
  NEIGHBOR_CELLS = new MKL_LONG* [N_cells];
  BEYOND_BOX = new MKL_LONG** [N_cells];
  
  for(MKL_LONG i=0; i<N_cells; i++)
    {
      // printf("%d\n", i);
      // CELL[i].initial(MAX_IN_CELL, 1, 0.);
      TOKEN[i] = -1; // -1 is default value when there is no index (note that the particle index is started with 0, which is the reason to avoid using 0 identifier
      // CELL[i] = (MKL_LONG*)mkl_malloc(MAX_IN_CELL*sizeof(MKL_LONG), BIT);
      CELL[i] = new MKL_LONG [MAX_IN_CELL];
      for(MKL_LONG j=0; j<MAX_IN_CELL; j++)
        CELL[i][j] = -1;
      // NEIGHBOR_CELLS[i] = (MKL_LONG*)mkl_malloc(N_neighbor_cells*sizeof(MKL_LONG), BIT);
      // BEYOND_BOX[i] = (MKL_LONG**)mkl_malloc(N_neighbor_cells*sizeof(MKL_LONG*), BIT);
      NEIGHBOR_CELLS[i] = new MKL_LONG [N_neighbor_cells];
      BEYOND_BOX[i] = new MKL_LONG* [N_neighbor_cells];
      for(MKL_LONG j=0; j<N_neighbor_cells; j++)
        {
          // BEYOND_BOX[i][j] = (MKL_LONG*)mkl_malloc(N_dimension*sizeof(MKL_LONG), BIT);
          BEYOND_BOX[i][j] = new MKL_LONG [N_dimension];
          for(MKL_LONG k=0; k<N_dimension; k++)
            BEYOND_BOX[i][j][k] = 0;

          NEIGHBOR_CELLS[i][j] = -1;
        }
      /* get_neighbor_cell_list(i,  */
    }
  printf("\tallocating neighbor cell list information");
  allocate_index_neighbor_cell_list(); // this is safety for the flag, CELL_LIST_BOOST
  printf("\tDONE\n");
}



MKL_LONG CLIST::identify_cell_from_given_position(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG index_particle, MKL_LONG *index_vec_boost)
{
  MKL_LONG re = 0;  
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      index_vec_boost[k] = (MKL_LONG)(TRAJ(index_t_now, index_particle, k)/cell_length);
    }
  index_vec2sca(index_vec_boost, re);
  return re;
}

MKL_LONG CLIST::identify_cell_from_given_position(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t_now, MKL_LONG index_particle, MKL_LONG *index_vec_boost)
{
  MKL_LONG re = 0;  
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      index_vec_boost[k] =(MKL_LONG)(TRAJ.data[index_t_now][index_particle][k]/cell_length);
    }
  index_vec2sca(index_vec_boost, re);
  return re;
}


MKL_LONG CLIST::allocate_cells_from_positions(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG *index_vec_boost)
{
  for(MKL_LONG i=0; i<N_cells; i++)
    TOKEN[i] = 0;
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG index = identify_cell_from_given_position(TRAJ, index_t_now, i, index_vec_boost);
      // if(index > N_cells)
      // 	printf("cell_index[%d] = %d with (%d, %d, %d)\n", i, index, index_vec_boost[0], index_vec_boost[1], index_vec_boost[2]);
      cell_index[i] = index;
      CELL[index][TOKEN[index]++] = i;
      // CELL[index](TOKEN[index]++) = i;
    }
  return 0;
}

MKL_LONG CLIST::allocate_cells_from_positions(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t_now, MKL_LONG *index_vec_boost)
{
  for(MKL_LONG i=0; i<N_cells; i++)
    TOKEN[i] = 0;
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG index = identify_cell_from_given_position(TRAJ, index_t_now, i, index_vec_boost);
      cell_index[i] = index;
      CELL[index][TOKEN[index]++] = i;
      // printf("cell_index[%d] = %d\n", i, index);
      // CELL[index](TOKEN[index]++) = i;
    }
  return 0;
}

// TRAJ(index_t_now, 
// for(MKL_LONG i=0; i<TRAJ.Np; i++)
//   {
//     for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
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
      re += index_vec[n]*(MKL_LONG)pow(N_div, N_dimension - (n+1));
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
      index_vec[n] = ((MKL_LONG)(index_sca/pow(N_div, N_dimension - (n+1))))%N_div;
    }
  return index_sca;
}

MKL_LONG CLIST::allocate_index_neighbor_cell_list()
{
  if(CELL_LIST_BOOST)
    {
      // MKL_LONG* index_vec_boost = (MKL_LONG*)mkl_malloc(N_dimension*sizeof(MKL_LONG), BIT);
      // MKL_LONG* sf_vec_boost = (MKL_LONG*)mkl_malloc(3*sizeof(MKL_LONG), BIT); // 3 means number of shift factors {-1, 0, +1}
      MKL_LONG* index_vec_boost = new MKL_LONG [N_dimension];
      MKL_LONG* sf_vec_boost = new MKL_LONG [3];
      // for(MKL_LONG i=0; i<N_neighbor_cells; i++) // this was bug.
      for(MKL_LONG i=0; i<N_cells; i++) // this is right. 
        {
          // note that the number of cells is N_div^N_dimension while number of neighbor cells is N_sf^N_dimension
          // for instance, if the PBC is divided by 5 cells for each axis, N_cells = 5^3 = 125 while the N_neighbor_cells = 3^3 = 9
          get_neighbor_cell_list(i, NEIGHBOR_CELLS[i], index_vec_boost, sf_vec_boost);
        }
      // mkl_free(index_vec_boost);
      // mkl_free(sf_vec_boost);
      delete[] index_vec_boost;
      delete[] sf_vec_boost;
    }
  else
    {
      // it prevent to violate the interface
      NEIGHBOR_CELLS[0][0] = 0;
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
  // MKL_LONG shift_factors[3] = {-1, 0, 1};
  MKL_LONG N_sf = 3; // {-1, 0, +1}
  // printf("=========\n");
  
  for(MKL_LONG nsf=0; nsf<pow(N_sf, N_dimension); nsf++)
    {
      UTILITY::index_sca2vec(nsf, sf_vec_boost, N_dimension, N_sf);
      for(MKL_LONG n=0; n<N_dimension; n++)
        {
          sf_vec_boost[n] += self_index_vec_boost[n] - 1;
          // the following function related with the periodic boundary condition
          // whenever the neighbor beyond the PBC box, it re-mapped inside of the box
          if(sf_vec_boost[n] < 0)
            {
              // this will reduce the overhead during relative distance between particle with PBC boundary condition	      
              BEYOND_BOX[index_sca][nsf][n] = -1;  // when the neighbor cell is beyond PBC boundary in left side
              sf_vec_boost[n] += N_div;
            }
          else if(sf_vec_boost[n] >= N_div)
            {
              BEYOND_BOX[index_sca][nsf][n] = +1; // when the neighbor cell is beyond PBC boundary in right side
              sf_vec_boost[n] -= N_div;
            }
          
          // sf_vec_boost[n] has the values in [0, 1, 2] since the N_sf is set with 3.
          // sf_vec_boost[n] -1 becomes [-1, 0, +1] which is exactly the same with shift factor array
          // self_index_vec[n] + (sf_vec_boost[n] - 1) -> sf_vec_boost[n]
        }
      // UTILITY::index_vec2sca(sf_vec_boost, index_neighbor_cells[nsf], N_dimension, N_sf); // this is the big bug
      UTILITY::index_vec2sca(sf_vec_boost, index_neighbor_cells[nsf], N_dimension, N_div);
      // printf("self=%ld(%ld, %ld, %ld), cell[%ld,%ld,%ld]=%ld\n", index_sca, self_index_vec_boost[0], self_index_vec_boost[1], self_index_vec_boost[2], sf_vec_boost[0], sf_vec_boost[1], sf_vec_boost[2], index_neighbor_cells[nsf]);
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


