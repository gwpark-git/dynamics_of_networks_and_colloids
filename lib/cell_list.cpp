#include "cell_list.h"
#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

CLIST::
CLIST(COND& given_condition)
{
  printf("\tCLIST initialization");
  INITIALIZATION = TRUE;
  N_dimension = atol(given_condition("N_dimension").c_str());

  box_dimension = new double [N_dimension];
  double min_box_dimension = atof(given_condition("box_dimension").c_str());
  double max_box_dimension = atof(given_condition("box_dimension").c_str());
  for (MKL_LONG k=0; k<N_dimension; k++)
    {
      // box_dimension[k] = atof(given_condition("box_dimension").c_str());
      box_dimension[k] = HANDLE_COND::get_LBk(given_condition, k);
      min_box_dimension = min(min_box_dimension, box_dimension[k]);
      max_box_dimension = max(max_box_dimension, box_dimension[k]);
    }
  if(min_box_dimension != max_box_dimension)
    {
      printf("\nAsymmetric PBC box is setted: ");
      for(MKL_LONG k=0; k<N_dimension; k++)
        {
          printf("%lf ", box_dimension[k]);
        }
      printf("\n");
    }  
  // if (atof(given_condition("LBx").c_str()) > atof(given_condition("box_dimension").c_str())):
  //   box_dimension[0] = atof(given_condition("LBx").c_str());
  // if (atof(given_condition("LBy").c_str()) > atof(given_condition("box_dimension").c_str())):
  //   box_dimension[1] = atof(given_condition("LBy").c_str());
  // if (atof(given_condition("LBz").c_str()) > atof(given_condition("box_dimension").c_str())):
  //   box_dimension[2] = atof(given_condition("LBz").c_str());
  
  
  cut_off_radius = atof(given_condition("cutoff_connection").c_str());
  Np = atol(given_condition("Np").c_str());
  MAX_IN_CELL = Np; // default
  CELL_LIST_BOOST = FALSE;
  if (given_condition("cell_list") == "TRUE")
    CELL_LIST_BOOST = TRUE;
  printf("\tcheck cut-off scheme");

  if((cut_off_radius <= 0 || cut_off_radius > min_box_dimension/3.) || (!CELL_LIST_BOOST))
    {
      // it prevent to use wrong value of cut_off_radius
      // when it is set with 0 or -1, the cut_off will not be applied
      // when the cut_off sectionize the box as number 3, the cell_list is not benefit to use it
      // at least, the number of cells per axis should be 4, so we have benefit to compute
      // note that it is of importance to prevent the cell_list functionality is turned on/off.
      cut_off_radius = min_box_dimension;
      if(CELL_LIST_BOOST)
        {
          CELL_LIST_BOOST = FALSE; 
          printf("\nWARNING:cell_list=TRUE but the cut_off_radius is given by %4.3f while min_box_dimension is %4.3f\n", cut_off_radius, min_box_dimension);
          printf("\tcell list is turned off\n");
        }
    }


  
  // N_div = 1;
  // cell_length = box_dimension;
  cell_length = new double [N_dimension];
  // N_cells = new MKL_LONG [N_dimension];
  for (MKL_LONG k=0; k<N_dimension; k++)
    {
      cell_length[k] = atof(given_condition("box_dimension").c_str());
      N_cells = 1;
    }
  N_neighbor_cells = 1;
  N_div = new MKL_LONG [N_dimension];
  for (MKL_LONG k=0; k<N_dimension; k++)
    N_div[k] = 1;
  N_cells = 1;
  if(CELL_LIST_BOOST)
    {
      for(MKL_LONG k=0; k<N_dimension; k++)
        {
          N_div[k] = (MKL_LONG) box_dimension[k]/cut_off_radius; // floor refine the given divisor
          cell_length[k] = (double) box_dimension[k]/(double)N_div[k]; // real length scale of cells in each axis
          N_cells *= (MKL_LONG)N_div[k];
          // N_cells[k] = (MKL_LONG)pow(N_div[k], N_dimension); // N_cells = N_div^N_dimension
        }
      N_neighbor_cells = (MKL_LONG)pow(3, N_dimension); // 3 means number for shift factor {-1, 0, +1}
    }
  printf("\tdynamic allocating CLIST member variables");

  cell_index = new MKL_LONG [Np];
  for(MKL_LONG i=0; i<Np; i++)
    cell_index[i] = -1;

  TOKEN = new MKL_LONG [N_cells];
  CELL = new MKL_LONG* [N_cells];
  NEIGHBOR_CELLS = new MKL_LONG* [N_cells];
  BEYOND_BOX = new MKL_LONG** [N_cells];
  
  for(MKL_LONG i=0; i<N_cells; i++)
    {
      TOKEN[i] = -1; // -1 is default value when there is no index (note that the particle index is started with 0, which is the reason to avoid using 0 identifier
      CELL[i] = new MKL_LONG [MAX_IN_CELL];
      for(MKL_LONG j=0; j<MAX_IN_CELL; j++)
        CELL[i][j] = -1;
      NEIGHBOR_CELLS[i] = new MKL_LONG [N_neighbor_cells];
      BEYOND_BOX[i] = new MKL_LONG* [N_neighbor_cells];
      for(MKL_LONG j=0; j<N_neighbor_cells; j++)
        {
          BEYOND_BOX[i][j] = new MKL_LONG [N_dimension];
          for(MKL_LONG k=0; k<N_dimension; k++)
            BEYOND_BOX[i][j][k] = 0;

          NEIGHBOR_CELLS[i][j] = -1;
        }
    }

  // related with simple shear
  SIMPLE_SHEAR = FALSE;
  STEP_SHEAR = FALSE;
  // if(given_condition("STEP_SHEAR") == "TRUE")
  if(given_condition("MECHANICAL_PERTURBATION") == "STEP_SHEAR")
    {
      STEP_SHEAR = TRUE;
      shear_axis = atoi(given_condition("shear_axis").c_str());
      shear_grad_axis = atoi(given_condition("shear_grad_axis").c_str());
      map_to_central_box_image = 0.; // started with zero (equilibrium PBC box)

    }
  // NEIGHBOR_CELLS_OFFSET = NULL;
  // if(given_condition("SIMPLE_SHEAR") == "TRUE")
  if(given_condition("MECHANICAL_PERTURBATION") == "SIMPLE_SHEAR")
    {
      SIMPLE_SHEAR = TRUE;
      shear_axis = atoi(given_condition("shear_axis").c_str());
      shear_grad_axis = atoi(given_condition("shear_grad_axis").c_str());
      map_to_central_box_image = 0.; // started with zero (equilibrium PBC box)

    //   if(CELL_LIST_BOOST)
	// {

	//   NEIGHBOR_CELLS_OFFSET = new MKL_LONG** [2]; // it is related with left and right boundaries
	  
	//   /*
	//     It is of important to ware that the NEIGHBOR_CELLS_OFFSET have "N_cells + 1" in rows.
	//     This is due to the fact that the value of map_to_central_box_image covers -0.5 box_dimension to 0.5 box_dimension, which support index function from "-0.5 box_dimension mod L_c" to "0.5 box_dimension mod L_c". Note that it gave us "N_cells + 1". For instance, box_dimension is 10 while L_c = 2. In this case, -5 into index of cell as -2 and +5 - delta into index of cell as 2. Therefore, we have 6 related cells: -2, -3, -4, 0, +1, +2 (minus signs are index for left-box).
	//    */	  

	//   /*
	//     Number of offset is related with how many neighbor cells are related with the given index function.
	//     For 2-dimensional case, the offset is obviously 2 because of +1 depth for searching cells both of left and right in the box in beyond shear_grad_axis.
	//     For 3-dimensional case, the offset is 6 since it is possibly three through on the third aixs: cross of shear and shear_grad axes.
	//    */
	  
	//   N_offset = 2;
	//   if(N_dimension == 3)
	//     N_offset = 2*3;

	  
	//   for(MKL_LONG i=0; i<2; i++)
	//     {
	//       NEIGHBOR_CELLS_OFFSET[i] = new MKL_LONG* [N_cells + 1];

	  
	//       for(MKL_LONG i=0; i<N_cells + N_offset; i++)
	// 	{
	// 	  /*
	// 	    Note that all the possible index for neighbor_cells even for offset is pre-determined rather than generated during simulation in order to remove overhead.
	// 	    So, the first dimension, "N_cells + 1", is the seed to recover possible neighbor cells in this case. 
	// 	  */
	      
	// 	  NEIGHBOR_CELLS_OFFSET[i] = new MKL_LONG [N_neighbor_cells + N_offset];
	// 	}
	//     }
	// }
      // note that 2*N_dimension is only valid for 3-dimensional space
      
    }

  
  printf("\tallocating neighbor cell list information");
  allocate_index_neighbor_cell_list(); // this is safety for the flag, CELL_LIST_BOOST
  printf("\tDONE\n");
}



MKL_LONG
CLIST::
identify_cell_from_given_position(TRAJECTORY& TRAJ, MKL_LONG index_t_now,
				  MKL_LONG index_particle, MKL_LONG *index_vec_boost)
{
  MKL_LONG re = 0;  
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      index_vec_boost[k] = (MKL_LONG)(TRAJ(index_t_now, index_particle, k)/cell_length[k]);
    }
  index_vec2sca(index_vec_boost, re);
  return re;
}

MKL_LONG
CLIST::
identify_cell_from_given_position(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t_now,
				  MKL_LONG index_particle, MKL_LONG *index_vec_boost)
{
  MKL_LONG re = 0;  
  for(MKL_LONG k=0; k<TRAJ.N_dimension; k++)
    {
      index_vec_boost[k] =(MKL_LONG)(TRAJ.data[index_t_now][index_particle][k]/cell_length[k]);
    }
  index_vec2sca(index_vec_boost, re);
  return re;
}


MKL_LONG
CLIST::
allocate_cells_from_positions(TRAJECTORY& TRAJ, MKL_LONG index_t_now,
			      MKL_LONG *index_vec_boost)
{
  for(MKL_LONG i=0; i<N_cells; i++)
    TOKEN[i] = 0;
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG index = identify_cell_from_given_position(TRAJ, index_t_now, i, index_vec_boost);
      cell_index[i] = index;
      CELL[index][TOKEN[index]++] = i;
    }
  return 0;
}

MKL_LONG
CLIST::
allocate_cells_from_positions(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t_now,
			      MKL_LONG *index_vec_boost)
{
  for(MKL_LONG i=0; i<N_cells; i++)
    TOKEN[i] = 0;
  for(MKL_LONG i=0; i<TRAJ.Np; i++)
    {
      MKL_LONG index = identify_cell_from_given_position(TRAJ, index_t_now, i, index_vec_boost);
      cell_index[i] = index;
      CELL[index][TOKEN[index]++] = i;
    }
  return 0;
}

MKL_LONG
UTILITY::
index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca,
	      const MKL_LONG N_dimension, const MKL_LONG *N_div)
{
  /*
    This is the original index mapping function for general N-dimensional case. 
  */
  MKL_LONG re = 0;
  for(MKL_LONG n=0; n<N_dimension; n++)
    {
      /*
	must be checked on here
      */
      double multiplier = 1.;
      for(MKL_LONG k=n + 1; k<N_dimension; k++)
        {
          multiplier *= N_div[k];
        }
      re += index_vec[n]*(MKL_LONG)multiplier;
      // re += index_vec[n]*(MKL_LONG)pow(N_div, N_dimension - (n+1));
    }
  index_sca = re;
  return re;
}

MKL_LONG
UTILITY::
index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec,
	      const MKL_LONG N_dimension, const MKL_LONG *N_div)
{
  /*
    This inverse map applied to 
  */
  for(MKL_LONG n=0; n<N_dimension; n++)
    {
      double multiplier = 1.;
      for(MKL_LONG k=n+1; k<N_dimension; k++)
        {
          multiplier *= N_div[k];
        }
      index_vec[n] = ((MKL_LONG)(index_sca/multiplier))%N_div[n];
      // index_vec[n] = ((MKL_LONG)(index_sca/pow(N_div, N_dimension - (n+1))))%N_div;
    }
  return index_sca;
}

MKL_LONG
CLIST::
allocate_index_neighbor_cell_list()
{
  if(CELL_LIST_BOOST)
    {
      MKL_LONG* index_vec_boost = new MKL_LONG [N_dimension];
      MKL_LONG* sf_vec_boost = new MKL_LONG [3];
      // for(MKL_LONG i=0; i<N_neighbor_cells; i++) // this was bug.
      for(MKL_LONG i=0; i<N_cells; i++) // this is right. 
        {
          // note that the number of cells is N_div^N_dimension while number of neighbor cells is N_sf^N_dimension
          // for instance, if the PBC is divided by 5 cells for each axis, N_cells = 5^3 = 125 while the N_neighbor_cells = 3^3 = 9
          get_neighbor_cell_list(i, NEIGHBOR_CELLS[i], index_vec_boost, sf_vec_boost);
        }
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

MKL_LONG
CLIST::
get_neighbor_cell_list(const MKL_LONG& index_sca,
		       MKL_LONG* index_neighbor_cells,
		       MKL_LONG* self_index_vec_boost,
		       MKL_LONG* sf_vec_boost)
{
  /*
    get_neighbor_cell_list will store the index of neighbor_list into the index_neighbor_cells based on the scalar values.
    However, to call this function during the simulation make big overhead to compute index sets.
    Hence, we need to generate the array for the list when it is initialized, then re-used the generated list during simulation.
    For non-equilibrium simulation such as simple shear flow, the re-usability will decrease which eventually adds some overhead to compute additional lists.
  */
  CLIST::index_sca2vec(index_sca, self_index_vec_boost);
  MKL_LONG *N_sf = new MKL_LONG [N_dimension];
  MKL_LONG N_blocks = 1.;
  for(MKL_LONG n=0; n<N_dimension; n++)
    {
      N_sf[n] = 3; // {-1, 0, +1} for each direction
      N_blocks *= N_sf[n];
    }
  // MKL_LONG N_sf[] = {3, 3, 3}; // {-1, 0, +1} for each direction
  
  for(MKL_LONG nsf=0; nsf<N_blocks; nsf++)
    // nsf have value from 0 to number of all possible neighbor cells
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
              sf_vec_boost[n] += N_div[n];
            }
          else if(sf_vec_boost[n] >= N_div[n])
            {
              BEYOND_BOX[index_sca][nsf][n] = +1; // when the neighbor cell is beyond PBC boundary in right side
              sf_vec_boost[n] -= N_div[n];
            }
          
        }
      // it indicate that the sf_vec_boost have index vector for one of neighbor cell
      // while index_neighbor_cells[nsf] should have the scalar values (index_vec2sca: R^3 -> R)
      // of the given index vector
      UTILITY::index_vec2sca(sf_vec_boost, index_neighbor_cells[nsf], N_dimension, N_div);
    }
  return 0;
}

// MKL_LONG
// CLIST::
// get_neighbor_cell_list_dynamic_offset(const MKL_LONG& index_sca,
// 				      MKL_LONG* index_neighbor_cells,
// 				      MKL_LONG* self_index_vec_boost,
// 				      MKL_LONG* sf_vec_boost)
// {
//   CLIST::index_sca2vec(index_sca, self_index_vec_boost);
//   MKL_LONG N_sf = 3;

//   for(MKL_LONG nsf=0; nsf< pow(N_sf, N_dimension) + N_offset; nsf++)
//     {

//     }

// }



MKL_LONG
CLIST::
index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca)
{
  return UTILITY::index_vec2sca(index_vec, index_sca, N_dimension, N_div);
}

MKL_LONG
CLIST::
index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec)
{
  return UTILITY::index_sca2vec(index_sca, index_vec, N_dimension, N_div);
}


