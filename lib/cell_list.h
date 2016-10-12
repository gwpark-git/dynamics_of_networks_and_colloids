
#ifndef CELL_LIST_H
#define CELL_LIST_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <math.h>
#include "matrix.h" 
#include "trajectory.h" // the converting interface will be implemented
/* #include "read_file_condition.h" */
#include "file_IO.h"

class CLIST
{
 private:
 public:

  bool CELL_LIST_BOOST;
  MKL_LONG N_dimension; // spatial dimension
  double box_dimension; // box dimension for one axis
  bool INITIALIZATION;
  /*
    Based on the given critical radius, rc, the real length for cell is computed based on proper divisor of whole box. For instance, if the given box length is 10 and the critical length scale is 0.9, then the length for cell becomes 1.0, which is the slight larger than the critical length scale.
  */
  double cut_off_radius; // cut_off_radius
  double cell_length; // length for cell
  MKL_LONG N_div; // number of cells per axis
  MKL_LONG N_cells; // number of cells
  MKL_LONG N_neighbor_cells; // including itself
  MKL_LONG MAX_IN_CELL; // it state how many particles are allowed for one cell. Set as all the number of particles as default.
  
  MKL_LONG **CELL;
  MKL_LONG *cell_index;
  MKL_LONG **NEIGHBOR_CELLS;


  
  MKL_LONG ***BEYOND_BOX;
  MKL_LONG *TOKEN;
  MKL_LONG Np;

  // following are related with mechanical perturbations
  bool SIMPLE_SHEAR;
  bool STEP_SHEAR;
  MKL_LONG shear_axis;
  MKL_LONG shear_grad_axis;
  double map_to_central_box_image;
  
  // when simple shear is implemented
  // MKL_LONG N_offset;
  // MKL_LONG ***NEIGHBOR_CELLS_OFFSET; // the first dimensionality suggest that up and bottom check (right and bottom in general description)
  

  
  // mapping function
  MKL_LONG
    index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca);
  MKL_LONG
    index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec);
  MKL_LONG
    get_neighbor_cell_list(const MKL_LONG& index_sca,
			   MKL_LONG* index_neighbor_cells,
			   MKL_LONG* self_index_vec_boost,
			   MKL_LONG* sf_vec_boost);
  MKL_LONG
    allocate_index_neighbor_cell_list();
  MKL_LONG
    identify_cell_from_given_position(TRAJECTORY& TRAJ, MKL_LONG index_t_now,
				      MKL_LONG index_particle, MKL_LONG *index_vec_boost);
  MKL_LONG
    allocate_cells_from_positions(TRAJECTORY& TRAJ, MKL_LONG index_t_now,
				  MKL_LONG *index_vec_boost);
  MKL_LONG
    identify_cell_from_given_position(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t_now,
				      MKL_LONG index_particle, MKL_LONG *index_vec_boost);
  MKL_LONG
    allocate_cells_from_positions(TRAJECTORY_HDF5& TRAJ, MKL_LONG index_t_now,
				  MKL_LONG *index_vec_boost);

    
  // operator overloading
  MKL_LONG&
    operator()
    (MKL_LONG i, MKL_LONG j)
    {
      return CELL[i][j];
    }
  
  // constructor
  CLIST(){}
  CLIST(COND& given_condition);

  // destructor
  virtual
    ~CLIST()
    {
      if(INITIALIZATION)
        {
          delete[] TOKEN;
          for(MKL_LONG i=0; i<N_cells; i++)
            {
              delete[] CELL[i];
              delete[] NEIGHBOR_CELLS[i];
              for(MKL_LONG j=0; j<N_neighbor_cells; j++)
                delete[] BEYOND_BOX[i][j];
              delete[] BEYOND_BOX[i];
            }
          delete[] CELL;
          delete[] NEIGHBOR_CELLS;
          delete[] BEYOND_BOX;
          delete[] cell_index;

          // if(SIMPLE_SHEAR)
          //   {
          //     for(MKL_LONG i=0; i<2; i++)
          //       {
          //         for(MKL_LONG j=0; j<N_cells + 1; j++)
          //           {
          //             for(MKL_LONG k=0; k<N_neighbor_cells + N_offset; k++)
          //               {
          //                 delete[] NEIGHBOR_CELLS_OFFSET[i][j][k];
          //               }
          //             delete[] NEIGHBOR_CELLS_OFFSET[i][j];
          //           }
          //         delete[] NEIGHBOR_CELLS_OFFSET[i];
          //       }
          //   }
        }
    }
};

namespace UTILITY
{
  MKL_LONG
    index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca,
		  const MKL_LONG N_dimension, const MKL_LONG N_div);
  MKL_LONG
    index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec,
		  const MKL_LONG N_dimension, const MKL_LONG N_div);
}

#endif
