
#ifndef CELL_LIST_H
#define CELL_LIST_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <math.h>
#include "matrix.h" 
#include "trajectory.h" // the converting interface will be implemented
#inlcude "read_file_condition.h"

class CLIST
{
 private:
 public:

  MKL_LONG N_dimension; // spatial dimension
  double box_dimension; // box dimension for one axis
  /*
    Based on the given critical radius, rc, the real length for cell is computed based on proper divisor of whole box. For instance, if the given box length is 10 and the critical length scale is 0.9, then the length for cell becomes 1.0, which is the slight larger than the critical length scale.
  */
  double cut_off_radius; // cut_off_radius
  double cell_length; // length for cell
  MKL_LONG N_div; // number of cells per axis
  MKL_LONG N_cells; // number of cells
  MKL_LONG N_neighbor_cells; // including itself

  MATRIX *CELL;
  MKL_LONG **NEIGHBOR_CELLS;
  MKL_LONG *TOKEN;
  MKL_LONG Np;
  // mapping function
  MKL_LONG index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca);
  MKL_LONG index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec);
  MKL_LONG get_neighbor_cell_list(const MKL_LONG& index_sca, MKL_LONG* index_neighbor_cells, MKL_LONG* self_index_vec_boost, MKL_LONG* sf_vec_boost);
  MKL_LONG allocate_index_neighbor_cell_list();
  MKL_LONG identify_cell_from_given_position(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG index_particle, MKL_LONG *index_vec_boost);
  MKL_LONG allocate_cells_from_positions(TRAJECTORY& TRAJ, MKL_LONG index_t_now, MKL_LONG *index_vec_boost);
  // constructor
  CLIST(){}
  /* CLIST(MKL_LONG given_N_dimension, double given_cut_off_radius, double given_box_dimension, MKL_LONG max_particle_in_cell) */
  CLIST(COND& given_condition)
    {
      /* cut_off_radius = given_cut_off_radius; */
      /* N_dimension = given_N_dimension; */
      /* box_dimension = given_box_dimension; */
      cut_off_radius = atof(given_condition("cutoff_connection").c_str());
      N_dimension = atoi(given_condition("N_dimension").c_str());
      box_dimension = atof(given_condition("box_dimension").c_str());
      Np = atoi(given_condition("Np").c_str());
      N_div = (MKL_LONG) box_dimension/cut_off_radius; // floor refine the given divisor
      cell_length = (double) box_dimension/(double)N_div; // real length scale of cells in each axis
      N_cells = pow(N_div, N_dimension); // N_cells = N_div^N_dimension
      N_neighbor_cells = pow(3, N_dimension); // 3 means number for shift factor {-1, 0, +1}
      MAX_IN_CELL = max_particle_in_cell;
      TOKEN = (MKL_LONG*)mkl_malloc(N_cells*sizeof(MKL_LONG), BIT);
      CELL = (MATRIX*)mkl_malloc(N_cells*sizeof(MATRIX), BIT);
      NEIGHBOR_CELLS = (MKL_LONG**)mkl_malloc(N_cells*sizeof(MKL_LONG*), BIT);
      for(MKL_LONG i=0; i<N_cells; i++)
        {
          CELL.initial(MAX_IN_CELL, 1, 0.);
          TOKEN[i] = -1; // -1 is default value when there is no index (note that the particle index is started with 0, which is the reason to avoid using 0 identifier
          NEIGHBOR_CELLS[i] = (MKL_LONG*)mkl_malloc(N_neighbor_cells*sizeof(MKL_LONG), BIT);
          /* get_neighbor_cell_list(i,  */
        }
      allocate_index_neighbor_cell_list();
    }
  virtual ~CLIST()
    {
      mkl_free(CELL);
      mkl_free(TOKEN);
    }
}

namespace UTILITY
{
  MKL_LONG index_vec2sca(const MKL_LONG* index_vec, MKL_LONG& index_sca, const MKL_LONG N_dimension, const MKL_LONG N_div);
  MKL_LONG index_sca2vec(const MKL_LONG& index_sca, MKL_LONG* index_vec, const MKL_LONG N_dimension, const MKL_LONG N_div);
}

#endif
