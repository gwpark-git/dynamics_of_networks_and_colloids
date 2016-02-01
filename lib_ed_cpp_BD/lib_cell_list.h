
#ifndef LIB_CELL_LIST_H
#define LIB_CELL_LIST_H

#define TRUE 1
#define FALSE 0
#define BIT 64

#include <iostream>
#include <math.h>
#include "matrix_ed.h" 
#include "lib_traj.h" // the converting interface will be implemented


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

  MATRIX *CELL;
  MKL_LONG *TOKEN;
  
  // mapping function
  /* MKL_LONG operator(MKL_LONG* index_vec); // index_vec has Nd components */
  MKL_LONG cell_index(MKL_LONG* index_vec);
  // constructor
  CLIST(){}
  CLIST(MKL_LONG given_N_dimension, double given_cut_off_radius, double given_box_dimension, MKL_LONG max_particle_in_cell)
    {
      cut_off_radius = given_cut_off_radius;
      N_dimension = given_N_dimension;
      box_dimension = given_box_dimension;
      N_div = (MKL_LONG) box_dimension/cut_off_radius; // floor refine the given divisor
      cell_length = (double) box_dimension/(double)N_div; // real length scale of cells in each axis
      N_cells = pow(N_div, N_dimension); // N_cells = N_div^N_dimension
      MAX_IN_CELL = max_particle_in_cell;
      TOKEN = (MKL_LONG*)mkl_malloc(N_cells*sizeof(MKL_LONG), BIT);
      CELL = (MATRIX*)mkl_malloc(N_cells*sizeof(MATRIX), BIT);
      for(MKL_LONG i=0; i<N_cells; i++)
	{
	  CELL.initial(MAX_IN_CELL, 1, 0.);
	  TOKEN[i] = -1; // -1 is default value when there is no index (note that the particle index is started with 0, which is the reason to avoid using 0 identifier
	}
    }
  virtual ~CLIST()
    {
      mkl_free(CELL);
      mkl_free(TOKEN);
    }
}


#endif
