
#ifndef POTENTIAL_H
#define POTENTIAL_H

#include <iostream>
#include <math.h>
#include <mkl.h>
/* #include "read_file_condition.h" */
#include "file_IO.h"
#include "matrix.h"
#include <string>
#define BIT 64

class POTENTIAL_SET
{
 public:
  
  // member variables
  double *force_variables;
  double (*f_connector)(double distance, double* given_variables);
  double (*e_connector)(double distance, double* given_variables);
  double (*PDF_connector)(double distance, double* given_variables);
  double (*f_repulsion)(double distance, double* given_variables);
  double (*e_repulsion)(double distance, double* given_variables);
  double (*scale_random)(MATRIX& basic_random_var_unity, double* given_variables);
  double (*w_function)(double distance, double tension, double* given_variables);
  double (*transition)(double distance, double tension, double* given_varialbes);

  POTENTIAL_SET()
    {
    }
  virtual ~POTENTIAL_SET()
    {
      if(force_variables)
        mkl_free(force_variables);
    }
};

namespace FORCE
{
  namespace GAUSSIAN
  {
    double spring_force(double distance, double N_dimension);
    double spring_potential(double distance, double N_dimension);
    double Boltzmann_distribution(double distance, double N_dimension);
  }
  namespace MODIFIED_GAUSSIAN
  {
    double spring_force(double distance, double N_dimension, double scale_factor);
    double spring_potential(double distance, double N_dimension, double scale_factor);
    double Boltzmann_distribution(double distance, double N_dimension, double scale_factor);
    double cutoff_Boltzmann_distribution(double distance, double N_dimension, double scale_factor, double cutoff_radius); 
  }
  namespace FENE
  {
    double non_Gaussian_factor(double distance, double N_dimension, double ratio_RM_R0);
    double spring_force(double distance, double N_dimension, double ratio_RM_R0);
    double spring_potential(double distance, double N_dimension, double ratio_RM_R0);
    double Boltzmann_distribution(double distance, double N_dimension, double ratio_RM_R0);
  }
  namespace DEFAULT
  {
    double EMPTY_force_contribution(double distance, double *given_variables);
    double basic_random(MATRIX& given_basic_random, double* given_variables);
    double time_scaling_random(MATRIX& given_basic_random, double scale_factor);
    MKL_LONG EMPTY_force_set(POTENTIAL_SET& given_POT, COND& given_condition);

  }
  namespace NAPLE
  {
    double excluded_volume_force(double distance, double effective_distance);
    double excluded_volume_potential(double distance, double effective_distance);
    
    namespace SIMPLE_REPULSION
    {
      double MAP_excluded_volume_force(double distance, double* given_variables);
      double MAP_excluded_volume_potential(double distance, double* given_variables);
      double MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables);
      MKL_LONG MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond);
    }
    namespace MC_ASSOCIATION
    {
      double MAP_Gaussian_spring_force(double distance, double* given_variables);
      double MAP_Gaussian_spring_potential(double distance, double* given_variables);
      double MAP_Gaussian_Boltzmann(double distance, double* given_variables);

      double MAP_modified_Gaussian_spring_force(double distance, double* given_variables);
      double MAP_modified_Gaussian_spring_potential(double distance, double* given_variables);
      double MAP_modified_Gaussian_Boltzmann(double distance, double* given_variables);
      double MAP_cutoff_modified_Gaussian_Boltzmann(double distance, double* given_variables);
      
      double MAP_FENE_spring_force(double distance, double* given_variables);
      double MAP_FENE_spring_potential(double distance, double* given_variables);
      double MAP_FENE_Boltzmann(double distance, double* given_variables);
      double MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables);
      MKL_LONG MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond);
    }

  }
}  

namespace KINETICS
{
  namespace UNIFORM
  {
    double detachment_weight(double distance, double tension, double* given_variables);
    double transition_probability(double distance, double tension, double* given_variables);
  }
  namespace WEIGHTED
  {
    double detachment_weight(double distance, double tension, double* given_variables);
    double transition_probability(double distance, double tension, double* given_variables);
  }

  namespace METROPOLIS
  {
    double detachment_weight(double distance, double tension, double* given_variables);
    double transition_probability(double distance, double tension, double* given_varialbes);
  }

  double dissociation_probability(double distance, double tension, double* given_variables);
  /*   double detachment_weight(double distance, double tension, double* given_variables); */
  /*   double transition_probability(double distance, double tension, double* given_varialbes); */
  /* } */

}



#endif
