
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
        /* mkl_free(force_variables); */
        delete[] force_variables;
    }
};


namespace FORCE
{
  namespace BROWNIAN
  {
    MKL_LONG MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond);
    double no_time_scaling_random(MATRIX& given_basic_random, double scale_factor);
    double MAP_no_time_scaling_random(MATRIX& given_basic_random, double *given_variables);
  }
  
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
  namespace FIRST_ORDER
  {
    double dissociation_probability(double distance, double tension, double* given_variables);
  }
  /*   double detachment_weight(double distance, double tension, double* given_variables); */
  /*   double transition_probability(double distance, double tension, double* given_varialbes); */
  /* } */

}


// inline functions
// note that all the inline function 'must be' defined on above for readability.

inline double KINETICS::WEIGHTED::detachment_weight(double distance, double tension, double* given_variables)
{
  // given_variables[3] == Nd
  // given_variables[4] == l_cap
  // in normalized scheme, the energy barrier is canceled out with normalization factor
  // therefore, it is only affected by the tension exerted on the chain
  return exp(tension*given_variables[4]);
}


inline double KINETICS::METROPOLIS::detachment_weight(double distance, double tension, double* given_variables)
{
  return 1.0;
}

inline double KINETICS::METROPOLIS::transition_probability(double distance, double tension, double* given_variables)
{
  double tpa = exp(tension*given_variables[4] - given_variables[6]);
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double KINETICS::WEIGHTED::transition_probability(double distance, double tension, double* given_variables)
{
  double tpa = exp(tension*given_variables[4]);
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double KINETICS::UNIFORM::detachment_weight(double distance, double tension, double *given_variables)
{
  return 1.0;
}

inline double KINETICS::UNIFORM::transition_probability(double distance, double tension, double* given_varialbes)
{
  return 1.0;
}

inline double KINETICS::dissociation_probability(double distance, double tension, double* given_variables)
{
  double tpa = given_variables[6]*exp(tension*given_variables[4]);
  if (tpa > 1.0)
    return 1.0;
  return tpa;
}

inline double KINETICS::FIRST_ORDER::dissociation_probability(double distance, double tension, double* given_variables)
{
  double Dt = given_variables[6];
  double lcap = given_variables[4];
  double beta = exp(tension*lcap);
  double tpa = 1.0 - exp(-beta*Dt);
  // double tpa = given_variables[6]*exp(tension*given_variables[4]);
  // if (tpa > 1.0)
  //   return 1.0;
  return tpa;
}


inline double FORCE::DEFAULT::EMPTY_force_contribution(double distance, double* given_varialbes)
{
  return 0;
}

inline double FORCE::DEFAULT::basic_random(MATRIX& given_basic_random, double* given_variables)
{
  return 0;
}


inline double FORCE::DEFAULT::time_scaling_random(MATRIX& given_basic_random, double scale_factor)
{
  matrix_mul(given_basic_random, scale_factor); // this is varied due to the given time scale
  return scale_factor;
}

inline double FORCE::NAPLE::SIMPLE_REPULSION::MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables)
{
  return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables)
{
  return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
}

inline double FORCE::NAPLE::excluded_volume_force(double distance, double effective_distance)
{
  if (distance < effective_distance)
    {
      return -(1. - pow(distance, 2.0));
    }
  return 0.;
}

inline double FORCE::NAPLE::excluded_volume_potential(double distance, double effective_distance)
{
  if (distance < effective_distance)
      return (1./3.)*pow(1. - distance, 2.)*(2. + distance);
  return 0.;
}

inline double FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force(double distance, double* given_variables)
{
  return FORCE::NAPLE::excluded_volume_force(distance, given_variables[1]);
}

inline double FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential(double distance, double* given_variables)
{
  return FORCE::NAPLE::excluded_volume_potential(distance, given_variables[1]);
}


inline double FORCE::GAUSSIAN::spring_force(double distance, double N_dimension)
{
  return N_dimension*distance;
}

inline double FORCE::GAUSSIAN::spring_potential(double distance, double N_dimension)
{
  return (0.5)*N_dimension*distance*distance;
}

inline double FORCE::GAUSSIAN::Boltzmann_distribution(double distance, double N_dimension)
{
  return exp(-spring_potential(distance, N_dimension));
}

inline double FORCE::MODIFIED_GAUSSIAN::spring_force(double distance, double N_dimension, double scale_factor)
{
  return FORCE::GAUSSIAN::spring_force(scale_factor*scale_factor*distance, N_dimension);
}

inline double FORCE::MODIFIED_GAUSSIAN::spring_potential(double distance, double N_dimension, double scale_factor)
{
  return FORCE::GAUSSIAN::spring_potential(scale_factor*distance, N_dimension);
}

inline double FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(double distance, double N_dimension, double scale_factor)
{
  return FORCE::GAUSSIAN::Boltzmann_distribution(scale_factor*distance, N_dimension);
}

inline double FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(double distance, double N_dimension, double scale_factor, double cutoff_radius)
{
  // note that this cut-off do not have any benefit at this moment
  // after implementaion for cell-list, the cut-off scheme becomes efficience
  if (distance < cutoff_radius)
    return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, N_dimension, scale_factor);
  return 0.;
}

inline double FORCE::FENE::non_Gaussian_factor(double distance, double N_dimension, double ratio_RM_R0)
{
  return 1.0/(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
}

inline double FORCE::FENE::spring_force(double distance, double N_dimension, double ratio_RM_R0)
{
  return non_Gaussian_factor(distance, N_dimension, ratio_RM_R0)*FORCE::GAUSSIAN::spring_force(distance, N_dimension);
}

inline double FORCE::FENE::spring_potential(double distance, double N_dimension, double ratio_RM_R0)
{
  return (-(double)N_dimension/2.0)*pow(ratio_RM_R0, 2.0)*log(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
}

inline double FORCE::FENE::Boltzmann_distribution(double distance, double N_dimension, double ratio_RM_R0)
{
  // Note that even without cut-off range based on the ratio_RM_R0, the code is properly working.
  // However, when we applied this scheme, the sorting procedure is more stable and normalized in naturally.
  if (distance < ratio_RM_R0)
    return exp(-FORCE::FENE::spring_potential(distance, N_dimension, ratio_RM_R0));
  return 0.;
}


inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force(double distance, double* given_variables)
{
  return FORCE::GAUSSIAN::spring_force(distance, given_variables[3]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_potential(double distance, double* given_variables)
{
  return FORCE::GAUSSIAN::spring_potential(distance, given_variables[3]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_Boltzmann(double distance, double* given_variables)
{
  return FORCE::GAUSSIAN::Boltzmann_distribution(distance, given_variables[3]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_force(double distance, double* given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::spring_force(distance, given_variables[3], given_variables[5]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential(double distance, double* given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::spring_potential(distance, given_variables[3], given_variables[5]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_Boltzmann(double distance, double* given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, given_variables[3], given_variables[5]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_cutoff_modified_Gaussian_Boltzmann(double distance, double* given_variables)
{
  return FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(distance, given_variables[3], given_variables[5], given_variables[7]);
}


inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_force(double distance, double* given_variables)
{
  return FORCE::FENE::spring_force(distance, given_variables[3], given_variables[5]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_potential(double distance, double* given_variables)
{
  return FORCE::FENE::spring_potential(distance, given_variables[3], given_variables[5]);
}

inline double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_Boltzmann(double distance, double* given_variables)
{
  return FORCE::FENE::Boltzmann_distribution(distance, given_variables[3], given_variables[5]);
}


#endif
