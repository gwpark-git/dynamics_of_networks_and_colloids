
#include "potential.h"

MKL_LONG FORCE::DEFAULT::EMPTY_force_set(POTENTIAL_SET& given_POT, COND& given_condition)
{
  given_POT.f_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.f_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.scale_random = FORCE::DEFAULT::basic_random;
  return 0;
}

MKL_LONG FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond)
{
  given_POT.force_variables = (double*) mkl_calloc(3, sizeof(double), BIT);
  given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);

  given_POT.f_repulsion = MAP_excluded_volume_force;
  given_POT.e_repulsion = MAP_excluded_volume_potential;

  given_POT.f_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_connector = FORCE::DEFAULT::EMPTY_force_contribution;

  given_POT.scale_random = MAP_time_scaling_random;

  return 0;
}


MKL_LONG FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_set(POTENTIAL_SET& given_POT, COND& given_cond)
{
  given_POT.force_variables = (double*) mkl_calloc(8, sizeof(double), BIT);
  given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);
  given_POT.force_variables[3] = atol(given_cond("N_dimension").c_str());
  given_POT.force_variables[4] = atof(given_cond("l_cap").c_str());

  given_POT.f_repulsion = FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force;
  given_POT.e_repulsion = FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential;

  if(given_cond("connector")=="FENE")
    {
      given_POT.force_variables[5] = atof(given_cond("ratio_RM_R0").c_str());
      given_POT.f_connector = MAP_FENE_spring_force;
      given_POT.e_connector = MAP_FENE_spring_potential;
      given_POT.PDF_connector = MAP_FENE_Boltzmann;
    }
  else if (given_cond("connector") == "Gaussian")
    {
      given_POT.f_connector = MAP_Gaussian_spring_force;
      given_POT.e_connector = MAP_Gaussian_spring_potential;
      given_POT.PDF_connector = MAP_Gaussian_Boltzmann;
    }
  else if (given_cond("connector") == "Modified_Gaussian")
    {
      given_POT.force_variables[5] = atof(given_cond("scale_factor_chain").c_str());
      given_POT.f_connector = MAP_modified_Gaussian_spring_force;
      given_POT.e_connector = MAP_modified_Gaussian_spring_potential;
      double cutoff_connection = atof(given_cond("cutoff_connection").c_str());
      if (cutoff_connection > 0.0 && cutoff_connection < atof(given_cond("box_dimension").c_str()))
        {
          given_POT.force_variables[7] = cutoff_connection;
          given_POT.PDF_connector = MAP_cutoff_modified_Gaussian_Boltzmann;
        }
      else
        {
          given_POT.PDF_connector = MAP_modified_Gaussian_Boltzmann;
        }
      // given_POT.PDF_connector = MAP_modified_Gaussian_Boltzmann;
    }
  else
    {
      printf("ERR: no avaliable connector inp.\n");
      return -1;
    }

  if(given_cond("transition_probability")=="UNIFORM")
    {
      given_POT.transition = KINETICS::UNIFORM::transition_probability;
    }
  else if (given_cond("transition_probability")=="METROPOLIS")
    {
      given_POT.transition = KINETICS::METROPOLIS::transition_probability;
    }
  else if (given_cond("transition_probability")=="WEIGHTED")
    {
      given_POT.transition = KINETICS::WEIGHTED::transition_probability;
    }
  else if (given_cond("transition_probability")=="DISSOCIATION")
    {
      given_POT.transition = KINETICS::dissociation_probability;
      // we have to check Dt for topological update
      // note that it is dt*N_steps_block
      given_POT.force_variables[6] = atof(given_cond("dt").c_str())*atof(given_cond("N_steps_block").c_str())/atof(given_cond("Rt").c_str());
    }
  else if (given_cond("transition_probability")=="FIRST_ORDER")
    {
      given_POT.transition = KINETICS::FIRST_ORDER::dissociation_probability;
      given_POT.force_variables[6] = atof(given_cond("dt").c_str())*atof(given_cond("N_steps_block").c_str())/atof(given_cond("Rt").c_str());

    }
  else
    {
      // this is default setting
      given_POT.transition = KINETICS::UNIFORM::transition_probability;
    }
  
  if(given_cond("chain_selection")=="METROPOLIS")
    {
      given_POT.force_variables[6] = atof(given_cond("energy_barrier").c_str());
      given_POT.w_function = KINETICS::METROPOLIS::detachment_weight;
    }
  else if (given_cond("chain_selection")=="WEIGHTED")
    {
      given_POT.w_function = KINETICS::WEIGHTED::detachment_weight;
    }
  else if (given_cond("chain_selection")=="UNIFORM")
    {
      given_POT.w_function = KINETICS::UNIFORM::detachment_weight;
    }
  else
    {
      printf("ERR: no avaliable kinetics inp.\n");
      return -1;
    }
  
  given_POT.scale_random = MAP_time_scaling_random;
  return 0;
}

// double FORCE::DEFAULT::EMPTY_force_contribution(double distance, double* given_varialbes)
// {
//   return 0;
// }

// double FORCE::DEFAULT::basic_random(MATRIX& given_basic_random, double* given_variables)
// {
//   return 0;
// }


// double FORCE::DEFAULT::time_scaling_random(MATRIX& given_basic_random, double scale_factor)
// {
//   matrix_mul(given_basic_random, scale_factor); // this is varied
//   return scale_factor;
// }

// double FORCE::NAPLE::SIMPLE_REPULSION::MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables)
// {
//   return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_time_scaling_random(MATRIX& given_basic_random, double* given_variables)
// {
//   return FORCE::DEFAULT::time_scaling_random(given_basic_random, given_variables[2]);
// }

// double FORCE::NAPLE::excluded_volume_force(double distance, double effective_distance)
// {
//   if (distance < effective_distance)
//     {
//       return -(1. - pow(distance, 2.0));
//     }
//   return 0.;
// }

// double FORCE::NAPLE::excluded_volume_potential(double distance, double effective_distance)
// {
//   if (distance < effective_distance)
//       return (1./3.)*pow(1. - distance, 2.)*(2. + distance);
//   return 0.;
// }

// double FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force(double distance, double* given_variables)
// {
//   return FORCE::NAPLE::excluded_volume_force(distance, given_variables[1]);
// }

// double FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential(double distance, double* given_variables)
// {
//   return FORCE::NAPLE::excluded_volume_potential(distance, given_variables[1]);
// }



// double FORCE::GAUSSIAN::spring_force(double distance, double N_dimension)
// {
//   return N_dimension*distance;
// }

// double FORCE::GAUSSIAN::spring_potential(double distance, double N_dimension)
// {
//   return (0.5)*N_dimension*distance*distance;
// }

// double FORCE::GAUSSIAN::Boltzmann_distribution(double distance, double N_dimension)
// {
//   return exp(-spring_potential(distance, N_dimension));
// }

// double FORCE::MODIFIED_GAUSSIAN::spring_force(double distance, double N_dimension, double scale_factor)
// {
//   return FORCE::GAUSSIAN::spring_force(scale_factor*scale_factor*distance, N_dimension);
// }

// double FORCE::MODIFIED_GAUSSIAN::spring_potential(double distance, double N_dimension, double scale_factor)
// {
//   return FORCE::GAUSSIAN::spring_potential(scale_factor*distance, N_dimension);
// }

// double FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(double distance, double N_dimension, double scale_factor)
// {
//   return FORCE::GAUSSIAN::Boltzmann_distribution(scale_factor*distance, N_dimension);
// }

// double FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(double distance, double N_dimension, double scale_factor, double cutoff_radius)
// {
//   // note that this cut-off do not have any benefit at this moment
//   // after implementaion for cell-list, the cut-off scheme becomes efficience
//   if (distance < cutoff_radius)
//     return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, N_dimension, scale_factor);
//   return 0.;
// }

// double FORCE::FENE::non_Gaussian_factor(double distance, double N_dimension, double ratio_RM_R0)
// {
//   return 1.0/(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
// }

// double FORCE::FENE::spring_force(double distance, double N_dimension, double ratio_RM_R0)
// {
//   return non_Gaussian_factor(distance, N_dimension, ratio_RM_R0)*FORCE::GAUSSIAN::spring_force(distance, N_dimension);
// }

// double FORCE::FENE::spring_potential(double distance, double N_dimension, double ratio_RM_R0)
// {
//   return (-(double)N_dimension/2.0)*pow(ratio_RM_R0, 2.0)*log(1.0 - pow(distance, 2.0)/pow(ratio_RM_R0, 2.0));
// }

// double FORCE::FENE::Boltzmann_distribution(double distance, double N_dimension, double ratio_RM_R0)
// {
//   // Note that even without cut-off range based on the ratio_RM_R0, the code is properly working.
//   // However, when we applied this scheme, the sorting procedure is more stable and normalized in naturally.
//   if (distance < ratio_RM_R0)
//     return exp(-FORCE::FENE::spring_potential(distance, N_dimension, ratio_RM_R0));
//   return 0.;
// }


// double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_force(double distance, double* given_variables)
// {
//   return FORCE::GAUSSIAN::spring_force(distance, given_variables[3]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_spring_potential(double distance, double* given_variables)
// {
//   return FORCE::GAUSSIAN::spring_potential(distance, given_variables[3]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_Gaussian_Boltzmann(double distance, double* given_variables)
// {
//   return FORCE::GAUSSIAN::Boltzmann_distribution(distance, given_variables[3]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_force(double distance, double* given_variables)
// {
//   return FORCE::MODIFIED_GAUSSIAN::spring_force(distance, given_variables[3], given_variables[5]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_spring_potential(double distance, double* given_variables)
// {
//   return FORCE::MODIFIED_GAUSSIAN::spring_potential(distance, given_variables[3], given_variables[5]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_modified_Gaussian_Boltzmann(double distance, double* given_variables)
// {
//   return FORCE::MODIFIED_GAUSSIAN::Boltzmann_distribution(distance, given_variables[3], given_variables[5]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_cutoff_modified_Gaussian_Boltzmann(double distance, double* given_variables)
// {
//   return FORCE::MODIFIED_GAUSSIAN::cutoff_Boltzmann_distribution(distance, given_variables[3], given_variables[5], given_variables[7]);
// }


// double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_force(double distance, double* given_variables)
// {
//   return FORCE::FENE::spring_force(distance, given_variables[3], given_variables[5]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_spring_potential(double distance, double* given_variables)
// {
//   return FORCE::FENE::spring_potential(distance, given_variables[3], given_variables[5]);
// }

// double FORCE::NAPLE::MC_ASSOCIATION::MAP_FENE_Boltzmann(double distance, double* given_variables)
// {
//   return FORCE::FENE::Boltzmann_distribution(distance, given_variables[3], given_variables[5]);
// }




// double KINETICS::WEIGHTED::detachment_weight(double distance, double tension, double* given_variables)
// {
//   // given_variables[3] == Nd
//   // given_variables[4] == l_cap
//   // in normalized scheme, the energy barrier is canceled out with normalization factor
//   // therefore, it is only affected by the tension exerted on the chain
//   return exp(tension*given_variables[4]);
// }


// double KINETICS::METROPOLIS::detachment_weight(double distance, double tension, double* given_variables)
// {
//   return 1.0;
// }

// double KINETICS::METROPOLIS::transition_probability(double distance, double tension, double* given_variables)
// {
//   double tpa = exp(tension*given_variables[4] - given_variables[6]);
//   if (tpa > 1.0)
//     return 1.0;
//   return tpa;
// }

// double KINETICS::WEIGHTED::transition_probability(double distance, double tension, double* given_variables)
// {
//   double tpa = exp(tension*given_variables[4]);
//   if (tpa > 1.0)
//     return 1.0;
//   return tpa;
// }

// double KINETICS::UNIFORM::detachment_weight(double distance, double tension, double *given_variables)
// {
//   return 1.0;
// }

// double KINETICS::UNIFORM::transition_probability(double distance, double tension, double* given_varialbes)
// {
//   return 1.0;
// }

// double KINETICS::dissociation_probability(double distance, double tension, double* given_variables)
// {
//   double tpa = given_variables[6]*exp(tension*given_variables[4]);
//   if (tpa > 1.0)
//     return 1.0;
//   return tpa;
// }

// double KINETICS::FIRST_ORDER::dissociation_probability(double distance, double tension, double* given_variables)
// {
//   double Dt = given_variables[6];
//   double lcap = given_variables[4];
//   double beta = exp(tension*lcap);
//   double tpa = 1.0 - exp(-beta*Dt);
//   // double tpa = given_variables[6]*exp(tension*given_variables[4]);
//   // if (tpa > 1.0)
//   //   return 1.0;
//   return tpa;
// }
