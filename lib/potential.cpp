
#include "potential.h"

MKL_LONG
FORCE::BROWNIAN::
MAP_potential_variable
(POTENTIAL_SET& given_POT, COND& given_cond)
{
  // given_POT.force_variables = NULL;
  given_POT.force_variables.N_dimension = atol(given_cond("N_dimension").c_str());
  return 0;
}

MKL_LONG
FORCE::DUMBBELL::
MAP_potential_variable
(POTENTIAL_SET& given_POT, COND& given_cond) 
{
  FORCE::BROWNIAN::MAP_potential_variable(given_POT, given_cond);  
 given_POT.force_variables.scale_factor_chain = atof(given_cond("scale_factor_chain").c_str());
 return 0;
}

MKL_LONG
FORCE::NAPLE::SIMPLE_REPULSION::
MAP_potential_variable
(POTENTIAL_SET& given_POT, COND& given_cond) 
{
  FORCE::BROWNIAN::MAP_potential_variable(given_POT, given_cond);
 
  given_POT.force_variables.repulsion_coefficient = atof(given_cond("repulsion_coefficient").c_str());
  given_POT.force_variables.inv_sqrt_repulsion_coefficient = 1./sqrt(given_POT.force_variables.repulsion_coefficient);  
  // given_POT.force_variables.C_rep_arr.initial(Np, 1, 1.);
  // note that initialization setting 1 means C_rep/C_rep=1 for all the particle.
  // this is due to the fact that the repulsive coefficient is invariance for repulsive Brownian dynamics (RBD)
  given_POT.force_variables.effective_distance = atof(given_cond("effective_distance").c_str());

  return 0;
}

MKL_LONG
FORCE::NAPLE::MC_ASSOCIATION::
MAP_potential_variable
(POTENTIAL_SET& given_POT, COND& given_cond) 
{
  FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_variable(given_POT, given_cond);
  given_POT.force_variables.scale_factor_chain = atof(given_cond("scale_factor_chain").c_str());
  given_POT.force_variables.l_cap = atof(given_cond("l_cap").c_str());
  return 0;
}


double
FORCE::BROWNIAN::
no_time_scaling_random
(MATRIX& given_basic_random, double scale_factor)
{
  return 1.; // indicate no changes on given_basic_random;
}

double
FORCE::BROWNIAN::
MAP_no_time_scaling_random
(MATRIX& given_basic_random, POTENTIAL_VARIABLE& given_variables)
{
  // the last argument *given_variables is related with POTs.force_variables. In the case for pure Brownian motion, however, the pointer for force_variables will be set with NULL. Therefore, it will be ignored, but the argument remains in order to use function pointer.
  return FORCE::BROWNIAN::no_time_scaling_random(given_basic_random, 1.); // anyhow, the unity will be ignored for performance
}


MKL_LONG
FORCE::BROWNIAN::
MAP_potential_set
(POTENTIAL_SET& given_POT, COND& given_cond)
{
  // given_POT.force_variables = NULL;
  FORCE::BROWNIAN::MAP_potential_variable(given_POT, given_cond);
  given_POT.f_repulsion = NULL;
  given_POT.f_repulsion_coupling = NULL;
  given_POT.e_repulsion = NULL;
  given_POT.f_connector = NULL;
  given_POT.e_connector = NULL;
  given_POT.scale_random = FORCE::BROWNIAN::MAP_no_time_scaling_random; // without any changes for time scaling since pure Brownian motion uses the basic time scale
  return 0;
}



MKL_LONG
FORCE::DUMBBELL::
MAP_potential_set
(POTENTIAL_SET& given_POT, COND& given_cond)
{
  /* 
     At this moment, the only Gaussian chain is implemented for dumbbell model.
     Note that Dumbbell model is not necessary any cut-off since there are only one permanent connection
  */
  // given_POT.force_variables = new double [3];
  // given_POT.force_variables[0] = atol(given_cond("N_dimension").c_str());
  // given_POT.force_variables[1] = atof(given_cond("scale_factor_chain").c_str());
  // given_POT.force_variables.N_dimension = atol(given_cond("N_dimension").c_str());
  // given_POT.force_variables.scale_factor_chain = atof(given_cond("scale_factor_chain").c_str());
  FORCE::DUMBBELL::MAP_potential_variable(given_POT, given_cond);
  given_POT.f_connector = FORCE::DUMBBELL::MAP_modified_Gaussian_spring_force;
  given_POT.e_connector = FORCE::DUMBBELL::MAP_modified_Gaussian_spring_potential;
  given_POT.scale_random = FORCE::BROWNIAN::MAP_no_time_scaling_random;
  
  return 0;
}

MKL_LONG
FORCE::DEFAULT::
EMPTY_force_set
(POTENTIAL_SET& given_POT, COND& given_condition)
{
  given_POT.f_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.f_repulsion_coupling = NULL;
  given_POT.e_repulsion = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.f_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.scale_random = FORCE::DEFAULT::basic_random;
  return 0;
}


MKL_LONG
FORCE::NAPLE::SIMPLE_REPULSION::
MAP_potential_set
(POTENTIAL_SET& given_POT, COND& given_cond)
{

  // given_POT.force_variables = new double [3];
  // given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  // given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  // given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);
  // given_POT.force_variables.repulsion_coefficient = atof(given_cond("repulsion_coefficient").c_str());
  // given_POT.force_variables.effective_distance = atof(given_cond("effective_distance").c_str());
  // given_POT.force_variables.inv_sqrt_repulsion_coefficient = 1./sqrt(given_POT.force_variables.repulsion_coefficient);
  FORCE::NAPLE::SIMPLE_REPULSION::MAP_potential_variable(given_POT, given_cond);
  
  given_POT.f_repulsion = MAP_excluded_volume_force;
  given_POT.e_repulsion = MAP_excluded_volume_potential;

  given_POT.f_repulsion_coupling = FORCE::NAPLE::SIMPLE_REPULSION::pre_averaged_repulsive_coefficient;
  
  given_POT.f_connector = FORCE::DEFAULT::EMPTY_force_contribution;
  given_POT.e_connector = FORCE::DEFAULT::EMPTY_force_contribution;

  given_POT.scale_random = MAP_time_scaling_random;

  return 0;
}


// MKL_LONG
// FORCE::NAPLE::MC_ASSOCIATION::
// MAP_potential_set_FENE(POTENTIAL_SET& given_POT)
// {
//   given_POT.force_variables.scale_factor_chain = atof(given_cond("scale_factor_chain").c_str());
//   given_POT.force_variables.ratio_RM_R0 = atof(given_cond("ratio_RM_R0").c_str());
  
// }


MKL_LONG
FORCE::NAPLE::MC_ASSOCIATION::
MAP_potential_set
(POTENTIAL_SET& given_POT, COND& given_cond)
{
  
  // given_POT.force_variables = new double [10];
  // given_POT.force_variables[0] = atof(given_cond("repulsion_coefficient").c_str());
  // given_POT.force_variables[1] = atof(given_cond("effective_distance").c_str());
  // given_POT.force_variables[2] = 1./sqrt(given_POT.force_variables[0]);
  // given_POT.force_variables[3] = atol(given_cond("N_dimension").c_str());
  // given_POT.force_variables[4] = atof(given_cond("l_cap").c_str());

  FORCE::NAPLE::MC_ASSOCIATION::MAP_potential_variable(given_POT, given_cond); 
  given_POT.f_repulsion = FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_force;
  given_POT.e_repulsion = FORCE::NAPLE::SIMPLE_REPULSION::MAP_excluded_volume_potential;

  if(given_cond("repulsion_coupling")=="GEOMETRIC_MEAN")
    {
      given_POT.f_repulsion_coupling = FORCE::NAPLE::SIMPLE_REPULSION::geometrical_mean_repulsive_coefficient;
    }
  else
    {
      /*
        default
       */
      given_POT.f_repulsion_coupling = FORCE::NAPLE::SIMPLE_REPULSION::pre_averaged_repulsive_coefficient;
    }
  
  // given_POT.f_repulsion_coupling = FORCE::NAPLE::SIMPLE_REPULSION::pre_averaged_repulsive_coefficient;
  
  if(given_cond("connector")=="FENE")
    {
      // given_POT.force_variables[5] = atof(given_cond("scale_factor_chain").c_str());
      // given_POT.force_variables[8] = atof(given_cond("ratio_RM_R0").c_str());
      given_POT.force_variables.scale_factor_chain = atof(given_cond("scale_factor_chain").c_str());
      given_POT.force_variables.ratio_RM_R0 = atof(given_cond("ratio_RM_R0").c_str());
      given_POT.f_connector = MAP_FENE_spring_force;
      given_POT.e_connector = MAP_FENE_spring_potential;
      given_POT.PDF_connector = MAP_FENE_Boltzmann;
      if(given_cond("association_probability") == "minimum_R0_Boltzmann")
        {
          given_POT.PDF_connector = MAP_minimum_R0_FENE_Boltzmann;
        }
    }
  else if (given_cond("connector") == "Gaussian")
    {
      given_POT.f_connector = MAP_Gaussian_spring_force;
      given_POT.e_connector = MAP_Gaussian_spring_potential;
      given_POT.PDF_connector = MAP_Gaussian_Boltzmann;
      if(given_cond("association_probability") == "minimum_R0_Boltzmann")
	{
	  given_POT.PDF_connector = MAP_minimum_R0_Gaussian_Boltzmann;
	}
      
    }
  else if (given_cond("connector") == "Modified_Gaussian")
    {
      // given_POT.force_variables[5] = atof(given_cond("scale_factor_chain").c_str());
      given_POT.force_variables.scale_factor_chain = atof(given_cond("scale_factor_chain").c_str());
      given_POT.f_connector = MAP_modified_Gaussian_spring_force;
      given_POT.e_connector = MAP_modified_Gaussian_spring_potential;
      double cutoff_connection = atof(given_cond("cutoff_connection").c_str());
      if (cutoff_connection > 0.0 && cutoff_connection < atof(given_cond("box_dimension").c_str()))
        {
          // given_POT.force_variables[7] = cutoff_connection;
          given_POT.force_variables.cutoff_connection = cutoff_connection;
          given_POT.PDF_connector = MAP_cutoff_modified_Gaussian_Boltzmann;
	  if(given_cond("association_probability") == "minimum_R0_Boltzmann")
	    {
	      given_POT.PDF_connector = MAP_minimum_R0_cutoff_modified_Gaussian_Boltzmann;
	    }
        }
      else
        {
          given_POT.PDF_connector = MAP_modified_Gaussian_Boltzmann;
	  if(given_cond("association_probability") == "minimum_R0_Boltzmann")
	    {
	      given_POT.PDF_connector = MAP_minimum_R0_modified_Gaussian_Boltzmann;
	    }
	  
        }
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
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str())/atof(given_cond("Rt").c_str());
      // from now on, the Rt = tau_0/tau_B instead of tau_0/tau_R
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      given_POT.force_variables.delta_t0 = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      
    }
  else if (given_cond("transition_probability")=="CUTOFF_DISSOCIATION")
    {
      given_POT.transition = KINETICS::cutoff_dissociation_probability;
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      given_POT.force_variables.delta_t0 = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      
    }
  else if (given_cond("transition_probability")=="MINIMUM_R0_DISSOCIATION")
    {
      given_POT.transition = KINETICS::minimum_R0_dissociation_probability;
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      given_POT.force_variables.delta_t0 = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));

    }
  else if (given_cond("transition_probability")=="FIRST_ORDER")
    {
      given_POT.transition = KINETICS::FIRST_ORDER::dissociation_probability;
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str())/atof(given_cond("Rt").c_str());
      // from now on, the Rt = tau_0/tau_B instead of tau_0/tau_R      
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      given_POT.force_variables.delta_t0 = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));
      
    }
  else if (given_cond("transition_probability")=="Modified_Gaussian")
    {
      given_POT.transition = KINETICS::dissociation_probability_equal_modified_gaussian;
      // given_POT.force_variables[6] = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));  
      given_POT.force_variables.delta_t0 = atof(given_cond("dt/tauR").c_str())*atof(given_cond("N_steps_block").c_str()) / (atof(given_cond("repulsion_coefficient").c_str())*atof(given_cond("Rt").c_str()));      
    }

  else
    {
      // this is default setting
      given_POT.transition = KINETICS::UNIFORM::transition_probability;
    }

  if(given_cond("association_probability") == "EQUAL_CUTOFF_RANGE")
    {
      double cutoff_connection = atof(given_cond("cutoff_connection").c_str());
      // given_POT.force_variables[7] = cutoff_connection;
      given_POT.force_variables.cutoff_connection = cutoff_connection;
      
      given_POT.PDF_connector = MAP_cutoff_equal_probability;
    }

  given_POT.zeta_particle = FORCE::DEFAULT::dimensionless_friction;
  if(given_cond("friction_junction")=="LOOP_DISSOCIATION_TIME")
    {
      // given_POT.force_variables[9] = atof(given_cond("Rt").c_str())*atof(given_cond("repulsion_coefficient").c_str());
      given_POT.force_variables.modified_friction_tau0 = atof(given_cond("Rt").c_str())*atof(given_cond("repulsion_coefficient").c_str());
      
      given_POT.zeta_particle = FORCE::NAPLE::MC_ASSOCIATION::friction_LOOP_DISSOCIATION_TIME;
    }
  
  
  if(given_cond("chain_selection")=="METROPOLIS")
    {
      // given_POT.force_variables[6] = atof(given_cond("energy_barrier").c_str());
      given_POT.force_variables.energy_barrier = atof(given_cond("energy_barrier").c_str());
      
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

  
  if(given_cond("repulsion_type")=="SOFT_REPULSION_P2")
    {
      /* 
         SOFT_REPULSION_P2 means C_rep(p) = C_0*p^2,
         which is the default function at this moment
       */
      given_POT.force_variables.repulsion_coefficient_base = atof(given_cond("repulsion_coefficient_base").c_str());
      /*
      // The following code is commented because the function is deatched from repulsion_type but it is related with the repulsion_coupling.
      // To be specific, SOFT_REPULSION_P2 means the repulsive coefficient is given by C_rep = C_0 p^2 where p is aggregation number (number of chain ends per micelle), which, however, has individual role compare with how coupled between different particles.
      
      given_POT.f_repulsion_coupling = FORCE::NAPLE::SIMPLE_REPULSION::geometrical_mean_repulsive_coefficient;
      */

    }
  else
    {
      /*
        default option. The definition will be called "SOFT_REPULSION". However, for compatibility with the previous inp files, if "repulsion_type" is not defined or is not given by "SOFT_REPULSION_P2", it will be regarded as default test.
       */
      given_POT.f_repulsion_coupling = FORCE::NAPLE::SIMPLE_REPULSION::pre_averaged_repulsive_coefficient;
    }

  if(given_cond("repulsion_coupling")=="PRE_AVERAGE")
    {
      // default option. no further change is necessary
    }
  else if(given_cond("repulsion_coupling")=="GEOMETRIC_MEAN")
    {
      
    }
  else
    {
      printf("ERR: no option for repulsion_coupling.\n");
      return -1;
    }
  
  
  return 0;
}

