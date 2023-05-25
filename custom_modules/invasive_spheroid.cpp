/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./invasive_spheroid.h"
#include "./extracellular_matrix.h"
#include <chrono>  // for high_resolution_clock - https://www.pluralsight.com/blog/software-development/how-to-measure-execution-time-intervals-in-c--
// Cell_Definition fibroblast; 
// Cell_Definition cancer_cell; 
ECM ecm;
unsigned long long int counter=0; // counter for calculating average for the ad hoc timing I am doing ... 
int time_total = 0;

// void create_cell_types( void )
// {
// 	// // use the same random seed so that future experiments have the 
// 	// // same initial histogram of oncoprotein, even if threading means 
// 	// // that future division and other events are still not identical 
// 	// // for all runs 
// 	// SeedRandom(0); 
	
// 	// housekeeping 
	
// 	initialize_default_cell_definition();
// 	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	
// 	// turn the default cycle model to live, 
// 	// so it's easier to turn off proliferation

// 	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
// 	// Make sure we're ready for 2D
	
// 	cell_defaults.functions.set_orientation = up_orientation; 
// 	cell_defaults.phenotype.geometry.polarity = 1.0; 
// 	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
// 	// use default proliferation and death 
	
// 	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
// 	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
// 	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
// 	int necrosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 

// 	// For strict ECM invasion testing, why have on death and birth at all??

// 	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.0;
//     cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0;
// 	cell_defaults.phenotype.death.rates[necrosis_index] = 0.0;

// 	cell_defaults.parameters.o2_proliferation_saturation = 38.0;
// 	cell_defaults.parameters.o2_reference = 38.0;
	
	
// 	// set default uptake and secretion 
// 	// oxygen 
// 	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
// 	cell_defaults.phenotype.secretion.uptake_rates[0] = parameters.doubles("oxygen_uptake"); 
// 	std::cout<<cell_defaults.phenotype.secretion.uptake_rates[0]<<std::endl;
// 	cell_defaults.phenotype.secretion.saturation_densities[0] = 38; 

// 	// Fields for fibroblast-cancer diffusing signal based communication. Commented out 03.11.19

// 	/*

// 	cell_defaults.phenotype.secretion.secretion_rates[1] = 0;
// 	cell_defaults.phenotype.secretion.uptake_rates[1] = 0;
// 	cell_defaults.phenotype.secretion.saturation_densities[1] = 1;

// 	cell_defaults.phenotype.secretion.secretion_rates[2] = 0;
// 	cell_defaults.phenotype.secretion.uptake_rates[2] = 0;
	
// 		cell_defaults.phenotype.secretion.saturation_densities[2] = 1;

// 	*/	

// 	// For phenotype switching if using
	
// 	// cell_defaults.functions.update_phenotype = switching_phenotype_model;
	
// 	cell_defaults.name = "cancer cell"; 
// 	cell_defaults.type = 0; 
	
// 	// set default motility parameters (even for when off)
	
// 	cell_defaults.phenotype.motility.is_motile = true;
// 	// consider what "best" persistence time would be, given the voxel dimensions. 
// 	cell_defaults.phenotype.motility.persistence_time = parameters.doubles("default_persistence_time"); //10.0; // Voxels are 20 um in all dimensions. Given a top speed of 0.5 um/min, cells will likely be in one voxel for 10 minutes or more. So update of 10 isn't bad. Should consider "best" number later. 
// 	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("default_cell_speed");
// 	cell_defaults.phenotype.motility.restrict_to_2D = true; 
// 	cell_defaults.phenotype.motility.migration_bias = 1.0;// completely random - setting in update_migration_bias - might wnat to call that immediately thing

// 	// add custom data 
// 	cell_defaults.custom_data.add_variable( "min ECM motility density", "dimensionless", parameters.doubles( "rho_L") );  // Minimum ECM density required for cell motility
// 	cell_defaults.custom_data.add_variable( "max ECM motility density", "dimensionless", parameters.doubles( "rho_H") );  // Maximum ECM density allowing cell motility
// 	cell_defaults.custom_data.add_variable( "ideal ECM motility density", "dimensionless", parameters.doubles( "rho_I") );  // Ideal ECM density cell motility
// 	cell_defaults.custom_data.add_variable( "max speed", "micron/min" , parameters.doubles( "default_cell_speed") ); // Maximum migration speed
// 	cell_defaults.custom_data.add_variable( "chemotaxis bias", "dimensionless", parameters.doubles( "default_chemotaxis_bias") ); 
// 	cell_defaults.custom_data.add_variable( "ECM sensitivity", "dimensionless", parameters.doubles("default_ECM_sensitivity") );
// 	cell_defaults.custom_data.add_variable( "hypoxic switch value" , "mmHg", 10 );
// 	cell_defaults.custom_data.add_variable( "target ECM density", "dimensionless", parameters.doubles( "default_ECM_density_target") ); 
// 	cell_defaults.custom_data.add_variable( "ECM_production_rate", "1/min", parameters.doubles( "default_ECM_production_rate") );
// 	cell_defaults.custom_data.add_variable( "Base hysteresis bias", "dimensionless", parameters.doubles( "default_hysteresis_bias") );
// 	cell_defaults.custom_data.add_variable( "previous anisotropy", "dimensionless", 0 );
// 	cell_defaults.custom_data.add_variable( "Anisotropy increase rate", "1/min", parameters.doubles( "anisotropy_increase_rate") );
// 	cell_defaults.custom_data.add_variable( "Fiber realignment rate", "1/min", parameters.doubles( "fiber_realignment_rate") );
	

// 	// <unit_test_setup description="Specifies cell parameters for consistent unit tests of ECM influenced mechanics and mechanics influence on ECM - sets adhesion to 1.25, repulsion to 25, and speed to 1.0" type="bool">cells at left boundary/march</unit_test_setup>

// 	if( parameters.ints("unit_test_setup") == 1)
// 	{
// 		cell_defaults.phenotype.motility.persistence_time = 10.0; 
// 		cell_defaults.phenotype.motility.migration_speed = 1.0;
// 		cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
// 		cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
// 		cell_defaults.phenotype.secretion.uptake_rates[0] = 0.0;
// 		cell_defaults.custom_data.add_variable( "max speed", "micron/min" , 1.0 ); // Maximum migration speed
// 		cell_defaults.custom_data.add_variable( "min ECM motility density", "dimensionless", 0.0 );  // Minimum ECM density required for cell motility
// 		cell_defaults.custom_data.add_variable( "max ECM motility density", "dimensionless", 1.0 );  // Maximum ECM density allowing cell motility
// 		cell_defaults.custom_data.add_variable( "ideal ECM motility density", "dimensionless", 0.5 );  // Ideal ECM density cell motility
// 		cell_defaults.custom_data.add_variable( "target ECM density", "dimensionless", 0.5 ); 

// 		std::cout<< "running unit test setup cancer cell"<<std::endl;

// 	}

// 	else if ( parameters.ints("unit_test_setup") == 0)
// 	{
// 		std::cout<<"not in unit test mode"<<std::endl;
// 	}

// 	else
// 	{
// 		std::cout<<"WARNING!!!!! Cell parameters not set correctly - unit test set up must either be true or false!!!!"<<std::endl;
// 	}
	
	
// 	// fibroblast cells 
	
// 	fibroblast = cell_defaults;
// 	fibroblast.name = "fibroblast cell"; 
// 	fibroblast.type = 1; 

// 	// Temperarily eliminating fibroblast/cancer cell signal

// 	// Obviously missing - add later
    
// 	// 10% proliferation 
//     // fibroblast.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.10;

// 	/*******************************************For "march" simulation****************************************/

// 	fibroblast.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.0;
//     fibroblast.phenotype.death.rates[apoptosis_index] = 0.0;
// 	fibroblast.phenotype.death.rates[necrosis_index] = 0.0;
	
    
// 	// Temperarily eliminating fibroblast/cancer cell signal	
	
// 	// Obviously missing - add later

// 	// turn on motility 
// 	fibroblast.phenotype.motility.is_motile = parameters.bools("fibroblast_motility_mode"); 
	
//     fibroblast.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("fibroblast_adhesion");
// 	// std::cout<<fibroblast.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
    
// 	fibroblast.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("fibroblast_repulsion");
//     // std::cout<<fibroblast.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;

// 	// Temperarily eliminating fibroblast/cancer cell signal	

// 	//    fibroblast.phenotype.secretion.secretion_rates[1] = 50; // fibroblast signal

//     // modify ECM
    
// 	if ( parameters.strings("ecm_update_model") == "ecm_update_from_cell_motility_vector")
//     	fibroblast.functions.custom_cell_rule = ecm_update_from_cell_motility_vector; // Only fibroblasts can modify ECM (phenotype -> ECM)

// 	else if( parameters.strings("ecm_update_model") == "ecm_update_from_cell_velocity_vector")
// 	{
// 		fibroblast.functions.custom_cell_rule = ecm_update_from_cell_velocity_vector; // Only fibroblasts can modify ECM (phenotype -> ECM)
// 	}
	
// 	else
// 	{
// 		std::cout<<"no reorientation model specified!!@!!!!!! Halting!!!!!!"<<std::endl;
// 		abort();
// 		return;
// 	}
	

// 	// set functions
	
// 	fibroblast.functions.update_migration_bias = chemotaxis_oxygen;//rightward_deterministic_cell_march; Use rightward deterministic march for march test. Set fibroblast fraction to 1.0.
	
//     fibroblast.functions.update_phenotype = NULL; // fibroblast_phenotype_model;

// 	if( parameters.ints("unit_test_setup") == 1)
// 	{
// 		fibroblast.phenotype.motility.persistence_time = 10.0;
// 		fibroblast.phenotype.motility.migration_speed = 0.50;
// 		fibroblast.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
// 		fibroblast.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
// 		fibroblast.custom_data.add_variable( "max speed", "micron/min" , 1.0 ); // Maximum migration speed
// 		fibroblast.custom_data.add_variable( "min ECM motility density", "dimensionless", 0.0 );  // Minimum ECM density required for cell motility
// 		fibroblast.custom_data.add_variable( "max ECM motility density", "dimensionless", 1.0 );  // Maximum ECM density allowing cell motility
// 		fibroblast.custom_data.add_variable( "ideal ECM motility density", "dimensionless", 0.5 );  // Ideal ECM density cell motility
// 		fibroblast.custom_data.add_variable( "target ECM density", "dimensionless", 0.5 ); 
// 		if (parameters.ints("march_unit_test_setup") == 1){
// 		fibroblast.functions.update_migration_bias = rightward_deterministic_cell_march;
// 		}

// 		// std::cout<< "running unit test setup fibroblast"<<std::endl;

// 	}
	
// 	// cancer cells

// 	cancer_cell = cell_defaults;
// 	cancer_cell.name = "cancer cell"; 
// 	cancer_cell.type = 2;
    
//     cancer_cell.functions.update_phenotype = cancer_cell_phenotype_model;

// 	cancer_cell.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("cancer_cell_adhesion");
// 	std::cout<<cancer_cell.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
// 	cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("cancer_cell_repulsion");
//    	std::cout<<cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
// 	cancer_cell.phenotype.motility.is_motile = parameters.bools("cancer_cell_motility_mode");
	
// 	// Selecting cell-ECM interaction wrt to hyteriss
	
// 	if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "cancer cell chemotaxis/no cancer cell hysteresis" || parameters.ints("unit_test_setup") == 1)
// 	{
// 		cancer_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis;
// 		std::cout<<"I selected cancer cell chemotaxsis" << std::endl;
// 	}

// 	else if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "cancer cell hysteresis/no cancer cell chemotaxis")
// 	{
// 		cancer_cell.functions.update_migration_bias = ECM_informed_motility_update_model_w_memory;
// 		// cancer_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
// 		// void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
// 		std::cout<<"I selected cancer cell hysteresis" << std::endl;
// 	}

// 	else if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "cancer cell chemotaxis with variable cancer cell speed")
// 	{
// 		cancer_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
// 		// cancer_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
// 		// void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
// 		std::cout<<"I selected cancer cell chemotaxis with variable cancer cell speed" << std::endl;
// 	}

// 	else
// 	{
// 		std::cout<<"WARNING: NO CELL-ECM MODEL SPECIFIED. FIX THIS!!!"<<std::endl;
// 		std::cout<<"Halting program!!!"<<std::endl;
// 		abort();
// 		return;
// 	}

// 	// Why do these lines not overwrite the cell defaults?

// 	if( parameters.ints("unit_test_setup") == 1)
// 	{
// 		cancer_cell.phenotype.motility.persistence_time = 10.0;
// 		cancer_cell.phenotype.motility.migration_speed = 1.0;
// 		cancer_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
// 		cancer_cell.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
// 		cancer_cell.custom_data.add_variable( "max speed", "micron/min" , 1.0 ); // Maximum migration speed
// 		cancer_cell.custom_data.add_variable( "min ECM motility density", "dimensionless", 0.0 );  // Minimum ECM density required for cell motility
// 		cancer_cell.custom_data.add_variable( "max ECM motility density", "dimensionless", 1.0 );  // Maximum ECM density allowing cell motility
// 		cancer_cell.custom_data.add_variable( "ideal ECM motility density", "dimensionless", 0.5 );  // Ideal ECM density cell motility
// 		cancer_cell.custom_data.add_variable( "target ECM density", "dimensionless", 0.5 ); 

// 		std::cout<< "running unit test setup cancer cell"<<std::endl;

// 	}
    
// 	std::cout<<"cancer cell cell migration speed "<<cancer_cell.phenotype.motility.migration_speed <<std::endl;

// 	// Temperarily eliminating fibroblast/cancer cell signal

//    	cancer_cell.phenotype.secretion.secretion_rates[1] = 50; // cancer cell signal
// 	cancer_cell.phenotype.secretion.saturation_densities[1] = 1; // cancer cell signal
    
// 	// Temperarily eliminating fibroblast/cancer cell signal
// 	return; 
// }	

void create_cell_types( void )
{
	// set the random seed 
	// SeedRandom( parameters.ints("random_seed") );  
	// SeedRandom(0);  
	
	/* 
	   Put any modifications to default cell definition here if you 
	   want to have "inherited" by other cell types. 
	   
	   This is a good place to set default functions. 
	*/ 
	
	initialize_default_cell_definition(); 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.volume_update_function = standard_volume_update_function;
	cell_defaults.functions.update_velocity = standard_update_cell_velocity;

	cell_defaults.functions.update_migration_bias = NULL; 
	cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based; 
	cell_defaults.functions.custom_cell_rule = NULL; 
	cell_defaults.functions.contact_function = NULL; 
	
	cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	cell_defaults.functions.calculate_distance_to_membrane = NULL; 

    cell_defaults.phenotype.motility.migration_speed = parameters.doubles("default_cell_speed");  //rwh
	
	/*
	   This parses the cell definitions in the XML config file. 
	*/
	
	initialize_cell_definitions_from_pugixml(); 

	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	build_cell_definitions_maps(); 

	/*
	   This intializes cell signal and response dictionaries 
	*/

	setup_signal_behavior_dictionaries(); 	

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 

    Cell_Definition* fibroblast = find_cell_definition("fibroblast");	
	Cell_Definition* cancer_cell = find_cell_definition("cancer cell");	

	if ( parameters.strings("ecm_update_model") == "ecm_update_from_cell_motility_vector")
    {
    	fibroblast->functions.custom_cell_rule = ecm_update_from_cell_motility_vector; // Only leaders can modify ECM (phenotype -> ECM)
    }
	else if( parameters.strings("ecm_update_model") == "ecm_update_from_cell_velocity_vector")
	{
		fibroblast->functions.custom_cell_rule = ecm_update_from_cell_velocity_vector; // Only leaders can modify ECM (phenotype -> ECM)
	}
	else
	{
		std::cout<<"no reorientation model specified!!@!!!!!! Halting!!!!!!"<<std::endl;
		abort();
		return;
	}
	
	fibroblast->functions.update_migration_bias = chemotaxis_oxygen;//rightward_deterministic_cell_march; Use rightward deterministic march for march test. Set leader fraction to 1.0.
	
    fibroblast->functions.update_phenotype = NULL; // leader_cell_phenotype_model;

    //--------- now follower:

    cancer_cell->functions.update_phenotype = NULL;// follower_cell_phenotype_model;

// <cell_motility_ECM_interaction_model_selector type="string" units="" description="follower chemotaxis/no follower hysteresis, follower hysteresis/no follower chemotaxis">follower chemotaxis/no follower hysteresis<

    // rwh: doing this one:
    if ( parameters.strings("cell_motility_ECM_interaction_model_selector") == "follower chemotaxis/no follower hysteresis" || parameters.ints("unit_test_setup") == 1)
	{
		cancer_cell->functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis;
		std::cout<<"I selected follower chemotaxsis" << std::endl;   // <------ rwh
	}
	else if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "follower hysteresis/no follower chemotaxis")
	{
		cancer_cell->functions.update_migration_bias = ECM_informed_motility_update_model_w_memory;
		// follower_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
		// void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
		std::cout<<"I selected follower hysteresis" << std::endl;
	}
	else if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "follower chemotaxis with variable follower speed")
	{
		cancer_cell->functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
		// follower_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
		// void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
		std::cout<<"I selected follower chemotaxis with variable follower speed" << std::endl;
	}
	else
	{
		std::cout<<"WARNING: NO CELL-ECM MODEL SPECIFIED. FIX THIS!!!"<<std::endl;
		std::cout<<"Halting program!!!"<<std::endl;
		abort();
		return;
	}


	/*
	   This builds the map of cell definitions and summarizes the setup. 
	*/
		
	display_cell_definitions( std::cout ); 
	
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

// void setup_microenvironment( void )
// {

// 	if(parameters.ints("unit_test_setup")==1)
// 	{
// 		default_microenvironment_options.calculate_gradients = false; 
// 	}
	
// 	else
// 	{
// 		default_microenvironment_options.calculate_gradients = true; 
// 	}

// 	// let BioFVM use oxygen as the default 
	
// 	default_microenvironment_options.use_oxygen_as_first_field = true; 
	

// 	// Temperarily eliminating fibroblast/cancer cell signal (except here)
    
// 	// 50 micron length scale 
//     // microenvironment.add_density( "fibroblast signal", "dimensionless", 1e5 , 1 );
//     microenvironment.add_density( "inflammatory_signal", "dimensionless", 1e5 , 1 );

// 	// Temperarily eliminating fibroblast/cancer cell signal	
	
// 	// set Dirichlet conditions 
	
// 	default_microenvironment_options.outer_Dirichlet_conditions = true;

// 	std::vector<double> bc_vector; 
// 	bc_vector = { 38.0, 0.0};
// 	// bc_vector = { 38.0 , 0.0, 0.0};  // 5% o2 , fibroblast signal, cancer cell signal

// 	// if(parameters.ints("unit_test_setup")==1 && parameters.ints("march_unit_test_setup") == 0)
// 	// {

// 	// 	bc_vector = { 38.0 }; // 5% o2 , fibroblast signal, cancer cell signal
// 	// 	default_microenvironment_options.X_range[0] = -500.0;
// 	// 	default_microenvironment_options.X_range[1] = 500.0;
// 	// 	default_microenvironment_options.Y_range[0] = -500.0;
// 	// 	default_microenvironment_options.Y_range[1] = 500.0;

// 	// }

// 	// else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
// 	// {
// 	// 	bc_vector = { 38.0}; // 5% o2 , fibroblast signal, cancer cell signal
// 	// 	default_microenvironment_options.X_range[0] = -500.0;
// 	// 	default_microenvironment_options.X_range[1] = 500.0;
// 	// 	default_microenvironment_options.Y_range[0] = -500.0;
// 	// 	default_microenvironment_options.Y_range[1] = 500.0;
// 	// }

// 	// else if(parameters.ints("unit_test_setup") == 0 && parameters.ints("march_unit_test_setup") == 0)
// 	// {
// 	// 	bc_vector = { 38.0 };  // 5% o2 , fibroblast signal, cancer cell signal
// 	// }

// 	// else
// 	// {
// 	// 	std::cout<<"ECM density and anisotropy not set correctly!!!! FIX!!!!!!!!!"<<std::endl;
// 	// 	std::cout<<"Halting!"<<std::endl;
// 	// 	abort();
// 	// 	return;
// 	// }
	
	
// 	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
    
// 	// Temperarily eliminating fibroblast/cancer cell signal	
// 	// default_microenvironment_options.Dirichlet_condition_vector[1] = 0; // normoxic conditions
// 	// default_microenvironment_options.Dirichlet_condition_vector[2] = 0; // normoxic conditions
    
// 	initialize_microenvironment(); 

// 	// ecm.ecm_mesh.resize(default_microenvironment_options.X_range[0], default_microenvironment_options.X_range[1] , 
// 	// 	default_microenvironment_options.Y_range[0], default_microenvironment_options.Y_range[1],default_microenvironment_options.Z_range[0], default_microenvironment_options.Z_range[1], \
// 	// 	parameters.doubles("ECM_dx"), parameters.doubles("ECM_dy"),parameters.doubles("ECM_dz"));
// 	// ecm.resize_ecm_units_from_ecm_mesh();


// 	// ecm.ecm_mesh.display_information(std::cout );

// 	// std::cout<<ecm.ecm_mesh.nearest_voxel_index(position)<<std::endl;
// 	// std::cout<<microenvironment.mesh.nearest_voxel_index(position)<<std::endl;
// 	// std::cout<<" hit Enter to continue:"<<std::flush;
// 	// std::cin.get();

// 	microenvironment.decay_rates[0] = parameters.doubles("chemotactic_substrate_decay_rate");

// 	// Trying to set the chemical gradient to be a starburst. Using the same code snippets as the ECM orientation. 

// 	// Catching if chemical field gradient not specified for unit testing - if chemical_field_setup_specified_bool is true, there is a string match in the chemical field setup, else, the field is not specified. Program halts in this condition if unit testing IS specified and displays message requesting user specify the chemical field set up.

// 	bool chemical_field_setup_specified_bool = (parameters.strings("chemical_field_setup")== "starburst" ||  parameters.strings("chemical_field_setup")== "vertical up" || parameters.strings("chemical_field_setup") == "horizontal right" || parameters.strings("chemical_field_setup") == "angle" || parameters.strings("chemical_field_setup") == "none");

// 	// std::cout<<"Chemical field specified? "<<chemical_field_setup_specified_bool<<std::endl;

// 	if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "starburst")
// 	{

// 		for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
// 		{
// 			std::vector<double> position = microenvironment.mesh.voxels[i].center; 
// 			microenvironment.gradient_vector(i)[0] = { position[0],position[1],0}; 
// 			normalize(&microenvironment.gradient_vector(i)[0]);
// 			// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
// 		}
// 	}

// 	if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "vertical up")
// 	{

// 		for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
// 		{
// 			// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
// 			microenvironment.gradient_vector(i)[0] = { 0,1,0}; 
// 			normalize(&microenvironment.gradient_vector(i)[0]);
// 			// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
// 		}
// 	}

// 	if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "horizontal right")
// 	{

// 		for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
// 		{
// 			// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
// 			microenvironment.gradient_vector(i)[0] = { 1,0,0}; 
// 			normalize(&microenvironment.gradient_vector(i)[0]);
// 			// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
// 		}
// 	}

// 	if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "angle")
// 	{
		

// 		for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
// 		{
// 			// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
// 			microenvironment.gradient_vector(i)[0] = { cos ( parameters.doubles("angle_of_chemical_field_gradient") * PhysiCell_constants::pi/180) , sin ( parameters.doubles("angle_of_chemical_field_gradient") * PhysiCell_constants::pi/180),0}; 
// 			normalize(&microenvironment.gradient_vector(i)[0]);
// 			// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
// 		}
// 	}

// 	if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "none")
// 	{

// 		for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
// 		{
// 			// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
// 			microenvironment.gradient_vector(i)[0] = { 0,0,0}; 
// 		}

// 	}

// 	else if( parameters.ints("unit_test_setup") == 1 && chemical_field_setup_specified_bool == 0)
// 	{
// 		std::cout<<"WARNING: NO CHEMICAL FIELD ORIENTATION SPECIFIED for unit testing. FIX THIS!!!"<<std::endl;
// 		std::cout<<"Halting program!!!"<<std::endl;
// 		abort();
// 		return;
// 	}


// 	// run to get a decent starting conditoin (refers to code no longer present but will be added back for fibroblast-cancer cell signaling models)
	
// 	// now, let's set the fibroblast signal to 1, so we don't hvae early swiching 
// 	/*
// 	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
// 	{
// 		microenvironment.density_vector(i)[1] = 1.0; 
// 	}
// 	*/

// 	// set up ECM density and anisotropy profile as needed
// 	// int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
// 	// int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 

// 	/*for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
// 	{
// 		std::vector<double> position = microenvironment.mesh.voxels[n].center; 
// 		if( fabs( position[0] ) > 200 || fabs( position[1] ) > 200 )
// 		{
// 			microenvironment(n)[ECM_density_index] = 0.0; 
// 			microenvironment(n)[ECM_anisotropy_index] = 1.0; 
// 		}
// 	}*/
	
// 	return; 
// }	

void setup_extracellular_matrix( void )
{
	// DEPENDS ON MICROENVIRONMENT - CALL SETUP MICROENVIRONEMNT FIRST!!!!!

	ecm.ecm_mesh.resize(default_microenvironment_options.X_range[0], default_microenvironment_options.X_range[1] , 
	default_microenvironment_options.Y_range[0], default_microenvironment_options.Y_range[1],default_microenvironment_options.Z_range[0], default_microenvironment_options.Z_range[1], \
	parameters.doubles("ECM_dx"), parameters.doubles("ECM_dy"),parameters.doubles("ECM_dz"));
	ecm.resize_ecm_units_from_ecm_mesh();

	ecm.ecm_mesh.display_information(std::cout );

	// set up ECM alignment 

	// <ECM_orientation_setup description="Specifies the initial ECM orientation: random, circular, starburt, oriented to the right, or oriented to the top" type="string" units="NA">circular</ECM_orientation_setup> parameters.string( "ECM_orientation_setup")
	
	for( int n = 0; n < ecm.ecm_mesh.voxels.size() ; n++ )
	{
		
		// ############################# Density and anisotropy ####################################
		if(parameters.ints("unit_test_setup")==1 && parameters.ints("march_unit_test_setup") == 0)
		{

			ecm.ecm_voxels[n].density = 0.5;
			ecm.ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");

		}

		else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
		{
			
			ecm.ecm_voxels[n].density = 0.5;
			ecm.ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");
			
		}

		else if(parameters.ints("unit_test_setup") == 0 && parameters.ints("march_unit_test_setup") == 0)
		{
			ecm.ecm_voxels[n].density = parameters.doubles("initial_ECM_density");
			ecm.ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");
			
		}

		else
		{
			std::cout<<"ECM density and anisotropy not set correctly!!!! FIX!!!!!!!!!"<<std::endl;
			std::cout<<"Halting!"<<std::endl;
			abort();
			return;
		}
		
		// ############################# Alignment ####################################
		// ############## CAN ALSO BE USED TO MAKE MORE COMPLEX PATTERNS in density and anisotropy ##############
		// For random 2-D initalization 
		if(parameters.strings( "ECM_orientation_setup") == "random")
		{
			double theta = 6.2831853071795864769252867665590 * uniform_random(); 
			// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
			// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
			// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
			ecm.ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0};
		}

		else if(parameters.strings( "ECM_orientation_setup") == "hard_line")
		{
			// std::vector<double> position = ecm.ecm_mesh.voxels[n].center; 
			if(ecm.ecm_mesh.voxels[n].center[1] < -parameters.doubles("tumor_radius") && ecm.ecm_mesh.voxels[n].center[1] > -parameters.doubles("tumor_radius") - 60)
			{
				ecm.ecm_voxels[n].ecm_fiber_alignment = {1.0, 0.0, 0.0};
				ecm.ecm_voxels[n].density = 1.0;
			}

			else
			{
				double theta = 6.2831853071795864769252867665590 * uniform_random(); 
				// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
				// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
				// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
				ecm.ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0};
			}
		}

		// for starburst initialization 
		else if(parameters.strings( "ECM_orientation_setup") == "starburst")
		{
			std::vector<double> position = ecm.ecm_mesh.voxels[n].center; 
			normalize( &position ); 
			ecm.ecm_voxels[n].ecm_fiber_alignment =  { position[0],position[1],0}; // oriented out (perpindeicular to concentric circles)
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}

		// for circular initialization 
		else if(parameters.strings( "ECM_orientation_setup") == "circular")
		{
			std::vector<double> position = ecm.ecm_mesh.voxels[n].center;; 
			normalize( &position );
			ecm.ecm_voxels[n].ecm_fiber_alignment =  { position[1],-position[0],0}; // oriented in cirlce
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}

		else if(parameters.strings( "ECM_orientation_setup") == "horizontal")
		{
			ecm.ecm_voxels[n].ecm_fiber_alignment =  { 1.0, 0.0, 0.0}; 
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}

		else if(parameters.strings( "ECM_orientation_setup") == "vertical")
		{
			ecm.ecm_voxels[n].ecm_fiber_alignment=  { 0.0, 1.0, 0.0}; 
			normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
		}

		else if(parameters.strings( "ECM_orientation_setup") == "split")
		{
			std::vector<double> position = ecm.ecm_mesh.voxels[n].center; 
			normalize( &position ); 

			if(position[1]<=0)
			{
				ecm.ecm_voxels[n].ecm_fiber_alignment =  { 1,-1,0}; // oriented out (perpindeicular to concentric circles)
				normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);
			}

			else
			{
				ecm.ecm_voxels[n].ecm_fiber_alignment =  { 1,1,0}; // oriented out (perpindeicular to concentric circles)
				normalize(&ecm.ecm_voxels[n].ecm_fiber_alignment);				
			}

		}

		else
		{
			std::cout<<"WARNING: NO ECM ORIENTATION SPECIFIED. FIX THIS!!!"<<std::endl;
			std::cout<<"Halting program!!!"<<std::endl;
			abort();
			return;
		}


		
	}


}

void run_biotransport( double t_max ) // used to set up initial chemical conditions to prevent unwanted phenotype switching due to starting the simulation
{
	std::cout << "working on initial conditions .. " << std::endl; 
	double t = 0.0;
	
	// make sure the associated cell has the correct rate vectors 
	
	for( int i=0 ; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i];
		
		// Commented out for currently unknown reason.

		/* pCell->secretion_rates = &secretion_rates; 
		pCell->uptake_rates = &uptake_rates; 
		pCell->saturation_densities = &saturation_densities;  */
			
		pCell->set_total_volume( pCell->phenotype.volume.total ); 
		pCell->set_internal_uptake_constants( diffusion_dt );
	}
	
	while( t < t_max )
	{
		microenvironment.simulate_diffusion_decay( diffusion_dt );
		#pragma omp parallel for 
		for( int i=0 ; i < (*all_cells).size() ; i++ )
		{
			(*all_cells)[i]->phenotype.secretion.advance( (*all_cells)[i] , (*all_cells)[i]->phenotype , diffusion_dt ) ;
		}
		t += diffusion_dt; 
	}
	
	std::cout << "done!" << std::endl; 
	return; 
}

void set_cell_motility_vectors( void )
{
	for( int i=0 ; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i];
		pCell->update_motility_vector( 20 );
		std::cout<<"Initial Motility vector "<<pCell->phenotype.motility.motility_vector<<std::endl;

	}
}

void alter_cell_uptake_secretion_saturation ( void )
{
	for( int i=0 ; i < (*all_cells).size() ; i++ )
		{
			Cell* pCell = (*all_cells)[i];
			
			// This assumes that the first substrate is the substrate of interest.
			pCell->phenotype.secretion.secretion_rates[0] = 0; 
			pCell->phenotype.secretion.uptake_rates[0] = 0; 
			pCell->phenotype.secretion.saturation_densities[0] = 0; 
		}

}
void setup_tissue( void )
{	
	// Setting seed so cells always start with same initial configuration
	SeedRandom(0);
	static Cell_Definition* fibroblast = find_cell_definition("fibroblast");	
	static Cell_Definition* cancer_cell = find_cell_definition("cancer cell");	

	if (parameters.ints("march_unit_test_setup") == 0)
	{

		if(parameters.strings("cell_setup") == "single")
		{
			Cell* pC;
			pC = create_cell(*fibroblast);
			pC->assign_position(0.0, 0.0, 0.0);
		}

		else if(parameters.strings("cell_setup") == "random")
		{
			std::cout<<"string worked"<<std::endl;
		
			/*******************************************Random initialization****************************************/
			Cell* pC;
			
			for( int n = 0 ; n < 200 ; n++ )
			{
				pC = create_cell(); 
				pC->assign_position( -450 + 900*UniformRandom() , -450 + 900*UniformRandom() , 0.0 );
			}

			std::cout<<"Cell's placed randomly on domain - this function uses a HARD CODED domain size!!!! WARNING!!!!!"<<std::endl;
			std::cout<<" hit Enter to continue:"<<std::flush;
			std::cin.get();

		}
		
		/******************************************2D Spheroid initialization***************************************/
		
		else if(parameters.strings("cell_setup") == "lesion")
		{
			// place a cluster of tumor cells at the center
		
			//Get tumor radius from XML parameters

			double tumor_radius; 

			if( parameters.ints("unit_test_setup") == 1)
			{
				tumor_radius = 150;

			}

			else
			{
				tumor_radius = parameters.doubles("tumor_radius");
			}
			
			

			// these lines produce automatically calcuated equilibirum spacing for intiailizing cells, even when changing adh-rep parameters.

			double cell_radius = cell_defaults.phenotype.geometry.radius;
			double relative_maximum_adhesion_distance = cell_defaults.phenotype.mechanics.relative_maximum_adhesion_distance;
			double sqrt_adhesion_to_repulsion_ratio;
			
			if(parameters.doubles("cancer_cell_repulsion") == 0)
			{
				sqrt_adhesion_to_repulsion_ratio = 0.632455; // value for adhesion = 10 and repulsion = 25.0 - "parameter set 21"
			}

			else
			{
				sqrt_adhesion_to_repulsion_ratio = sqrt(parameters.doubles("cancer_cell_adhesion")/parameters.doubles("cancer_cell_repulsion"));
			} 

			double cell_spacing = (1 - sqrt_adhesion_to_repulsion_ratio);
			cell_spacing /= (0.5 * 1/cell_radius - 0.5 * sqrt_adhesion_to_repulsion_ratio/(relative_maximum_adhesion_distance * cell_radius));
			
			Cell* pCell = NULL; 
			
			double x = 0.0;
			double x_outer = tumor_radius; 
			double y = 0.0;

			double fibroblast_fraction;
			
			if( parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 0)
			{
				fibroblast_fraction = 0.0;
				cell_spacing = 1.90 * cell_radius;
			}

			else
			{
				fibroblast_fraction = parameters.doubles("initial_fibroblast_fraction"); // 0.2;
			}

			int n = 0; 
			while( y < tumor_radius )
			{
				x = 0.0; 
				if( n % 2 == 1 )
				{ x = 0.5*cell_spacing; }
				x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
				
				while( x < x_outer )
				{
					if( UniformRandom() < fibroblast_fraction )
					{ pCell = create_cell(*fibroblast); }
					else
					{ pCell = create_cell(*cancer_cell);}
						
					pCell->assign_position( x , y , 0.0 );
					
					if( fabs( y ) > 0.01 )
					{
						if( UniformRandom() < fibroblast_fraction )
						{ pCell = create_cell(*fibroblast); }
						else
						{ pCell = create_cell(*cancer_cell); }
						pCell->assign_position( x , -y , 0.0 );
					}
					
					if( fabs( x ) > 0.01 )
					{ 
						if( UniformRandom() < fibroblast_fraction )
						{ pCell = create_cell(*fibroblast); }
						else
						{ pCell = create_cell(*cancer_cell); }
						pCell->assign_position( -x , y , 0.0 );
						
						if( fabs( y ) > 0.01 )
						{
							if( UniformRandom() < fibroblast_fraction )
							{ pCell = create_cell(*fibroblast); }
							else
							{ pCell = create_cell(*cancer_cell); }
							
							pCell->assign_position( -x , -y , 0.0 );
						}
					}
					x += cell_spacing; 
					
				}
				
				y += cell_spacing * sqrt(3.0)/2.0; 
				n++; 
			}

			std::cout<<"Cell's placed in 2-lesion at center of domain"<<std::endl;

		}


		/************************************Spheroid with fibroblasts***************************************/		

		else if(parameters.strings("cell_setup") == "invasive_spheroid")
		{
			// place a cluster of tumor cells at the center and fibroblasts below it
		
			//Get tumor radius from XML parameters

			double tumor_radius; 

			if( parameters.ints("unit_test_setup") == 1)
			{
				tumor_radius = 150;

			}

			else
			{
				tumor_radius = parameters.doubles("tumor_radius");
			}

			// these lines produce automatically calcuated equilibirum spacing for intiailizing cells, even when changing adh-rep parameters.

			double cell_radius = cell_defaults.phenotype.geometry.radius;
			double relative_maximum_adhesion_distance = cell_defaults.phenotype.mechanics.relative_maximum_adhesion_distance;
			double sqrt_adhesion_to_repulsion_ratio;
			
			if(parameters.doubles("cancer_cell_repulsion") == 0)
			{
				sqrt_adhesion_to_repulsion_ratio = 0.632455; // value for adhesion = 10 and repulsion = 25.0 - "parameter set 21"
			}

			else
			{
				sqrt_adhesion_to_repulsion_ratio = sqrt(parameters.doubles("cancer_cell_adhesion")/parameters.doubles("cancer_cell_repulsion"));
			} 

			double cell_spacing = (1 - sqrt_adhesion_to_repulsion_ratio);
			cell_spacing /= (0.5 * 1/cell_radius - 0.5 * sqrt_adhesion_to_repulsion_ratio/(relative_maximum_adhesion_distance * cell_radius));
			
			Cell* pCell = NULL; 
			
			double x = 0.0;
			double x_outer = tumor_radius; 
			// double center_offset = 0.0;
			double y = 0.0;

			double fibroblast_fraction;
			
			if( parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 0)
			{
				fibroblast_fraction = 0.0;
				cell_spacing = 1.90 * cell_radius;
			}

			else
			{
				fibroblast_fraction = parameters.doubles("initial_fibroblast_fraction"); // 0.2;
			}

			int n = 0; 
			while( y < tumor_radius )
			{
				x = 0.0; 
				if( n % 2 == 1 )
				{ x = 0.5*cell_spacing; }
				x_outer = sqrt( tumor_radius*tumor_radius - y*y ); 
				
				while( x < x_outer )
				{
					if( UniformRandom() < fibroblast_fraction )
					{ pCell = create_cell(*fibroblast); }
					else
					{ pCell = create_cell(*cancer_cell);}
						
					pCell->assign_position( x , y , 0.0 );
					
					if( fabs( y ) > 0.01 )
					{
						if( UniformRandom() < fibroblast_fraction )
						{ pCell = create_cell(*fibroblast); }
						else
						{ pCell = create_cell(*cancer_cell); }
						pCell->assign_position( x , -y , 0.0 );
					}
					
					if( fabs( x ) > 0.01 )
					{ 
						if( UniformRandom() < fibroblast_fraction )
						{ pCell = create_cell(*fibroblast); }
						else
						{ pCell = create_cell(*cancer_cell); }
						pCell->assign_position( -x , y , 0.0 );
						
						if( fabs( y ) > 0.01 )
						{
							if( UniformRandom() < fibroblast_fraction )
							{ pCell = create_cell(*fibroblast); }
							else
							{ pCell = create_cell(*cancer_cell); }
							
							pCell->assign_position( -x , -y , 0.0 );
						}
					}
					x += cell_spacing; 
					
				}
				
				y += cell_spacing * sqrt(3.0)/2.0; 
				n++; 
			}

			std::cout<<"Cell's placed in 2-lesion at center of domain"<<std::endl;

			/***********************************Add in the fibroblasts**************************************/

			int number_of_fibroblasts = parameters.ints("number_of_fibroblasts");
			n =-500.0;
			for(int i=0; i<number_of_fibroblasts; i++)
			{
				pCell = create_cell(*fibroblast);
				pCell->assign_position( n , -600.0, 0.0 );
				std::cout<<"Fibroblast placed at "<<pCell->position<<std::endl;
				n += 100.0;
			}

		}		

		/******************************************3D Spheroid initialization***************************************/

		/*To come later*/

		/************************************Circle of cells at R = 300 initialization***************************************/

		else if(parameters.strings("cell_setup") == "circle of cells")
		{
			double theta2 = 0.0;
			for (int a = 0; a<42; a++)
			{
				Cell* pCell = NULL;
				pCell = create_cell(*cancer_cell); 
				pCell->assign_position( 300 * cos(theta2) , 300 * sin(theta2) , 0.0 );
				theta2 += 0.14959952;
			}

		}

		/************************************Line of cells at y = 0, x > 0 initialization***************************************/

		else if(parameters.strings("cell_setup") == "cells at y = 0")
		{

			Cell* pCell = NULL; 
			int n = default_microenvironment_options.X_range[0] + 10.0; 
			while( n <= default_microenvironment_options.X_range[1] )
			{
				pCell = create_cell(*cancer_cell); 

				// To prevent droping cells in areas of high ECM curvature. 
				while(abs(n) < 70)
				{n = n + 30.0;}

				pCell->assign_position( n , 0.0 , 0.0 );
				n = n + 30.0;
			}
			std::cout<<"Cell's placed at y = 0"<<std::endl;
		}

		/******************************************Line of cells at x = left boundary + 10 initialization***************************************/
		else if(parameters.strings("cell_setup") == "cells at left boundary/march")
		{
			
			Cell* pCell = NULL; 
			int n = default_microenvironment_options.X_range[0] + 10.0; 
			while( n <= default_microenvironment_options.X_range[1] - 10.0 )
			{
				if (parameters.ints("march_unit_test_setup") == 1)
				{pCell = create_cell(*fibroblast);}

				else 
				{pCell = create_cell(*cancer_cell);}
				pCell->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
				n = n + 10.0;
			}
			std::cout<<"Cell's placed at left boundary for march test"<<std::endl;
		}

		else
		{
			std::cout<<"WARNING!!! NO CELL SETUP SPECIFIED. SEE DOCUMENTATION and FIX"<<std::endl;
			std::cout<<"Halting program!!!"<<std::endl;
			abort();
			return;
		}
	}

	else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
	{
		Cell* pCell = NULL; 
		int n = default_microenvironment_options.X_range[0] + 10.0; 
		while( n <= default_microenvironment_options.X_range[1] - 10.0 )
		{
			pCell = create_cell(*fibroblast); 
			pCell->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
			n = n + 10.0;
		}
		std::cout<<"Cell's placed at left boundary for march test"<<std::endl;

	}

	else
	{
		std::cout<<"RUN MODE (TESTING OR NOT TESTING) NOT SPECIFIED!!!!! WARNING!!!!"<<std::endl;
		std::cout<<"Halting program!!!"<<std::endl;
		abort();
		return;
	}
	

	return; 
}

double dot_product( const std::vector<double>& v , const std::vector<double>& w )
{
	double out = 0.0; 
	for( unsigned int i=0 ; i < v.size() ; i++ )
	{ out += ( v[i] * w[i] ); }

	if( fabs(out) < 1e-10)
	{out = 0.0;}

	return out; 
}

double sign_function (double number)
{
	// double sign = 0.0
	if (number<0)
	{ return -1.0;}

	else
	{ return 1.0;}

}

/* To eliminate chemotaxis, set chemotaxis bias to zero. To eliminate ECM influence, set a to 0 (permanently) or ECM senstiivity to zero */

void ECM_informed_motility_update_w_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	
	if(phenotype.death.dead == true)
	{
		
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		std::cout<<2<<std::endl;
		std::cout<<"Cell is dead"<<std::endl;
	}
	// Updates cell bias vector and cell speed based on the ECM density, anisotropy, and fiber direction
	
	// find location of variables and base parameter values
	// static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	// static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 

	    // rwh: use underscores now that they are in the .xml as tags
	static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max_speed" ); 
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis_bias");
	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
	static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min_ECM_motility_density");
    if (min_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int min_ECM_mot_den_index = " <<min_ECM_mot_den_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
	static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max_ECM_motility_density");
    if (max_ECM_mot_den_index < 0) std::exit(-1);
	static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal_ECM_motility_density");
    if (ideal_ECM_mot_den_index  < 0) std::exit(-1);
	
	// sample ECM - only changes for decoupling **should** be here as nothign gets written to the ECM...
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   

	double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 
	double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; 
	std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
	
	
	/****************************************Begin new migration direction update****************************************/

	// Select random direction (for random portion of motility vector) and begin building updated motility direction vector
	// (note - there is NO memeory of previous direction in this model - previous ECM-based motility used the current 
	// velocity vector to build off, not a random one - this could produce divergent behaviors between models)

	// See lab note book for more notes - MUST start with random vector. In the old method I defintiely used the previous motility vector in the method, but makes no sense here!

	// get random vector - cell's "intended" or chosen random direction
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_random = { cos(angle) , sin(angle) , 0.0 };

	// std::cout<<"D random "<<d_random<<std::endl;

	// get vector for chemotaxis (sample uE)
	std::vector<double> chemotaxis_grad = pCell->nearest_gradient(o2_index);

	// std::cout<<"D chemo"<<chemotaxis_grad<<std::endl;

	normalize( &chemotaxis_grad ); 

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)

	// New bias - bias such that the agents can more closely follow the gradient IF the written signals are stronger. 

	std::vector<double> d_motility;

	// d_motility = {0,0,0};
	
	if (parameters.ints("link_anisotropy_and_bias") == 0)
	{
		d_motility = (1-a) * d_random + a * chemotaxis_grad;
	}

	else if (parameters.ints("link_anisotropy_and_bias") == 1)
	{
		// NON-ECM linked way to signal 
		d_motility = (1-pCell->custom_data[chemotaxis_bias_index])*d_random + pCell->custom_data[chemotaxis_bias_index]*chemotaxis_grad;
	}

	else
	{
		std::cout<<"Must specify reader chemotaxis modeling mode - see XML parameter \"link_anisotropy_and_bias\" Halting!!!!!!"<<std::endl;
		abort();
		return;
	}

	normalize( &d_motility ); 


	// std::cout<<"D motility "<<d_motility<<std::endl;

	// to determine direction along f, find part of d_choice that is perpendicular to f; 
	std::vector<double> d_perp = d_motility - dot_product(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product( d_motility , d_perp ); 
	double c_2 = dot_product( d_motility, f ); 

	// std::cout<<"D_mot dot d_perp c_1 = "<<c_1<<std::endl;
	// std::cout<<"D_mot dot f c_2 = "<<c_2<<std::endl;

	// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy

	double gamma = pCell->custom_data[ECM_sensitivity_index] * a; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.
	// std::cout<<"anisotropy = "<<a<<std::endl;
	// std::cout<<"ECM sensitivity index = "<<pCell->custom_data[ECM_sensitivity_index]<<std::endl;
	// std::cout<<"gamma = "<< gamma <<std::endl;
	// std::cout<<"(1.0-gamma)*c_1*d_perp "<<(1.0-gamma)*c_1*d_perp<<std::endl;
	// std::cout<<"c_2*f"<<c_2*f<<std::endl;

	phenotype.motility.migration_bias_direction = (1.0-gamma)*c_1*d_perp + c_2*f;
	// std::cout<<"migration_bias_direction before normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
	if(parameters.bools("normalize_ECM_influenced_motility_vector") == true)
	{
		// normalize( &phenotype.motility.migration_bias_direction ); // only needed if not running through the update_migration_bias code/bias not set to 1.0
		// std::cout<<"migration_bias_direction after normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
		pCell->phenotype.motility.migration_speed = 1.0;
	}

	else
	{
		pCell->phenotype.motility.migration_speed = norm( phenotype.motility.migration_bias_direction);
		//  std::cout<<"Magnitutude of motility vector is "<< pCell->phenotype.motility.migration_speed<<std::endl;
	}
	
	
	phenotype.motility.migration_bias = 1.0; // MUST be set at 1.0 so that standard update_motility function doesn't add random motion. 

	// double magnitude = norm( phenotype.motility.motility_vector);	

	// std::cout<<"Magnitutude of motility vector is "<< magnitude<<std::endl;

	// if(magnitude > 0.00000001)
	// {
	// 	std::cout<<"Cell is moving!!!!"<<std::endl;
	// }

	/****************************************END new migration direction update****************************************/


	/*********************************************Begin speed update***************************************************/
	
	// New speed update (06.18.19) - piece wise continous
	
	double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];

	// std::cout<<"ECM_density = "<<ECM_density<<std::endl;

	if (ECM_density <= rho_low)
	{
		pCell->phenotype.motility.migration_speed = 0.0;

	}

	else if (rho_low < ECM_density && ECM_density <= rho_ideal)
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
		// std::cout<<"max_cell_speed = "<<pCell->custom_data[max_cell_speed_index]<<std::endl;
		// std::cout<<"speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	}

	else //if (ECM_density >= rho_high)
	{
		pCell->phenotype.motility.migration_speed = 0.0;
	}

	// int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	// int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	// int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	// int necrosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 

	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// std::cout<<"cell adhesion = "<<pCell->phenotype.mechanics.cell_cell_adhesion_strength <<std::endl;
	// std::cout<<"cell repulsion = "<<pCell->phenotype.mechanics.cell_cell_repulsion_strength <<std::endl;
	// std::cout<<"cell persistence time ="<<pCell->phenotype.motility.persistence_time <<std::endl;
	// std::cout<<"cell transition rates = "<<phenotype.death.rates[apoptosis_index] <<std::endl;
	// std::cout<<"cell death rates = "<<phenotype.death.rates[necrosis_index] <<std::endl;

    // fibroblast.phenotype.death.rates[apoptosis_index] = 0.0;
	// fibroblast.phenotype.death.rates[necrosis_index] = 0.0;
	
	// END New speed update 

	// Old parabolic update  - just the one line!

	// pCell->phenotype.motility.migration_speed = (-(4.0)*pow((ECM_density-0.5),2.0) + 1.0) * pCell->custom_data[max_cell_speed_index];
	//   std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;


	/*********************************************END speed update***************************************************/

	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
	return; 
}

/* To eliminate memory, set b_h to zero. To eliminate ECM influence, set a to 0 (permanently) or ECM senstiivity to zero,
   This model uses memory from previous motility vector AND previous anisotropy location. */

void ECM_informed_motility_update_model_w_memory ( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	if(phenotype.death.dead == true)
	{
		
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		std::cout<<2<<std::endl;
		std::cout<<"Cell is dead"<<std::endl;
	}
	// Updates cell bias vector and cell speed based on the ECM density, anisotropy, and fiber direction
	
	// find location of variables and base parameter values
	static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 

	static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max speed" ); 
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis bias");
	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM sensitivity");
	static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min ECM motility density");
	static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max ECM motility density");
	static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal ECM motility density");
	static int base_motility_hysteresis_bias = pCell->custom_data.find_variable_index( "Base hysteresis bias");
	static int previous_anistropy = pCell->custom_data.find_variable_index( "previous anisotropy" );

	// sample ECM - only changes for decoupling **should** be here as nothign gets written to the ECM...
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   

	double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 
	double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; 
	std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
	
	
	/****************************************Begin new migration direction update****************************************/

	// Select random direction (for random portion of motility vector) and begin building updated motility direction vector
	// (note - there is NO memeory of previous direction in this model - previous ECM-based motility used the current 
	// velocity vector to build off, not a random one - this could produce divergent behaviors between models)

	// See lab note book for more notes - MUST start with random vector. In the old method I defintiely used the previous motility vector in the method, but makes no sense here!

	// get random vector - cell's "intended" or chosen random direction
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_random = { cos(angle) , sin(angle) , 0.0 };

	// std::cout<<"D random "<<d_random<<std::endl;

	// get vector for chemotaxis (sample uE)
	std::vector<double> chemotaxis_grad = pCell->nearest_gradient(o2_index);

	// std::cout<<"D chemo"<<chemotaxis_grad<<std::endl;

	normalize( &chemotaxis_grad ); 

	// get vector for hysteris (get previous motility vector)

	// hysterises bias (should put in custom data ***)
	// double b_h = 1.0;

	std::vector<double> previous_motility_direction = phenotype.motility.migration_bias_direction; // double check this ***

	// std::vector<double> d_random_plus_previous = b_h * previous_motility_direction + pCell->custom_data[chemotaxis_bias_index] * chemotaxis_grad;

	// I am not sure how to still incorporatea chemotaxis - ask Paul. For testing, I am moving on and just replace chemotaxis with the other. 

	// Scale the cells sensitivity to previous direction by the amount of previous direction (a) there is available to remember!!!!
	double delta = pCell->custom_data[base_motility_hysteresis_bias]* pCell->custom_data[previous_anistropy];

	// std::cout<<"Previous a = "<< pCell->custom_data[previous_anistropy]<<std::endl;
	// std::cout<<"Current a = "<< a<<std::endl;

	// update memory to the newest a for next motiltiy update
	pCell->custom_data[previous_anistropy] = a;

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)
	std::vector<double> d_motility = (1-delta)*d_random + delta*previous_motility_direction;
	normalize( &d_motility ); 

	// std::cout<<"D motility "<<d_motility<<std::endl;

	// to determine direction along f, find part of d_choice that is perpendicular to f; 
	std::vector<double> d_perp = d_motility - dot_product(d_motility,f)*f; 
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product( d_motility , d_perp ); 
	double c_2 = dot_product( d_motility, f ); 

	// std::cout<<"D_mot dot d_perp c_1 = "<<c_1<<std::endl;
	// std::cout<<"D_mot dot f c_2 = "<<c_2<<std::endl;

	// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy

	double gamma = pCell->custom_data[ECM_sensitivity_index] * a; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.
	// std::cout<<"anisotropy = "<<a<<std::endl;
	// std::cout<<"ECM sensitivity index = "<<pCell->custom_data[ECM_sensitivity_index]<<std::endl;
	// std::cout<<"gamma = "<< gamma <<std::endl;
	// std::cout<<"(1.0-gamma)*c_1*d_perp "<<(1.0-gamma)*c_1*d_perp<<std::endl;
	// std::cout<<"c_2*f"<<c_2*f<<std::endl;

	phenotype.motility.migration_bias_direction = (1.0-gamma)*c_1*d_perp + c_2*f;
	// std::cout<<"migration_bias_direction before normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
	if(parameters.bools("normalize_ECM_influenced_motility_vector") == true)
	{
		// normalize( &phenotype.motility.migration_bias_direction ); // only needed if not running through the update_migration_bias code/bias not set to 1.0
		// std::cout<<"migration_bias_direction after normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
		pCell->phenotype.motility.migration_speed = 1.0;
	}

	else
	{
		pCell->phenotype.motility.migration_speed = norm( phenotype.motility.migration_bias_direction);
		//  std::cout<<"Magnitutude of motility vector is "<< pCell->phenotype.motility.migration_speed<<std::endl;
	}
	
	
	phenotype.motility.migration_bias = 1.0; // MUST be set at 1.0 so that standard update_motility function doesn't add random motion. 

	// double magnitude = norm( phenotype.motility.motility_vector);	

	// std::cout<<"Magnitutude of motility vector is "<< magnitude<<std::endl;

	// if(magnitude > 0.00000001)
	// {
	// 	std::cout<<"Cell is moving!!!!"<<std::endl;
	// }

	/****************************************END new migration direction update****************************************/


	/*********************************************Begin speed update***************************************************/
	
	// New speed update (06.18.19) - piece wise continous
	
	double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];

	// std::cout<<"ECM_density = "<<ECM_density<<std::endl;

	if (ECM_density <= rho_low)
	{
		pCell->phenotype.motility.migration_speed = 0.0;

	}

	else if (rho_low < ECM_density && ECM_density <= rho_ideal)
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	}

	else //if (ECM_density >= rho_high)
	{
		pCell->phenotype.motility.migration_speed = 0.0;
	}

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	int necrosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 

	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// std::cout<<"cell adhesion = "<<pCell->phenotype.mechanics.cell_cell_adhesion_strength <<std::endl;
	// std::cout<<"cell repulsion = "<<pCell->phenotype.mechanics.cell_cell_repulsion_strength <<std::endl;
	// std::cout<<"cell persistence time ="<<pCell->phenotype.motility.persistence_time <<std::endl;
	// std::cout<<"cell transition rates = "<<phenotype.death.rates[apoptosis_index] <<std::endl;
	// std::cout<<"cell death rates = "<<phenotype.death.rates[necrosis_index] <<std::endl;

    // fibroblast.phenotype.death.rates[apoptosis_index] = 0.0;
	// fibroblast.phenotype.death.rates[necrosis_index] = 0.0;
	
	// END New speed update 

	// Old parabolic update  - just the one line!

	// pCell->phenotype.motility.migration_speed = (-(4.0)*pow((ECM_density-0.5),2.0) + 1.0) * pCell->custom_data[max_cell_speed_index];
	//   std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;


	/*********************************************END speed update***************************************************/

	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
	return; 
}

// Sets speed based on voxel anisotropy (instead of magnitude of motility vector). With chemotaxis senstivity set to 1, cells will chemotax, but have speed vary with only anisotropy. Direction WILL NOT be effeced by ECM, only speed. 
// Only fiber alignemnt (self alignemnt) impacts the motility. 

void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	
	if(phenotype.death.dead == true)
	{
		
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		std::cout<<2<<std::endl;
		std::cout<<"Cell is dead"<<std::endl;
	}
	// Updates cell bias vector and cell speed based on the ECM density, anisotropy, and fiber direction
	
	// find location of variables and base parameter values
	// static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	// static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	
	static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max speed" ); 
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis bias");
	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM sensitivity");
	static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min ECM motility density");
	static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max ECM motility density");
	static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal ECM motility density");
	
	// sample ECM - only changes for decoupling **should** be here as nothign gets written to the ECM...
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   

	double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 
	double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; 
	std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
	
	
	/****************************************Begin new migration direction update****************************************/

	// Select random direction (for random portion of motility vector) and begin building updated motility direction vector
	// (note - there is NO memeory of previous direction in this model - previous ECM-based motility used the current 
	// velocity vector to build off, not a random one - this could produce divergent behaviors between models)

	// See lab note book for more notes - MUST start with random vector. In the old method I defintiely used the previous motility vector in the method, but makes no sense here!

	// get random vector - cell's "intended" or chosen random direction
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_random = { cos(angle) , sin(angle) , 0.0 };

	// std::cout<<"D random "<<d_random<<std::endl;

	// get vector for chemotaxis (sample uE)
	std::vector<double> chemotaxis_grad = pCell->nearest_gradient(o2_index);

	// std::cout<<"D chemo"<<chemotaxis_grad<<std::endl;

	normalize( &chemotaxis_grad ); 

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)
	std::vector<double> d_motility = (1-pCell->custom_data[chemotaxis_bias_index])*d_random + pCell->custom_data[chemotaxis_bias_index]*chemotaxis_grad;
	normalize( &d_motility ); 


	// std::cout<<"D motility "<<d_motility<<std::endl;

	// to determine direction along f, find part of d_choice that is perpendicular to f; 
	std::vector<double> d_perp = d_motility - dot_product(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product( d_motility , d_perp ); 
	double c_2 = dot_product( d_motility, f ); 

	// std::cout<<"D_mot dot d_perp c_1 = "<<c_1<<std::endl;
	// std::cout<<"D_mot dot f c_2 = "<<c_2<<std::endl;

	// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy

	double gamma = pCell->custom_data[ECM_sensitivity_index] * a; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.
	// std::cout<<"anisotropy = "<<a<<std::endl;
	// std::cout<<"ECM sensitivity index = "<<pCell->custom_data[ECM_sensitivity_index]<<std::endl;
	// std::cout<<"gamma = "<< gamma <<std::endl;
	// std::cout<<"(1.0-gamma)*c_1*d_perp "<<(1.0-gamma)*c_1*d_perp<<std::endl;
	// std::cout<<"c_2*f"<<c_2*f<<std::endl;

	phenotype.motility.migration_bias_direction = (1.0-gamma)*c_1*d_perp + c_2*f;
	// std::cout<<"migration_bias_direction before normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
	if(parameters.bools("normalize_ECM_influenced_motility_vector") == true)
	{
		// normalize( &phenotype.motility.migration_bias_direction ); // only needed if not running through the update_migration_bias code/bias not set to 1.0
		// std::cout<<"migration_bias_direction after normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
		pCell->phenotype.motility.migration_speed = 1.0;
	}

	else
	{
		pCell->phenotype.motility.migration_speed = a*norm( phenotype.motility.migration_bias_direction);
		//  std::cout<<"Magnitutude of motility vector is "<< pCell->phenotype.motility.migration_speed<<std::endl;
	}
	
	
	phenotype.motility.migration_bias = 1.0; // MUST be set at 1.0 so that standard update_motility function doesn't add random motion. 

	// double magnitude = norm( phenotype.motility.motility_vector);	

	// std::cout<<"Magnitutude of motility vector is "<< magnitude<<std::endl;

	// if(magnitude > 0.00000001)
	// {
	// 	std::cout<<"Cell is moving!!!!"<<std::endl;
	// }

	/****************************************END new migration direction update****************************************/


	/*********************************************Begin speed update***************************************************/
	
	// New speed update (06.18.19) - piece wise continous
	
	double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];

	// std::cout<<"ECM_density = "<<ECM_density<<std::endl;

	if (ECM_density <= rho_low)
	{
		pCell->phenotype.motility.migration_speed = 0.0;

	}

	else if (rho_low < ECM_density && ECM_density <= rho_ideal)
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	}

	else //if (ECM_density >= rho_high)
	{
		pCell->phenotype.motility.migration_speed = 0.0;
	}

	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	int necrosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 

	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// std::cout<<"cell adhesion = "<<pCell->phenotype.mechanics.cell_cell_adhesion_strength <<std::endl;
	// std::cout<<"cell repulsion = "<<pCell->phenotype.mechanics.cell_cell_repulsion_strength <<std::endl;
	// std::cout<<"cell persistence time ="<<pCell->phenotype.motility.persistence_time <<std::endl;
	// std::cout<<"cell transition rates = "<<phenotype.death.rates[apoptosis_index] <<std::endl;
	// std::cout<<"cell death rates = "<<phenotype.death.rates[necrosis_index] <<std::endl;

    // fibroblast.phenotype.death.rates[apoptosis_index] = 0.0;
	// fibroblast.phenotype.death.rates[necrosis_index] = 0.0;
	
	// END New speed update 

	// Old parabolic update  - just the one line!

	// pCell->phenotype.motility.migration_speed = (-(4.0)*pow((ECM_density-0.5),2.0) + 1.0) * pCell->custom_data[max_cell_speed_index];
	//   std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;


	/*********************************************END speed update***************************************************/

	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
	return; 
}

void rightward_deterministic_cell_march (Cell* pCell , Phenotype& phenotype , double dt )
{

	pCell->phenotype.motility.migration_bias_direction[0] = 1.0;
   	pCell->phenotype.motility.migration_bias_direction[1] = 0.0;
    pCell->phenotype.motility.migration_bias_direction[2] = 0.0;

	// std::cout<<"Am I running?"<<std::endl;

	return;
}

void reset_cell_position( void ) // for cell mark/ECM change test
{
	int n = default_microenvironment_options.X_range[0] + 10.0;

	std::cout<<1<<std::endl;

	for(int i = 0; i < (*all_cells).size(); i++)
	{
		Cell* pCell = (*all_cells)[i];
		pCell->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
		n = n + 10.0;
	}
		
}

void chemotaxis_oxygen( Cell* pCell , Phenotype& phenotype , double dt )
{
	/*********************************************Chemotaxsis update***************************************************/
	
	// sample uE
	static int inflam_sig_index = microenvironment.find_density_index( "inflammatory_signal" ); 
	
	phenotype.motility.is_motile = true; 
	phenotype.motility.migration_bias = parameters.doubles("inflam_sig_migration_bias_for_fibroblasts"); //0.95;
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(inflam_sig_index);

	normalize( &( phenotype.motility.migration_bias_direction ) );

	/*********************************************Begin speed update***************************************************/

	// sample ECM - only changes for decoupling **should** be here as nothign gets written to the ECM...
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 

	static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max speed" ); 
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis bias");
	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM sensitivity");
	static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min ECM motility density");
	static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max ECM motility density");
	static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal ECM motility density");
	
	// New speed update (06.18.19) - piece wise continous
	
	double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];

	// std::cout<<"ECM_density = "<<ECM_density<<std::endl;
	pCell->phenotype.motility.migration_speed = 1.0;
	if (ECM_density <= rho_low)
	{
		pCell->phenotype.motility.migration_speed = 0.0;

	}

	else if (rho_low < ECM_density && ECM_density <= rho_ideal)
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))
		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
		// std::cout<<"max_cell_speed = "<<pCell->custom_data[max_cell_speed_index]<<std::endl;
		// std::cout<<"rho_ideal = "<<rho_ideal<<std::endl;
		// std::cout<<"rho_low = "<<rho_low<<std::endl;
		// std::cout<<"rho_high = "<<rho_high<<std::endl;
		// std::cout<<"speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
		// pCell->phenotype.motility.migration_speed = 10.0;
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	}

	else //if (ECM_density >= rho_high)
	{
		pCell->phenotype.motility.migration_speed = 0.0;
	}

	if(phenotype.death.dead == true)
	{
		phenotype.motility.is_motile=false;
		pCell->functions.update_phenotype = NULL;
		pCell->functions.update_migration_bias = NULL;
		std::cout<<"Cell is dead"<<std::endl;


	}	
	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
   	// std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;
	
	return; 
}

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell , Phenotype& phenotype , double dt ) 
{	
	// o2-based birth and death, nothing more 
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	return; 
}


// void cancer_cell_phenotype_model0( Cell* pCell , Phenotype& phenotype , double dt )
// {
//     // o2-based birth and death, nothing more
//     tumor_cell_phenotype_with_oncoprotein(pCell,phenotype,dt);

//     return;
// }

/* // not being used currently (03.11.19)

std::vector<std::string> ECM_anisotropy_coloring_function( Cell* pCell)
{
    std::vector< std::string > output( 4, "black" );
    char color [1024];
    int ecm_index = pCell->get_current_voxel_index();
   	
    double anisotropy = ecm.ecm_data[ecm_index].anisotropy;
    sprintf( color , "rgb(%d,%d,%d)" , int(anisotropy*255.0), int(anisotropy* 255.0), int(255-anisotropy*255));
   	// std::cout<<color<<std::endl;
   	// rgb(0,0,255)
    output[0] = color;
    output[2] = color;
    return output;
} */

std::vector<std::string> AMIGOS_invasion_coloring_function( Cell* pCell )
{
	// fibroblasts are blue 
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		output[0] = "blue"; 
		output[2] = "blue"; 
		return output; 
	} 

	// cancer cells are yellow except for the marker cells when running testing mode (then 20 % of cells are red)
    
	if( pCell->type == 2 )
    {
        output[2] = "yellow";	
		output[0] = "yellow";
		
		// Return yellow for cancer cells and exit statement

		if(parameters.ints("unit_test_setup")==0)
		{return output;}

		// Return red for 20% of cancer cells if unit test is called for	

		if( pCell->ID % 5 == 0)
		{
			output[0] = "red";
        	output[2] = "red";
        	return output;
		}

		// If doing unit testing AND cell not selected as marker cell, return blue. 
		output[2] = "blue";	
		output[0] = "blue";


        return output;
    }
	
	// dead are red
	if( pCell->phenotype.death.dead == false )
	{
		output[0] = "red";
		output[2] = "red";
		
		return output; 
	}

	// if not, dead colors 
	
	if (pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::apoptotic )  // Apoptotic - Red
	{
		output[0] = "rgb(255,0,0)";
		output[2] = "rgb(125,0,0)";
	}
	
	// Necrotic - Brown
	if( pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_swelling || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic_lysed || 
		pCell->phenotype.cycle.current_phase().code == PhysiCell_constants::necrotic )
	{
		output[0] = "rgb(250,138,38)";
		output[2] = "rgb(139,69,19)";
	}	
	
	return output; 
}

// void no_birth_or_death_phenotype_model (Cell* pCell , Phenotype& phenotype , double dt) - well this is a smart idea to have for testing ... 
// {
// 	return;
// }

void fibroblast_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int hypoxic_i = pCell->custom_data.find_variable_index( "hypoxic switch value" ); 	
	static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" ); 
	
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	double pO2 = (pCell->nearest_density_vector())[oxygen_i]; // PhysiCell_constants::oxygen_index]; 
	
	// set death and birth 
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt); 

	// ALWAYS MOTILE

    // phenotype.motility.is_motile = true;
    
	/* if( pO2 > pCell->custom_data[hypoxic_i] )
	{
		// proliferate (don't overwrite)
	    // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) =
	    // 10.0 * pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
		// turn off motility
		phenotype.motility.is_motile = true ;
	}
	else
	{
		// don't proliferate,
	    // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
		// turn on motility
		phenotype.motility.is_motile = true;
	} */
	return;
}


void cancer_cell_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt )
{
    static int hypoxic_i = pCell->custom_data.find_variable_index( "hypoxic switch value" );
    static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" );
    
    
    int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
    int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
    
    double pO2 = (pCell->nearest_density_vector())[oxygen_i]; // PhysiCell_constants::oxygen_index];
    
    // set death and birth
    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
    
    // ALWAYS MOTILE
    // phenotype.motility.is_motile = true;


	if(phenotype.death.dead == true)
	{
		phenotype.motility.is_motile=false;
		pCell->functions.update_phenotype = NULL;
		pCell->functions.update_migration_bias = NULL;
		std::cout<<"cancer cell is dead"<<std::endl;
		// std::cout<<"Is it motile?"<<std::endl;
		// if(phenotype.motility.is_motile==true)
		// {
		// 	std::cout<<"Yes"<<std::endl;
		// 	phenotype.motility.is_motile==false;
		// 	pCell->functions.update_phenotype = NULL;
		// 	pCell->functions.update_migration_bias = NULL;
		// }

		// if(phenotype.motility.is_motile==false)
		// {
		// 	std::cout<<"No"<<std::endl;
		// }

	}

    
   	// phenotype.motility.migration_bias = 1.0;
    
	/* if( pO2 > pCell->custom_data[hypoxic_i] )
	{
		// proliferate (don't overwrite)
       	// phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) =
        // 10.0 * pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
		// turn off motility
		phenotype.motility.is_motile = true ;
	}
	else
	{
		// don't proliferate,
        // phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
		// turn on motility
		phenotype.motility.is_motile = true;
	} */
    return;
}

long fibonacci(unsigned n) // just being used for timing. 
{
    if (n < 2) return n;
    return fibonacci(n-1) + fibonacci(n-2);
}

// uses cell motility vector for fiber reorientation

void ecm_update_from_cell_motility_vector(Cell* pCell , Phenotype& phenotype , double dt) 
{

	// Find correct items
	std::vector<double> cell_position = pCell->position;
	// std::cout<<cell_position<<std::endl;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	// std::cout<<nearest_ecm_voxel_index<<std::endl;
	// std::cin.get();
    
	static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target_ECM_density");
	static int Cell_ECM_production_rate_index = pCell->custom_data.find_variable_index( "ECM_production_rate");
	static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy_increase_rate");
	static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber_realignment_rate");
    

    // Cell-ECM density interaction

    double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;
    double r = pCell->custom_data[Cell_ECM_production_rate_index];
	// std::cout<<"ECM_production_rate: "<<r<<std::endl;
    
    ecm.ecm_voxels[nearest_ecm_voxel_index].density = ECM_density + r * dt  * (pCell->custom_data[Cell_ECM_target_density_index] - ECM_density);
    // END Cell-ECM density interaction

	// Cell-ECM Fiber realingment - continous then discrete

	if( parameters.ints("discrete_ECM_remodeling") == 1)
	{

		// Get index for accessing the ECM_fiber_alignment data structure and then copy the correct value
		// int n = pCell->get_current_voxel_index();
		std::vector<double> ECM_orientation = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment; 

		double anisotropy = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy;
		double migration_speed = pCell->phenotype.motility.migration_speed;
		double r_0 = pCell->custom_data[Cell_fiber_realignment_rate_index]*migration_speed; // 1/10.0 // min-1 // NOTE!!! on 08.06.18 run - this wasn't multiplied by migration_speed!!! should be the same but worth noting!!!!
		// std::cout<<r_0<<std::endl;
		double r_realignment = r_0 * (1-anisotropy);

		double ddotf;
		std::vector<double> norm_cell_motility = phenotype.motility.motility_vector;
		// norm_cell_motility.resize(3,0.0);
		// norm_cell_motility = phenotype.motility.motility_vector;
		normalize(&norm_cell_motility);

		ddotf = dot_product(ECM_orientation, norm_cell_motility);

		ECM_orientation = sign_function(ddotf) * ECM_orientation; // flips the orientation vector so that it is aligned correctly with the moving cell for proper reoirentation later.
		
		std::vector<double> f_minus_d;
		f_minus_d.resize(3,0.0);

		// f_minus_d = ECM_orientation - norm_cell_motility; // Fix this later

		for(int i = 0; i < 3; i++)
		{
			if (ddotf<0.0)
			{
				ECM_orientation = -1.0 * ECM_orientation;
			}
			f_minus_d[i] = ECM_orientation[i] - norm_cell_motility[i]; // 06.05.19 - fixed 
			ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment[i] -= dt * r_realignment * f_minus_d[i]; 
		}
		
		
		normalize(&(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment)); // why by reference??

		// End Cell-ECM Fiber realingment
	
		// Cell-ECM Anisotrophy Modification
		
		double r_a0 = pCell->custom_data[Cell_anistoropy_rate_of_increase_index] ; // min-1
		
		double r_anisotropy = r_a0 * migration_speed;
		
		ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy = anisotropy + r_anisotropy * dt  * (1- anisotropy);
		
		// END Cell-ECM Anisotropy Modification


	}

	else if (parameters.ints("discrete_ECM_remodeling") == 0)
	{
		if (ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy == 1)
		{
			// std::cout<<"pass code"<<std::endl;
		}
		
		else
		{
		if (norm(phenotype.motility.motility_vector) == 0)
		{std::cout<<"Motility vector norm = 0"<<std::endl;}

 		ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment = phenotype.motility.motility_vector;
		
		if (norm(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment) == 0)
		{std::cout<<"ECM orientation vector norm = 0"<<std::endl;}

		ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy = 1;
		}
	}

	else
	{
		std::cout<<"Must specify ECM remodeling mode - see XML parameter \"discrete_ECM_remodeling\" Halting!!!!!!"<<std::endl;
		abort();
		return;
	}
	
	
    
    return;

}


// uses cell velocity vector for fiber reorientation

void ecm_update_from_cell_velocity_vector(Cell* pCell , Phenotype& phenotype , double dt)
{

	// Find correct fields
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target ECM density");
	static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy increase rate");
	static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber realignment rate");
    
    // Cell-ECM density interaction
    double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;
    double r = 1.0;
    
    ecm.ecm_voxels[nearest_ecm_voxel_index].density = ECM_density + r * dt  * (pCell->custom_data[Cell_ECM_target_density_index] - ECM_density);
    
    // END Cell-ECM density interaction
    
    // Cell-ECM Fiber realingment

	// Get index for accessing the ECM_fiber_alignment data structure and then copy the correct value
	// int n = pCell->get_current_voxel_index();
	std::vector<double> ECM_orientation = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment; 

	double anisotropy = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy;
    double migration_speed = pCell->phenotype.motility.migration_speed;
	
    double r_0 = pCell->custom_data[Cell_fiber_realignment_rate_index]*migration_speed; // 1/10.0 // min-1 // NOTE!!! on 08.06.18 run - this wasn't multiplied by migration_speed!!! should be the same but worth noting!!!!
	// std::cout<<r_0<<std::endl;
    double r_realignment = r_0 * (1-anisotropy);

	double ddotf;
	std::vector<double> norm_cell_velocity = pCell->velocity;
	// norm_cell_velocity.resize(3,0.0);
	// norm_cell_velocity = pCell->velocity;
	normalize(&norm_cell_velocity);

	ddotf = dot_product(ECM_orientation, norm_cell_velocity);
	
	ECM_orientation = sign_function(ddotf) * ECM_orientation; // flips the orientation vector so that it is aligned correctly with the moving cell for proper reoirentation later.
	
	std::vector<double> f_minus_d;
	f_minus_d.resize(3,0.0);

	// f_minus_d = ECM_orientation - norm_cell_motility; // Fix this later

	for(int i = 0; i < 3; i++)
	{
		f_minus_d[i] = ECM_orientation[i] - norm_cell_velocity[i]; // 06.05.19 - fixed 
		ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment[i] -= dt * r_realignment * f_minus_d[i]; 
	}
	
	
    normalize(&(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment)); // why by reference??


    // End Cell-ECM Fiber realingment
    
    // Cell-ECM Anisotrophy Modification
    
    double r_a0 = pCell->custom_data[Cell_anistoropy_rate_of_increase_index] ; // min-1
	
    double r_anisotropy = r_a0 * migration_speed;
   	
    ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy = anisotropy + r_anisotropy * dt  * (1- anisotropy);
    
    // END Cell-ECM Anisotropy Modification
    
    return;

}

void SVG_plot_custom( std::string filename , Microenvironment& M, double z_slice , double time, std::vector<std::string> (*cell_coloring_function)(Cell*), std::string line_pattern )
{
	double X_lower = M.mesh.bounding_box[0];
	double X_upper = M.mesh.bounding_box[3];
 
	double Y_lower = M.mesh.bounding_box[1]; 
	double Y_upper = M.mesh.bounding_box[4]; 

	double plot_width = X_upper - X_lower; 
	double plot_height = Y_upper - Y_lower; 

	double font_size = 0.025 * plot_height; // PhysiCell_SVG_options.font_size; 
	double top_margin = font_size*(.2+1+.2+.9+.5 ); 

	// open the file, write a basic "header"
	std::ofstream os( filename , std::ios::out );
	if( os.fail() )
	{ 
		std::cout << std::endl << "Error: Failed to open " << filename << " for SVG writing." << std::endl << std::endl; 

		std::cout << std::endl << "Error: We're not writing data like we expect. " << std::endl
		<< "Check to make sure your save directory exists. " << std::endl << std::endl
		<< "I'm going to exit with a crash code of -1 now until " << std::endl 
		<< "you fix your directory. Sorry!" << std::endl << std::endl; 
		exit(-1); 
	} 
	
	Write_SVG_start( os, plot_width , plot_height + top_margin );

	// draw the background 
	Write_SVG_rect( os , 0 , 0 , plot_width, plot_height + top_margin , 0.002 * plot_height , "white", "white" );
	
	// bool Write_SVG_circle( std::ostream& os, double center_x, double center_y, double radius, double stroke_size, 
    //                    std::string stroke_color , std::string fill_color )

	// write the simulation time to the top of the plot
 
	char* szString; 
	szString = new char [1024]; 
 
	int total_cell_count = all_cells->size(); 
 
	double temp_time = time; 

	std::string time_label = formatted_minutes_to_DDHHMM( temp_time ); 
 
	sprintf( szString , "Current time: %s, z = %3.2f %s", time_label.c_str(), 
		z_slice , PhysiCell_SVG_options.simulation_space_units.c_str() ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1), 
		font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	sprintf( szString , "%u agents" , total_cell_count ); 
	Write_SVG_text( os, szString, font_size*0.5,  font_size*(.2+1+.2+.9), 
		0.95*font_size, PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	
	delete [] szString; 

	
	// add an outer "g" for coordinate transforms 
	
	os << " <g id=\"tissue\" " << std::endl 
	   << "    transform=\"translate(0," << plot_height+top_margin << ") scale(1,-1)\">" << std::endl; 
	   
	// prepare to do mesh-based plot (later)
	
	double dx_stroma = M.mesh.dx; 
	double dy_stroma = M.mesh.dy; 
	
	os << "  <g id=\"ECM\">" << std::endl; 
  
	int ratio = 1; 
	double voxel_size = dx_stroma / (double) ratio ; 
  
	double half_voxel_size = voxel_size / 2.0; 
	double normalizer = 78.539816339744831 / (voxel_size*voxel_size*voxel_size); 
 
 // color in the background ECM

	os << "  </g>" << std::endl; 
 
	// plot intersecting cells 
	os << "  <g id=\"cells\">" << std::endl; 
	for( int i=0 ; i < total_cell_count ; i++ )
	{
		Cell* pC = (*all_cells)[i]; // global_cell_list[i]; 
  
		static std::vector<std::string> Colors; 
		if( fabs( (pC->position)[2] - z_slice ) < pC->phenotype.geometry.radius )
		{
			double r = pC->phenotype.geometry.radius ; 
			double rn = pC->phenotype.geometry.nuclear_radius ; 
			double z = fabs( (pC->position)[2] - z_slice) ; 
   
			Colors = cell_coloring_function( pC ); 

			os << "   <g id=\"cell" << pC->ID << "\">" << std::endl; 
  
			// figure out how much of the cell intersects with z = 0 
   
			double plot_radius = sqrt( r*r - z*z ); 

			Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
				plot_radius , 0.5, Colors[1], Colors[0] ); 

			// plot the nucleus if it, too intersects z = 0;
			if( fabs(z) < rn && PhysiCell_SVG_options.plot_nuclei == true )
			{   
				plot_radius = sqrt( rn*rn - z*z ); 
			 	Write_SVG_circle( os, (pC->position)[0]-X_lower, (pC->position)[1]-Y_lower, 
					plot_radius, 0.5, Colors[3],Colors[2]); 
			}					  
			os << "   </g>" << std::endl;
		}
	}
	os << "  </g>" << std::endl; 

	// Plot guidelines the circle guides
	
	if (line_pattern == "concentric circles")
	{
		for (int i=0; i<plot_width/(2*parameters.doubles("ECM_dx")); i++)
		{
			double radius = plot_width/2 - (i *parameters.doubles("ECM_dx"));
			// std::cout<<"Index "<<i<<std::endl;
			// std::cout<<"Radius "<<radius<<std::endl;
			Write_SVG_circle( os, plot_width/2, plot_height/2, radius, 0.5, "black" , "none" );

		}
	}

	else if (line_pattern == "vertical lines")
	{
		for (int i=0; i<plot_width/parameters.doubles("ECM_dx"); i++)
		{
			double x_line_position = parameters.doubles("ECM_dx")*i;
			// std::cout<<"Index "<<i<<std::endl;
			// std::cout<<"X position "<<x_line_position<<std::endl;			
			Write_SVG_line(os, x_line_position, 0, x_line_position, plot_height, 0.5, "black");
			// bool Write_SVG_line( std::ostream& os , double start_x, double start_y, double end_x , double end_y, double thickness, 
            //         std::string stroke_color )
		}
	}

	else if (line_pattern == "horizontal lines")
	{
		for (int i=0; i<plot_height/parameters.doubles("ECM_dy"); i++)
		{
			double y_line_position = parameters.doubles("ECM_dy")*i;
			// std::cout<<"Index "<<i<<std::endl;
			// std::cout<<"Y position "<<y_line_position<<std::endl;			
			Write_SVG_line(os, 0, y_line_position, plot_width, y_line_position, 0.5, "black");
			// bool Write_SVG_line( std::ostream& os , double start_x, double start_y, double end_x , double end_y, double thickness, 
            //         std::string stroke_color )
		}
	}

	else if (line_pattern == "none") {} // Don't make lines!!!

	else if (line_pattern != "none" || "horizontal lines" || "vertical lines" || "concentric circles")
	{
		std::cout<<"Use of this custom SVG output function requires specifying \"none\" or a line pattern" <<std::endl;
		std::cout<<"Halting: see inputs to custom SVG function \"SVG_plot_custom\"" << std::endl;
		abort();
		return;

	}

	
	// end of the <g ID="tissue">
	os << " </g>" << std::endl; 
 
	// draw a scale bar
 
	double bar_margin = 0.025 * plot_height; 
	double bar_height = 0.01 * plot_height; 
	double bar_width = PhysiCell_SVG_options.length_bar; 
	double bar_stroke_width = 0.001 * plot_height; 
	
	std::string bar_units = PhysiCell_SVG_options.simulation_space_units; 
	// convert from micron to mm
	double temp = bar_width;  

	if( temp > 999 && std::strstr( bar_units.c_str() , PhysiCell_SVG_options.mu.c_str() )   )
	{
		temp /= 1000;
		bar_units = "mm";
	}
	// convert from mm to cm 
	if( temp > 9 && std::strcmp( bar_units.c_str() , "mm" ) == 0 )
	{
		temp /= 10; 
		bar_units = "cm";
	}
	
	szString = new char [1024];
	sprintf( szString , "%u %s" , (int) round( temp ) , bar_units.c_str() );
 
	Write_SVG_rect( os , plot_width - bar_margin - bar_width  , plot_height + top_margin - bar_margin - bar_height , 
		bar_width , bar_height , 0.002 * plot_height , "rgb(255,255,255)", "rgb(0,0,0)" );
	Write_SVG_text( os, szString , plot_width - bar_margin - bar_width + 0.25*font_size , 
		plot_height + top_margin - bar_margin - bar_height - 0.25*font_size , 
		font_size , PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() ); 
	
	delete [] szString; 

	// // plot runtime 
	// szString = new char [1024]; 
	// RUNTIME_TOC(); 
	// std::string formatted_stopwatch_value = format_stopwatch_value( runtime_stopwatch_value() );
	// Write_SVG_text( os, formatted_stopwatch_value.c_str() , bar_margin , top_margin + plot_height - bar_margin , 0.75 * font_size , 
	// 	PhysiCell_SVG_options.font_color.c_str() , PhysiCell_SVG_options.font.c_str() );
	// delete [] szString; 

	// draw a box around the plot window
	Write_SVG_rect( os , 0 , top_margin, plot_width, plot_height , 0.002 * plot_height , "rgb(0,0,0)", "none" );
	
	// close the svg tag, close the file
	Write_SVG_end( os ); 
	os.close();
 
	return; 
}

void write_ECM_Data_matlab( std::string filename )

{

    int number_of_data_entries = ecm.ecm_mesh.voxels.size();

    int size_of_each_datum = 8;

	
	// static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	// static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 

    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "ECM_Data" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )

    {

	    fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp ); // 1

        fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp ); // 2

        fwrite( (char*) &( ecm.ecm_mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp ); //3
		
		fwrite( (char*) &( ecm.ecm_voxels[i].anisotropy), sizeof(double) , 1 , fp ); // 4
	
        fwrite( (char*) &( ecm.ecm_voxels[i].density), sizeof(double) , 1 , fp ); // 5

        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[0]), sizeof(double) , 1 , fp ); // 6

        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[1]), sizeof(double) , 1 , fp ); // 7

        fwrite( (char*) &( ecm.ecm_voxels[i].ecm_fiber_alignment[2]), sizeof(double) , 1 , fp ); // 8

		// This will only work if the diffusion and ECM meshes are the same size. Commenting out for actualrunning. To do a direct comparison, leave them in and change length. Will have to change the vizualization to get this from the regular BioFVM outputs.

		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][0]), sizeof(double) , 1 , fp ); // 9

		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][1]), sizeof(double) , 1 , fp ); // 10

		// fwrite( (char*) &( microenvironment.gradient_vector(i)[0][2]), sizeof(double) , 1 , fp ); // 11

    }



    fclose( fp );



    return;

}
