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

#include "./AMIGOS-invasion.h"
// #include "./ECM.h"
#include <chrono>  // for high_resolution_clock - https://www.pluralsight.com/blog/software-development/how-to-measure-execution-time-intervals-in-c--
Cell_Definition leader_cell; 
Cell_Definition follower_cell; 
std::vector< std::vector<double> > ECM_fiber_alignment; 
unsigned long long int counter=0; // counter for calculating average for the ad hoc timing I am doing ... 
int time_total = 0;

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	SeedRandom(0); 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// turn the default cycle model to live, 
	// so it's easier to turn off proliferation
	
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	// Make sure we're ready for 2D
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// use default proliferation and death 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	int necrosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 

	// For strict ECM invasion testing, why have on death and birth at all??

	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.0;
    cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0;
	cell_defaults.phenotype.death.rates[necrosis_index] = 0.0;

	cell_defaults.parameters.o2_proliferation_saturation = 38.0;
	cell_defaults.parameters.o2_reference = 38.0;
	
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = parameters.doubles("oxygen_uptake"); 
	std::cout<<cell_defaults.phenotype.secretion.uptake_rates[0]<<std::endl;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 38; 

	// Fields for leader-follower diffusing signal based communication. Commented out 03.11.19

	/*

	cell_defaults.phenotype.secretion.secretion_rates[1] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[1] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[1] = 1;

	cell_defaults.phenotype.secretion.secretion_rates[2] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[2] = 0;
	
		cell_defaults.phenotype.secretion.saturation_densities[2] = 1;

	*/

	// For phenotype switching if using
	
	// cell_defaults.functions.update_phenotype = switching_phenotype_model;
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// set default motility parameters (even for when off)
	
	cell_defaults.phenotype.motility.is_motile = true;
	// consider what "best" persistence time would be, given the voxel dimensions. 
	cell_defaults.phenotype.motility.persistence_time = parameters.doubles("default_persistence_time"); //10.0; // Voxels are 20 um in all dimensions. Given a top speed of 0.5 um/min, cells will likely be in one voxel for 10 minutes or more. So update of 10 isn't bad. Should consider "best" number later. 
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("default_cell_speed");
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	cell_defaults.phenotype.motility.migration_bias = 1.0;// completely random - setting in update_migration_bias - might wnat to call that immediately thing
	cell_defaults.functions.update_migration_bias = ECM_informed_motility_update; // NEW!
	// std::cout<<cell_defaults.functions.update_migration_bias<<std::endl;
	// std::cin.get();

	// add custom data 
	cell_defaults.custom_data.add_variable( "min ECM motility density", "dimensionless", parameters.doubles( "rho_L") );  // Minimum ECM density required for cell motility
	cell_defaults.custom_data.add_variable( "max ECM motility density", "dimensionless", parameters.doubles( "rho_H") );  // Maximum ECM density allowing cell motility
	cell_defaults.custom_data.add_variable( "ideal ECM motility density", "dimensionless", parameters.doubles( "rho_I") );  // Ideal ECM density cell motility
	cell_defaults.custom_data.add_variable( "max speed", "micron/min" , parameters.doubles( "default_cell_speed") ); // Maximum migration speed
	cell_defaults.custom_data.add_variable( "chemotaxis bias", "dimensionless", parameters.doubles( "default_chemotaxis_bias") ); 
	cell_defaults.custom_data.add_variable( "ECM sensitivity", "dimensionless", parameters.doubles("default_ECM_sensitivity") );
	cell_defaults.custom_data.add_variable( "hypoxic switch value" , "mmHg", 38 );
	cell_defaults.custom_data.add_variable( "target ECM density", "dimensionless", parameters.doubles( "default_ECM_density_target") ); 

			// <unit_test_setup description="Specifies cell parameters for consistent unit tests of ECM influenced mechanics and mechanics influence on ECM - sets adhesion to 1.25, repulsion to 25, and speed to 1.0" type="bool">cells at left boundary/march</unit_test_setup>

	if( parameters.ints("unit_test_setup") == 1)
	{
		cell_defaults.phenotype.motility.persistence_time = 10.0; 
		cell_defaults.phenotype.motility.migration_speed = 1.0;
		cell_defaults.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
		cell_defaults.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
		cell_defaults.phenotype.secretion.uptake_rates[0] = 0.0;
		cell_defaults.custom_data.add_variable( "max speed", "micron/min" , 1.0 ); // Maximum migration speed
		cell_defaults.custom_data.add_variable( "min ECM motility density", "dimensionless", 0.0 );  // Minimum ECM density required for cell motility
		cell_defaults.custom_data.add_variable( "max ECM motility density", "dimensionless", 1.0 );  // Maximum ECM density allowing cell motility
		cell_defaults.custom_data.add_variable( "ideal ECM motility density", "dimensionless", 0.5 );  // Ideal ECM density cell motility
		cell_defaults.custom_data.add_variable( "target ECM density", "dimensionless", 0.5 ); 

		std::cout<< "running unit test setup follower"<<std::endl;

	}

	else if ( parameters.ints("unit_test_setup") == 0)
	{
		std::cout<<"not in unit test mode"<<std::endl;
	}

	else
	{
		std::cout<<"WARNING!!!!! Cell parameters not set correctly - unit test set up must either be true or false!!!!"<<std::endl;
	}
	
	
	// leader cells 
	
	leader_cell = cell_defaults;
	leader_cell.name = "leader cell"; 
	leader_cell.type = 1; 

	// Temperarily eliminating leader/follower signal

	// Obviously missing - add later
    
	// 10% proliferation 
    // leader_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.10;

	/*******************************************For "march" simulation****************************************/

	leader_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.0;
    leader_cell.phenotype.death.rates[apoptosis_index] = 0.0;
	leader_cell.phenotype.death.rates[necrosis_index] = 0.0;
	
    
	// Temperarily eliminating leader/follower signal	
	
	// Obviously missing - add later

	// turn on motility 
	leader_cell.phenotype.motility.is_motile = parameters.bools("leader_motility_mode"); 
	
    leader_cell.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("leader_adhesion");
	// std::cout<<leader_cell.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
    
	leader_cell.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("leader_repulsion");
    // std::cout<<leader_cell.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;

	// Temperarily eliminating leader/follower signal	

	//    leader_cell.phenotype.secretion.secretion_rates[1] = 50; // leader signal

    // modify ECM
    
    leader_cell.functions.custom_cell_rule = ecm_update_from_cell; // Only leaders can modify ECM (phenotype -> ECM)

	// set functions
	
	leader_cell.functions.update_migration_bias = chemotaxis_oxygen;//rightward_deterministic_cell_march; Use rightward deterministic march for march test. Set leader fraction to 1.0.
	
    leader_cell.functions.update_phenotype = leader_cell_phenotype_model; //NULL; // leader_cell_phenotype_model;

	if( parameters.ints("unit_test_setup") == 1)
	{
		leader_cell.phenotype.motility.persistence_time = 10.0;
		leader_cell.phenotype.motility.migration_speed = 1.0;
		leader_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
		leader_cell.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
		leader_cell.custom_data.add_variable( "max speed", "micron/min" , 1.0 ); // Maximum migration speed
		leader_cell.custom_data.add_variable( "min ECM motility density", "dimensionless", 0.0 );  // Minimum ECM density required for cell motility
		leader_cell.custom_data.add_variable( "max ECM motility density", "dimensionless", 1.0 );  // Maximum ECM density allowing cell motility
		leader_cell.custom_data.add_variable( "ideal ECM motility density", "dimensionless", 0.5 );  // Ideal ECM density cell motility
		leader_cell.custom_data.add_variable( "target ECM density", "dimensionless", 0.5 ); 
		if (parameters.ints("march_unit_test_setup") == 1){
		leader_cell.functions.update_migration_bias = rightward_deterministic_cell_march;
		}

		// std::cout<< "running unit test setup leader"<<std::endl;

	}
	
	// follower cells

	follower_cell = cell_defaults;
	follower_cell.name = "follower cell"; 
	follower_cell.type = 2;
    
    follower_cell.functions.update_phenotype = follower_cell_phenotype_model;//NULL;// follower_cell_phenotype_model;

	follower_cell.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("follower_adhesion");
	std::cout<<follower_cell.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
	follower_cell.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("follower_repulsion");
   	std::cout<<follower_cell.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
	follower_cell.phenotype.motility.is_motile = parameters.bools("follower_motility_mode");
	follower_cell.functions.update_migration_bias = ECM_informed_motility_update; // NEW!

	if( parameters.ints("unit_test_setup") == 1)
	{
		follower_cell.phenotype.motility.persistence_time = 10.0;
		follower_cell.phenotype.motility.migration_speed = 1.0;
		follower_cell.phenotype.mechanics.cell_cell_adhesion_strength = 0.0;
		follower_cell.phenotype.mechanics.cell_cell_repulsion_strength = 0.0;
		follower_cell.custom_data.add_variable( "max speed", "micron/min" , 1.0 ); // Maximum migration speed
		follower_cell.custom_data.add_variable( "min ECM motility density", "dimensionless", 0.0 );  // Minimum ECM density required for cell motility
		follower_cell.custom_data.add_variable( "max ECM motility density", "dimensionless", 1.0 );  // Maximum ECM density allowing cell motility
		follower_cell.custom_data.add_variable( "ideal ECM motility density", "dimensionless", 0.5 );  // Ideal ECM density cell motility
		follower_cell.custom_data.add_variable( "target ECM density", "dimensionless", 0.5 ); 

		std::cout<< "running unit test setup follower"<<std::endl;

	}
    
	std::cout<<"Follower cell migration speed "<<follower_cell.phenotype.motility.migration_speed <<std::endl;

	// Temperarily eliminating leader/follower signal

   	// follower_cell.phenotype.secretion.secretion_rates[2] = 50; // follower signal
    
	// Temperarily eliminating leader/follower signal
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters - see XML	
	
	// Add ECM structures
	microenvironment.add_density( "ECM", "dimensionless", 0.0 , 0.0 ); 
	microenvironment.add_density( "ECM anisotropy" , "dimensionless" , 0.0 , 0.0 ); 

	// Turn on gradients - oxygen chemotaxis
	
	default_microenvironment_options.calculate_gradients = false; 

	// Temperarily eliminating leader/follower signal (except here)
    
	// 50 micron length scale 
    microenvironment.add_density( "leader signal", "dimensionless", 1e5 , 1 );
    microenvironment.add_density( "follower signal", "dimensionless", 1e5 , 1 );

	// Temperarily eliminating leader/follower signal	
    
	// let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = true; 
	
	// set Dirichlet conditions 
	
	default_microenvironment_options.outer_Dirichlet_conditions = true;

	std::vector<double> bc_vector;

	if(parameters.ints("unit_test_setup")==1 && parameters.ints("march_unit_test_setup") == 0)
	{

		bc_vector = { 38.0 , 0.5, 1.0 , 0.0, 0.0 }; // 5% o2 , half max ECM , anisotropic, leader signal, follower signal
		default_microenvironment_options.X_range[0] = -500.0;
		default_microenvironment_options.X_range[1] = 500.0;
		default_microenvironment_options.Y_range[0] = -500.0;
		default_microenvironment_options.Y_range[1] = 500.0;

	}

	else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
	{
		bc_vector = { 38.0 , 0.5, parameters.doubles("initial_anisotropy"), 0.0, 0.0 }; // 5% o2 , half max ECM , anisotropic, leader signal, follower signal
		default_microenvironment_options.X_range[0] = -500.0;
		default_microenvironment_options.X_range[1] = 500.0;
		default_microenvironment_options.Y_range[0] = -500.0;
		default_microenvironment_options.Y_range[1] = 500.0;
	}

	else if(parameters.ints("unit_test_setup") == 0 && parameters.ints("march_unit_test_setup") == 0)
	{
		bc_vector = { 38.0 , parameters.doubles("initial_ECM_density") , parameters.doubles("initial_anisotropy") , 0.0, 0.0 };  // 5% o2 , half max ECM , isotropic, leader signal, follower signal
	}

	else
	{
		std::cout<<"ECM density and anisotropy not cell correctly!!!! WARNING!!!"<<std::endl;

		return;
	}
	
	
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
    
	// Temperarily eliminating leader/follower signal	
    
	// default_microenvironment_options.Dirichlet_condition_vector[1] = 0; // normoxic conditions
	// default_microenvironment_options.Dirichlet_condition_vector[2] = 0; // normoxic conditions
    
	initialize_microenvironment(); 

	// Trying to set the chemical gradient to be a starburst. Using the same code snippets as the ECM orientation. 

	if( parameters.ints("unit_test_setup") == 1)
	{

		for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
		{
			std::vector<double> position = microenvironment.mesh.voxels[i].center; 
			microenvironment.gradient_vector(i)[0] = { position[0],position[1],0}; 
			normalize(&microenvironment.gradient_vector(i)[0]);
			// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
		}
	}
	// run to get a decent starting conditoin (refers to code no longer present but will be added back for leader-follower signaling models)
	
	// now, let's set the leader signal to 1, so we don't hvae early swiching 
	/*
	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.density_vector(i)[1] = 1.0; 
	}
	*/

	// set up ECM density and anisotropy profile as needed
	int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 

	/*for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		std::vector<double> position = microenvironment.mesh.voxels[n].center; 
		if( fabs( position[0] ) > 200 || fabs( position[1] ) > 200 )
		{
			microenvironment(n)[ECM_density_index] = 0.0; 
			microenvironment(n)[ECM_anisotropy_index] = 1.0; 
		}
	}*/

	// set up ECM alignment 

	// <ECM_orientation_setup description="Specifies the initial ECM orientation: random, circular, starburt, oriented to the right, or oriented to the top" type="string" units="NA">circular</ECM_orientation_setup> parameters.string( "ECM_orientation_setup")
	std::vector<double> fiber_direction = { 1.0 , 0.0, 0.0 }; 
	ECM_fiber_alignment.resize( microenvironment.mesh.voxels.size() , fiber_direction );  

	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		// For random 2-D initalization 
		if(parameters.strings( "ECM_orientation_setup") == "random")
		{
			double theta = 6.2831853071795864769252867665590 * uniform_random(); 
			// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
			// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
			// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
			ECM_fiber_alignment[n] = {cos(theta), sin(theta), 0.0};
		}

		// for starburst initialization 
		else if(parameters.strings( "ECM_orientation_setup") == "starburt")
		{
			std::vector<double> position = microenvironment.mesh.voxels[n].center; 
			normalize( &position ); 
			ECM_fiber_alignment[n] =  { position[0],position[1],0}; // oriented out (perpindeicular to concentric circles)
			normalize(&ECM_fiber_alignment[n]);
		}

		// for circular initialization 
		else if(parameters.strings( "ECM_orientation_setup") == "circular")
		{
			std::vector<double> position = microenvironment.mesh.voxels[n].center; 
			normalize( &position );
			ECM_fiber_alignment[n] =  { position[1],-position[0],0}; // oriented in cirlce
			normalize(&ECM_fiber_alignment[n]);
		}

		else if(parameters.strings( "ECM_orientation_setup") == "oriented to the right")
		{
			ECM_fiber_alignment[n] =  { 1.0, 0.0, 0.0}; 
			normalize(&ECM_fiber_alignment[n]);
		}

		else if(parameters.strings( "ECM_orientation_setup") == "oriented to the top")
		{
			ECM_fiber_alignment[n] =  { 0.0, 1.0, 0.0}; 
			normalize(&ECM_fiber_alignment[n]);
		}

		else
		{
			std::cout<<"WARNING: NO ECM ORIENTATION SPECIFIED. FIX THIS!!!"<<std::endl;
			return;
		}
		
	}
	
	return; 
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

void setup_tissue( void )
{	

	if (parameters.ints("march_unit_test_setup") == 0)
	{
		if(parameters.strings("cell_setup") == "random")
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
			double sqrt_adhesion_to_repulsion_ratio = sqrt(parameters.doubles("follower_adhesion")/parameters.doubles("follower_repulsion"));
			
			double cell_spacing = (1 - sqrt_adhesion_to_repulsion_ratio);
			cell_spacing /= (0.5 * 1/cell_radius - 0.5 * sqrt_adhesion_to_repulsion_ratio/(relative_maximum_adhesion_distance * cell_radius));
			
			Cell* pCell = NULL; 
			
			double x = 0.0;
			double x_outer = tumor_radius; 
			double y = 0.0;

			double leader_cell_fraction;
			
			if( parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 0)
			{
				leader_cell_fraction = 0.0;
			}

			else
			{
				leader_cell_fraction = parameters.doubles("initial_leader_cell_fraction"); // 0.2;
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
					if( UniformRandom() < leader_cell_fraction )
					{ pCell = create_cell(leader_cell); }
					else
					{ pCell = create_cell(follower_cell);}
						
					pCell->assign_position( x , y , 0.0 );
					
					if( fabs( y ) > 0.01 )
					{
						if( UniformRandom() < leader_cell_fraction )
						{ pCell = create_cell(leader_cell); }
						else
						{ pCell = create_cell(follower_cell); }
						pCell->assign_position( x , -y , 0.0 );
					}
					
					if( fabs( x ) > 0.01 )
					{ 
						if( UniformRandom() < leader_cell_fraction )
						{ pCell = create_cell(leader_cell); }
						else
						{ pCell = create_cell(follower_cell); }
						pCell->assign_position( -x , y , 0.0 );
						
						if( fabs( y ) > 0.01 )
						{
							if( UniformRandom() < leader_cell_fraction )
							{ pCell = create_cell(leader_cell); }
							else
							{ pCell = create_cell(follower_cell); }
							
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

		

		/******************************************3D Spheroid initialization***************************************/

		/*To come later*/

		/************************************Line of cells at y = 0, x > 0 initialization***************************************/

		else if(parameters.strings("cell_setup") == "cells at y = 0")
		{

			Cell* pCell = NULL; 
			int n = default_microenvironment_options.X_range[0] + 10.0; 
			while( n <= default_microenvironment_options.X_range[1] )
			{
				pCell = create_cell(follower_cell); 
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
				pCell = create_cell(leader_cell); 
				pCell->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
				n = n + 10.0;
			}
			std::cout<<"Cell's placed at left boundary for march test"<<std::endl;
		}

		else
		{
			std::cout<<"WARNING!!! NO CELL SETUP SPECIFIED. SEE DOCUMENTATION"<<std::endl;
		}
	}

	else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
	{
		Cell* pCell = NULL; 
		int n = default_microenvironment_options.X_range[0] + 10.0; 
		while( n <= default_microenvironment_options.X_range[1] - 10.0 )
		{
			pCell = create_cell(leader_cell); 
			pCell->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
			n = n + 10.0;
		}
		std::cout<<"Cell's placed at left boundary for march test"<<std::endl;

	}

	else
	{
		std::cout<<"RUN MODE (TESTING OR NOT TESTING) NOT SPECIFIED!!!!! WARNING!!!!"<<std::endl;
	}
	

	return; 
}

double dot_product( const std::vector<double>& v , const std::vector<double>& w )
{
 double out = 0.0; 
 for( unsigned int i=0 ; i < v.size() ; i++ )
 { out += ( v[i] * w[i] ); }
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

void ECM_informed_motility_update( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	
	if(phenotype.death.dead == true)
	{
		
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		std::cout<<2<<std::endl;
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
	
	// sample ECM 
	double ECM_density = pCell->nearest_density_vector()[ECM_density_index]; 
	double a = pCell->nearest_density_vector()[ECM_anisotropy_index]; 
	int n = pCell->get_current_voxel_index();
	std::vector<double> f = ECM_fiber_alignment[n];
	
	
	/****************************************Begin new migration direction update****************************************/

	// Select random direction (for random portion of motility vector) and begin building updated motility direction vector
	// (note - there is NO memeory of previous direction in this model - previous ECM-based motility used the current 
	// velocity vector to build off, not a random one - this could produce divergent behaviors between models)

	// See lab note book for more notes - MUST start with random vector. In the old method I defintiely used the previous motility vector in the method, but makes no sense here!

	// get random vector - cell's "intended" or chosen random direction
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_random = { cos(angle) , sin(angle) , 0.0 };

	// get vector for chemotaxis (sample uE)
	std::vector<double> chemotaxis_grad = pCell->nearest_gradient(o2_index);

	std::cout<<"D chemo"<<chemotaxis_grad<<std::endl;

	normalize( &chemotaxis_grad ); 

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)
	std::vector<double> d_motility = (1-pCell->custom_data[chemotaxis_bias_index])*d_random + pCell->custom_data[chemotaxis_bias_index]*chemotaxis_grad;
	normalize( &d_motility ); 


	std::cout<<"D motility "<<d_motility<<std::endl;

	// to determine direction along f, find part of d_choice that is perpendicular to f; 
	std::vector<double> d_perp = d_motility - dot_product(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	std::cout<<"D perp"<<d_perp<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product( d_motility , d_perp ); 
	double c_2 = dot_product( d_motility, f ); 

	std::cout<<"D_mot dot d_perp c_1 = "<<c_1<<std::endl;
	std::cout<<"D_mot dot f c_2 = "<<c_2<<std::endl;

	// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy

	double gamma = pCell->custom_data[ECM_sensitivity_index] * a; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.

	phenotype.motility.migration_bias_direction = (1.0-gamma)*c_1*d_perp + c_2*f;
	std::cout<<"migration_bias_direction before normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
	if(parameters.bools("normalize_ECM_influenced_motility_vector") == true)
	{
		normalize( &phenotype.motility.migration_bias_direction ); 
	}
	std::cout<<"migration_bias_direction after normalization"<<phenotype.motility.migration_bias_direction<<std::endl;
	phenotype.motility.migration_bias = 1.0; // MUST be set at 1.0 so that standard update_motility function doesn't add random motion. 

	double magnitude = norm( phenotype.motility.motility_vector);

	std::cout<<"Cell speed is "<< magnitude<<std::endl;

	if(magnitude > 0.00000001)
	{
		std::cout<<"Cell is moving!!!!"<<std::endl;
	}

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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low));
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high));
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

    // leader_cell.phenotype.death.rates[apoptosis_index] = 0.0;
	// leader_cell.phenotype.death.rates[necrosis_index] = 0.0;
	
	// END New speed update 

	// Old parabolic update  - just the one line!

	// pCell->phenotype.motility.migration_speed = (-(4.0)*pow((ECM_density-0.5),2.0) + 1.0) * pCell->custom_data[max_cell_speed_index];
	//   std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;


	/*********************************************END speed update***************************************************/


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
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	
	phenotype.motility.is_motile = true; 
	phenotype.motility.migration_bias = parameters.doubles("oxygen_migration_bias_for_leaders"); //0.95;
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(o2_index);

	normalize( &( phenotype.motility.migration_bias_direction ) );

	if(phenotype.death.dead == true)
	{
		phenotype.motility.is_motile=false;
		pCell->functions.update_phenotype = NULL;
		pCell->functions.update_migration_bias = NULL;


	}	

   	// std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;
	
	return; 
}

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell , Phenotype& phenotype , double dt ) 
{	
	// o2-based birth and death, nothing more 
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	return; 
}


// void follower_cell_phenotype_model0( Cell* pCell , Phenotype& phenotype , double dt )
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
	// leaders are blue 
	std::vector< std::string > output( 4, "black" ); 
	
	if( pCell->type == 1 )
	{ 
		output[0] = "blue"; 
		output[2] = "blue"; 
		return output; 
	} 


	// followers are yellow
    if( pCell->type == 2 )
    {
        output[0] = "yellow";
        output[2] = "yellow";
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

void leader_cell_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt )
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


void follower_cell_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt )
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
		// std::cout<<"Follower is dead"<<std::endl;
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

// commented out by JPM 11.04.18 to eliminate chance of seg fault from currnetly non-exsistent fields.

/* void switching_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt )
{
   static int hypoxic_i = pCell->custom_data.find_variable_index( "hypoxic switch value" );
   static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" );


   int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
   int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );

   double pO2 = (pCell->nearest_density_vector())[oxygen_i]; // PhysiCell_constants::oxygen_index];

   // set death and birth
   update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);

   // if a leader and now happy, switch to follower
   if( pO2 > pCell->custom_data[hypoxic_i] && pCell->type == 1 )
   {
       pCell->convert_to_cell_definition( follower_cell );

   }

   // if a folower and now unhappy, switch to leader, if not supressed

   if( pO2 < pCell->custom_data[hypoxic_i] && pCell->type == 2 )
   {
       if( (pCell->nearest_density_vector())[1] < 0.1 )
       {
        pCell->convert_to_cell_definition( leader_cell );
       }
   }

   return;
} */

long fibonacci(unsigned n) // just being used for timing. 
{
    if (n < 2) return n;
    return fibonacci(n-1) + fibonacci(n-2);
}

void ecm_update_from_cell(Cell* pCell , Phenotype& phenotype , double dt)
{

	// Find correct fields
	static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target ECM density");
    
    // Cell-ECM density interaction
    double ECM_density = pCell->nearest_density_vector()[ECM_density_index]; 
    double r = 1.0;
    
    pCell->nearest_density_vector()[ECM_density_index] = ECM_density + r * dt  * (pCell->custom_data[Cell_ECM_target_density_index] - ECM_density);
    
    // END Cell-ECM density interaction
    
    // Cell-ECM Fiber realingment

	// Get index for accessing the ECM_fiber_alignment data structure and then copy the correct value
	int n = pCell->get_current_voxel_index();
	std::vector<double> ECM_orientation = ECM_fiber_alignment[n];

	double anisotropy = pCell->nearest_density_vector()[ECM_anisotropy_index]; 
    double migration_speed = pCell->phenotype.motility.migration_speed;

    double r_0 = 1/1.0*migration_speed; // 1/10.0 // min-1 // NOTE!!! on 08.06.18 run - this wasn't multiplied by migration_speed!!! should be the same but worth noting!!!!

    double r_realignment = r_0 * (1-anisotropy);

	double ddotf;
	std::vector<double> norm_cell_motility;
	norm_cell_motility.resize(3,0.0);
	norm_cell_motility = phenotype.motility.motility_vector;
	normalize(&norm_cell_motility);

	ddotf = dot_product(ECM_orientation, norm_cell_motility);
	
	ECM_orientation = sign_function(ddotf) * ECM_orientation; // flips the orientation vector so that it is aligned correctly with the moving cell for proper reoirentation later.
	
	std::vector<double> f_minus_d;
	f_minus_d.resize(3,0.0);

	// f_minus_d = ECM_orientation - norm_cell_motility; // Fix this later

	for(int i = 0; i < 3; i++)
	{
		f_minus_d[i] = ECM_orientation[i] - norm_cell_motility[i]; // 06.05.19 - fixed 
		ECM_fiber_alignment[n][i] -= dt * r_realignment * f_minus_d[i]; 
	}
	
	
    normalize(&(ECM_fiber_alignment[n])); // why by reference??


    // End Cell-ECM Fiber realingment
    
    // Cell-ECM Anisotrophy Modification
    
    double r_a0 = 1.0/1000.0; // min-1
    double r_anisotropy = r_a0 * migration_speed;
   	
    pCell->nearest_density_vector()[ECM_anisotropy_index] = anisotropy + r_anisotropy * dt  * (1- anisotropy);
    
    // END Cell-ECM Anisotropy Modification
    
    return;

}
    
void write_ECM_Data_matlab( std::string filename )

{

    int number_of_data_entries = microenvironment.number_of_voxels();

    int size_of_each_datum = 11;

	
	static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 

    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "ECM_Data" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )

    {

	    fwrite( (char*) &( microenvironment.mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp ); // 1

        fwrite( (char*) &( microenvironment.mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp ); // 2

        fwrite( (char*) &( microenvironment.mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp ); //3
		
		fwrite( (char*) &( microenvironment.density_vector(i)[ECM_anisotropy_index]), sizeof(double) , 1 , fp ); // 4
	
        fwrite( (char*) &( microenvironment.density_vector(i)[ECM_density_index]), sizeof(double) , 1 , fp ); // 5

        fwrite( (char*) &( ECM_fiber_alignment[i][0]), sizeof(double) , 1 , fp ); // 6

        fwrite( (char*) &( ECM_fiber_alignment[i][1]), sizeof(double) , 1 , fp ); // 7

        fwrite( (char*) &( ECM_fiber_alignment[i][2]), sizeof(double) , 1 , fp ); // 8

		fwrite( (char*) &( microenvironment.gradient_vector(i)[0][0]), sizeof(double) , 1 , fp ); // 9

		fwrite( (char*) &( microenvironment.gradient_vector(i)[0][1]), sizeof(double) , 1 , fp ); // 10

		fwrite( (char*) &( microenvironment.gradient_vector(i)[0][2]), sizeof(double) , 1 , fp ); // 11

    }



    fclose( fp );



    return;

}
