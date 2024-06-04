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

#include "./AMIGOS-invasion_uncoupled.h"
#include "./extracellular_matrix.h"
#include "./cell_ECM_interactions.h"
#include <chrono>  // for high_resolution_clock - https://www.pluralsight.com/blog/software-development/how-to-measure-execution-time-intervals-in-c--
// std::vector< std::vector<double> > ECM_fiber_alignment; 
unsigned long long int counter=0; // counter for calculating average for the ad hoc timing I am doing ... 
int time_total = 0;

// Cell_Definition leader_cell; 
// Cell_Definition follower_cell; 

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
	
	create_default_ECM_compatible_agent(); // in cell_ECM_interactions.cpp. Sets custom velocity function (cell-ECM motility interaction) and custom cell rule (ECM remodeling).
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
       Cell rule definitions 
	*/

	setup_cell_rules(); 

	/* 
	   Put any modifications to individual cell definitions here. 
	   
	   This is a good place to set custom functions. 
	*/ 
	
    Cell_Definition* leader_cell = find_cell_definition("leader cell");	
	Cell_Definition* follower_cell = find_cell_definition("follower cell");	
	

	leader_cell->functions.update_velocity = custom_update_cell_velocity;

	if ( parameters.strings("ecm_update_model") == "ecm_update_from_cell_motility_vector")
    {
    	leader_cell->functions.custom_cell_rule = ECM_remodeling_function;  // combined_ECM_remodeling_and_speed_update; // In cell_ECM_interactions.cpp
    }
	else if( parameters.strings("ecm_update_model") == "ecm_update_from_cell_velocity_vector")
	{
		leader_cell->functions.custom_cell_rule = ecm_update_from_cell_velocity_vector; // Only leaders can modify ECM (phenotype -> ECM)
	}
	else
	{
		std::cout<<"no reorientation model specified!!@!!!!!! Halting!!!!!!"<<std::endl;
		abort();
		return;
	}
	
	// leader_cell->functions.update_migration_bias = ECM_and_chemotaxis_based_cell_migration_update; //in cell_ECM_interactions.cpp
	
    leader_cell->functions.update_phenotype = NULL; // leader_cell_phenotype_model;

	// if( parameters.ints("unit_test_setup") == 1)
	// {
        if (parameters.ints("march_unit_test_setup") == 1)
        {
            leader_cell->functions.update_migration_bias = rightward_deterministic_cell_march;
		}
    // }

    //--------- now follower:

    follower_cell->functions.update_phenotype = NULL;// follower_cell_phenotype_model;
	// follower_cell->functions.custom_cell_rule = ECM_based_speed_update; // includes both speed and ECM remodeling

	follower_cell->functions.update_velocity = custom_update_cell_velocity;

// <cell_motility_ECM_interaction_model_selector type="string" units="" description="follower chemotaxis/no follower hysteresis, follower hysteresis/no follower chemotaxis">follower chemotaxis/no follower hysteresis<

    // rwh: doing this one:
    if ( parameters.strings("cell_motility_ECM_interaction_model_selector") == "follower chemotaxis/no follower hysteresis" || parameters.ints("unit_test_setup") == 1)
	{
		// follower_cell->functions.update_migration_bias = ECM_and_chemotaxis_based_cell_migration_update; // In cell_ECM_interactions.cpp
		std::cout<<"Selection not currently supported" << std::endl;   // <------ rwh
		std::cout<<"Using default chemotaxis" << std::endl;   // <------ rwh
		// exit(-1);
	}
	else if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "follower hysteresis/no follower chemotaxis")
	{
		follower_cell->functions.update_migration_bias = ECM_informed_motility_update_model_w_memory;
		// follower_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
		// void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
		std::cout<<"I selected follower hysteresis" << std::endl;
				std::cout<<"Selection not currently supported" << std::endl;   // <------ rwh
		exit(-1);
	}
	else if( parameters.strings("cell_motility_ECM_interaction_model_selector") == "follower chemotaxis with variable follower speed")
	{
		follower_cell->functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
		// follower_cell.functions.update_migration_bias = ECM_informed_motility_update_w_chemotaxis_w_variable_speed;
		// void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
		std::cout<<"I selected follower chemotaxis with variable follower speed" << std::endl;
				std::cout<<"Selection not currently supported" << std::endl;   // <------ rwh
		exit(-1);
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

ECM ecm;

// END ECM STUFF

void setup_microenvironment( void )
{
	// May reintroduce code to set up the chemical field later. ery useful for doing unit testing of the chemotaxis code.
	// bool chemical_field_setup_specified_bool = (parameters.strings("chemical_field_setup")== "starburst" ||  parameters.strings("chemical_field_setup")== "vertical up" || parameters.strings("chemical_field_setup") == "horizontal right" || parameters.strings("chemical_field_setup") == "angle" || parameters.strings("chemical_field_setup") == "none");

	// if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "starburst")
	// {

	// 	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	// 	{
	// 		std::vector<double> position = microenvironment.mesh.voxels[i].center; 
	// 		microenvironment.gradient_vector(i)[0] = { position[0],position[1],0}; 
	// 		normalize(&microenvironment.gradient_vector(i)[0]);
	// 		// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
	// 	}
	// }

	// if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "vertical up")
	// {
	// 	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	// 	{
	// 		// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
	// 		microenvironment.gradient_vector(i)[0] = { 0,1,0}; 
	// 		normalize(&microenvironment.gradient_vector(i)[0]);
	// 		// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
	// 	}
	// }

	// if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "horizontal right")
	// {
	// 	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	// 	{
	// 		// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
	// 		microenvironment.gradient_vector(i)[0] = { 1,0,0}; 
	// 		normalize(&microenvironment.gradient_vector(i)[0]);
	// 		// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
	// 	}
	// }

	// if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "angle")
	// {
	// 	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	// 	{
	// 		// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
	// 		microenvironment.gradient_vector(i)[0] = { cos ( parameters.doubles("angle_of_chemical_field_gradient") * PhysiCell_constants::pi/180) , sin ( parameters.doubles("angle_of_chemical_field_gradient") * PhysiCell_constants::pi/180),0}; 
	// 		normalize(&microenvironment.gradient_vector(i)[0]);
	// 		// std::cout<<microenvironment.gradient_vector(i)[0][2]<<std::endl;
	// 	}
	// }

	// if( parameters.ints("unit_test_setup") == 1 && parameters.strings("chemical_field_setup") == "none")
	// {
	// 	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	// 	{
	// 		// std::vector<double> position = microenvironment.mesh.voxels[i].center; 
	// 		microenvironment.gradient_vector(i)[0] = { 0,0,0}; 
	// 	}
	// }
	// else if( parameters.ints("unit_test_setup") == 1 && chemical_field_setup_specified_bool == 0)
	// {
	// 	std::cout<<"WARNING: NO CHEMICAL FIELD ORIENTATION SPECIFIED for unit testing. FIX THIS!!!"<<std::endl;
	// 	std::cout<<"Halting program!!!"<<std::endl;
	// 	abort();
	// 	return;
	// }

	// to note it, previously lots of effort was put into keeping simulations easy to run by using parameters to specify domain size and many other aspects
	// of a simulation. With the migration to the full XML to specificy simualtions and Studio, this need has drastically reduced. One maybe still see some of this around though

	// set domain parameters 
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_extracellular_matrix( void )
{
	// DEPENDS ON MICROENVIRONMENT - CALL SETUP MICROENVIRONEMNT FIRST!!!!!

	ecm.ecm_mesh.resize(default_microenvironment_options.X_range[0], default_microenvironment_options.X_range[1] , 
	default_microenvironment_options.Y_range[0], default_microenvironment_options.Y_range[1],default_microenvironment_options.Z_range[0], default_microenvironment_options.Z_range[1], \
	parameters.doubles("ECM_dx"), parameters.doubles("ECM_dy"),parameters.doubles("ECM_dz"));
	ecm.resize_ecm_units_from_ecm_mesh();

	ecm.ecm_mesh.display_information(std::cout );

	for( int n = 0; n < ecm.ecm_mesh.voxels.size() ; n++ )
	{
		// For random 2-D initalization 
		if(parameters.strings( "ECM_orientation_setup") == "random")
		{
			double theta = 6.2831853071795864769252867665590 * UniformRandom(); 
			ecm.ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0};
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

		else if(parameters.ints("unit_test_setup") == 0 && parameters.ints("march_unit_test_setup") == 1)
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
	}
}

void run_biotransport( double t_max ) // used to set up initial chemical conditions to prevent unwanted phenotype switching due to starting the simulation
{
	std::cout << "working on initial conditions .. " << std::endl; 
	double t = 0.0;
	
	for( int i=0 ; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i];
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
	
	std::cout << "Preconditioned chemical microenvironemnt and did agent secretion for " << t_max << " minutes." << std::endl; 
	return; 
}

void set_cell_motility_vectors( void )
{
	for( int i=0 ; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i];
		pCell->update_motility_vector( 100000000 );
		// pCell->custom_update_motility_vector( 100000000 );
		// std::cout<<"Initial Motility vector "<<pCell->phenotype.motility.motility_vector<<std::endl;
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
	// SeedRandom(0);
    static Cell_Definition* leader_cell = find_cell_definition("leader cell");	
	static Cell_Definition* follower_cell = find_cell_definition("follower cell");	

    Cell* pC;

	if (parameters.ints("march_unit_test_setup") == 0)
	{
		if(parameters.strings("cell_setup") == "single")
		{
			// pC = create_cell( *follower_cell );
			pC = create_cell( *leader_cell );
			pC->assign_position(0.0, 0.0, 0.0);
		}
		else if(parameters.strings("cell_setup") == "random")
		{
			std::cout<<"string worked"<<std::endl;
			/*******************************************Random initialization****************************************/
			for( int n = 0 ; n < 200 ; n++ )
			{
				pC = create_cell();   //rwh: what cell type does this create? Default? What's desired??
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
			if(parameters.doubles("follower_repulsion") == 0)
			{
				sqrt_adhesion_to_repulsion_ratio = 0.632455; // value for adhesion = 10 and repulsion = 25.0 - "parameter set 21"
			}
			else
			{
				sqrt_adhesion_to_repulsion_ratio = sqrt(parameters.doubles("follower_adhesion")/parameters.doubles("follower_repulsion"));
			} 

			double cell_spacing = (1 - sqrt_adhesion_to_repulsion_ratio);
			cell_spacing /= (0.5 * 1/cell_radius - 0.5 * sqrt_adhesion_to_repulsion_ratio/(relative_maximum_adhesion_distance * cell_radius));
			
			Cell* pC = NULL; 
			
			double x = 0.0;
			double x_outer = tumor_radius; 
			double y = 0.0;

			double leader_cell_fraction;
			if( parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 0)
			{
				leader_cell_fraction = 0.0;
				cell_spacing = 1.90 * cell_radius;
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
					{ pC = create_cell( *leader_cell ); }
					else
					{ pC = create_cell( *follower_cell );}
						
					pC->assign_position( x , y , 0.0 );
					
					if( fabs( y ) > 0.01 )
					{
						if( UniformRandom() < leader_cell_fraction )
						{ pC = create_cell( *leader_cell ); }
						else
						{ pC = create_cell( *follower_cell ); }
						pC->assign_position( x , -y , 0.0 );
					}
					
					if( fabs( x ) > 0.01 )
					{ 
						if( UniformRandom() < leader_cell_fraction )
						{ pC = create_cell( *leader_cell ); }
						else
						{ pC = create_cell( *follower_cell ); }
						pC->assign_position( -x , y , 0.0 );
						
						if( fabs( y ) > 0.01 )
						{
							if( UniformRandom() < leader_cell_fraction )
							{ pC = create_cell( *leader_cell ); }
							else
							{ pC = create_cell( *follower_cell ); }
							
							pC->assign_position( -x , -y , 0.0 );
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

		/************************************Circle of cells at R = 300 initialization***************************************/
		else if(parameters.strings("cell_setup") == "circle of cells")
		{
			double theta2 = 0.0;
			for (int a = 0; a<42; a++)
			{
				pC = create_cell( *follower_cell ); 
				pC->assign_position( 300 * cos(theta2) , 300 * sin(theta2) , 0.0 );
				theta2 += 0.14959952;
			}
		}

		/************************************Line of cells at y = 0, x > 0 initialization***************************************/
		else if(parameters.strings("cell_setup") == "cells at y = 0")
		{
			int n = default_microenvironment_options.X_range[0] + 10.0; 
			while( n <= default_microenvironment_options.X_range[1] )
			{
				pC = create_cell( *follower_cell );

				// To prevent droping cells in areas of high ECM curvature. 
				while(abs(n) < 70)
				{n = n + 30.0;}

				pC->assign_position( n , 0.0 , 0.0 );
				n = n + 30.0;
			}
			std::cout<<"Cell's placed at y = 0"<<std::endl;
		}

		/******************************************Line of cells at x = left boundary + 10 initialization***************************************/
		else if(parameters.strings("cell_setup") == "cells at left boundary/march")
		{
			
			int n = default_microenvironment_options.X_range[0] + 10.0; 
			while( n <= default_microenvironment_options.X_range[1] - 10.0 )
			{
				if (parameters.ints("march_unit_test_setup") == 1)
				{pC = create_cell( *leader_cell );}

				else 
				{pC = create_cell( *follower_cell );}
				pC->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
				n = n + 20.0;
			}
			std::cout<<"Cell's placed at left boundary for march test 000"<<std::endl;
		}

		else
		{
			std::cout<<"WARNING!!! NO CELL SETUP SPECIFIED. SEE DOCUMENTATION and FIX"<<std::endl;
			std::cout<<"Halting program!!!"<<std::endl;
			abort();
			return;
		}
	}

	// else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
	else if (parameters.ints("unit_test_setup") == 0 && parameters.ints("march_unit_test_setup") == 1)
	{
		int n = default_microenvironment_options.X_range[0] + 5.0; 
		while( n <= default_microenvironment_options.X_range[1] - 10.0 )
		{
			pC = create_cell( *leader_cell ); 
			pC->assign_position( default_microenvironment_options.X_range[0] + 10.0 , n , 0.0 );
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


//rwh: this is being called
void ECM_informed_motility_update_w_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	// std::cout<<"ECM_informed_motility_update_w_chemotaxis(): cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;  //rwh
	
	if(phenotype.death.dead == true)
	{
		phenotype.motility.is_motile = false;
		pCell->functions.update_migration_bias = NULL;
		pCell->functions.update_phenotype = NULL;
		// std::cout<<2<<std::endl;
		std::cout<<"         Cell is dead"<<std::endl;
	}
	// Updates cell bias vector and cell speed based on the ECM density, anisotropy, and fiber direction
	
	// find location of variables and base parameter values
	// static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	// static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	

            // <custom_data>
            //     <min_ECM_motility_density conserved="false" units="dimensionless" description="">0.0</min_ECM_motility_density>
            //     <max_ECM_motility_density conserved="false" units="dimensionless" description="">1.0</max_ECM_motility_density>
            //     <ideal_ECM_motility_density conserved="false" units="dimensionless" description="">0.5</ideal_ECM_motility_density>
            //     <max_speed conserved="false" units="micron/min" description="">0.5</max_speed>
            //     <chemotaxis_bias conserved="false" units="dimensionless" description="">0.05</chemotaxis_bias>
            //     <ECM_sensitivity conserved="false" units="dimensionless" description="">1.0</ECM_sensitivity>
            //     <hypoxic_switch_value conserved="false" units="mmHg" description="">38</hypoxic_switch_value>
            //     <target_ECM_density conserved="false" units="dimensionless" description="">0.5</target_ECM_density>
            //     <Base_hysteresis_bias conserved="false" units="dimensionless" description="">1.0</Base_hysteresis_bias>
            //     <previous_anisotropy conserved="false" units="dimensionless" description="">0</previous_anisotropy>
            //     <Anisotropy_increase_rate conserved="false" units="1/min" description="">0.004</Anisotropy_increase_rate>
            //     <Fiber_realignment_rate conserved="false" units="1/min" description="">4.0</Fiber_realignment_rate>
            // </custom_data>

    // std::cout << "------------  ECM_informed_motility_update_w_chemotaxis():\n";  //rwh

	static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max_speed" ); 
    // std::cout << "        static int max_cell_speed_index = " <<max_cell_speed_index << std::endl;
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis_bias");
	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
    // std::cout << "        static int ECM_sensitivity_index = " <<ECM_sensitivity_index << std::endl;
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

	std::cout<<"D random "<<d_random<<std::endl;

	// get vector for chemotaxis (sample uE)
	std::vector<double> chemotaxis_grad = pCell->nearest_gradient(o2_index);

	std::cout<<"D chemo"<<chemotaxis_grad<<std::endl;

	normalize( &chemotaxis_grad ); 

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)

	// New bias - bias such that the agents can more closely follow the gradient IF the written signals are stronger. 

	std::vector<double> d_motility;

	// d_motility = {0,0,0};
	
	if (parameters.ints("link_anisotropy_and_bias") == 0)   // rwh: this
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
	else  // rwh: this (bool is false); speed is the length of the bias direction
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
        std::cout << "rwh #1: setting speed = 0 due to ECM_density (" << ECM_density << ") <=  rho_low (" << rho_low << ")\n";
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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
        // std::cout << "rwh #2: setting speed = " << pCell->phenotype.motility.migration_speed << std::endl;
	}
	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{
		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
        std::cout << "rwh #3: setting speed = " << pCell->phenotype.motility.migration_speed << std::endl;
	}
	else //if (ECM_density >= rho_high)
	{
        std::cout << "rwh #4: setting speed = 0\n";
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

	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
	return; 
}

/* To eliminate memory, set b_h to zero. To eliminate ECM influence, set a to 0 (permanently) or ECM senstiivity to zero,
   This model uses memory from previous motility vector AND previous anisotropy location. */

//rwh: commenting out for now; one would need to replace spaces with underscores on custom data vars!

void ECM_informed_motility_update_model_w_memory ( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout<<"ECM_informed_motility_update_model_w_memory() called. Exiting until it is fixed!\n";  //rwh
    std::exit(-1);  //rwh

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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
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

	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
	return; 
}

// Sets speed based on voxel anisotropy (instead of magnitude of motility vector). With chemotaxis senstivity set to 1, cells will chemotax, but have speed vary with only anisotropy. Direction WILL NOT be effeced by ECM, only speed. 
// Only fiber alignemnt (self alignemnt) impacts the motility. 

void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout<<"ECM_informed_motility_update_w_chemotaxis_w_variable_speed() called. Exiting until it is fixed!\n";  //rwh
    std::exit(-1);  //rwh

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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
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
	int n = default_microenvironment_options.X_range[0] + 5.0;

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
		std::cout<<"Cell is dead"<<std::endl;


	}	
	// std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
   	// std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;
	
	return; 
}

std::vector<std::string> AMIGOS_invasion_coloring_function( Cell* pCell )
{
	// leaders are blue 
	std::vector< std::string > output( 4, "black" ); 
	

	if( pCell->type == 0 )
	{ 
		output[0] = "blue"; 
		output[2] = "blue"; 
		if(parameters.ints("unit_test_setup")==0)
		{return output;}
		
		// Return red for 20% of followers if unit test is called for	

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

	if( pCell->type == 1 )
    {
        output[2] = "yellow";	
		output[0] = "yellow";
		
		// Return yellow for followers and exit statement

		if(parameters.ints("unit_test_setup")==0)
		{return output;}

		// Return red for 20% of followers if unit test is called for	

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

long fibonacci(unsigned n) // just being used for timing. 
{
    if (n < 2) return n;
    return fibonacci(n-1) + fibonacci(n-2);
}

// uses cell motility vector for fiber reorientation
void ecm_update_from_cell_motility_vector(Cell* pCell , Phenotype& phenotype , double dt) 
{
	// Find correct items
	// static int ECM_density_index = microenvironment.find_density_index( "ECM" ); 
	// static int ECM_anisotropy_index = microenvironment.find_density_index( "ECM anisotropy" ); 
	std::vector<double> cell_position = pCell->position;
	// std::cout<<cell_position<<std::endl;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	// std::cout<<nearest_ecm_voxel_index<<std::endl;
	// std::cin.get();

    //rwh: replace spaces with underscores
	// static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target ECM density");
	// static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy increase rate");
	// static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber realignment rate");

	static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target_ECM_density");
	static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy_increase_rate");
	static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber_realignment_rate");
    
    // Cell-ECM density interaction

    double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;
    double r = 1.0;
    
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

    //rwh: replace spaces with underscores
	// static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target ECM density");
	// static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy increase rate");
	// static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber realignment rate");

	static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target_ECM_density");
	static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy_increase_rate");
	static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber_realignment_rate");
    
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
