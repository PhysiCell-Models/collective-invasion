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
#include "./cell_ECM_interactions.h"
#include <chrono>  // for high_resolution_clock - https://www.pluralsight.com/blog/software-development/how-to-measure-execution-time-intervals-in-c--
// Cell_Definition fibroblast; 
// Cell_Definition cancer_cell; 
ECM ecm;
unsigned long long int counter=0; // counter for calculating average for the ad hoc timing I am doing ... 
int time_total = 0;

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

	fibroblast->functions.custom_cell_rule = ECM_remodeling_function; // Only leaders can modify ECM (phenotype -> ECM)
	
	fibroblast->functions.update_migration_bias = fibroblast_ECM_informed_motility_update_w_chemotaxis;//rightward_deterministic_cell_march; Use rightward deterministic march for march test. Set leader fraction to 1.0.
	
    fibroblast->functions.update_phenotype = NULL; // leader_cell_phenotype_model;

	cancer_cell->functions.custom_cell_rule = custom_cancer_cell_ECM_remodeling_and_adhesion_function;

	cancer_cell->functions.update_migration_bias = cancer_cell_ECM_informed_motility_update_w_chemotaxis;

	cancer_cell->functions.update_phenotype = NULL;// follower_cell_phenotype_model;

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

			if(parameters.bools( "heterogeneous_invasive_spheroid") == true)
				{if(ecm.ecm_mesh.voxels[n].center[1] > -parameters.doubles("tumor_radius"))
					{
						double theta = 6.2831853071795864769252867665590 * uniform_random(); 
						// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
						// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
						// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
						ecm.ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0};
						ecm.ecm_voxels[n].density = 0.25;
					}
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

double dot_product_ext_old( const std::vector<double>& v , const std::vector<double>& w )
{
	double out = 0.0; 
	for( unsigned int i=0 ; i < v.size() ; i++ )
	{ out += ( v[i] * w[i] ); }

	if( fabs(out) < 1e-10)
	{out = 0.0;}

	return out; 
}

double sign_function_old (double number)
{
	// double sign = 0.0
	if (number<0)
	{ return -1.0;}

	else
	{ return 1.0;}

}

/* To eliminate chemotaxis, set chemotaxis bias to zero. To eliminate ECM influence, set a to 0 (permanently) or ECM senstiivity to zero */

void cancer_cell_ECM_informed_motility_update_w_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
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
	
	static int link_anisotropy_and_bias_index = pCell->custom_data.find_variable_index( "link_anisotropy_and_bias" );
	    if (link_anisotropy_and_bias_index < 0) 
    {
        std::cout << "        static int link_anisotropy_and_bias_index = " <<link_anisotropy_and_bias_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
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
	
	if (pCell->custom_data[link_anisotropy_and_bias_index] < 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 0
	{
		d_motility = (1-a) * d_random + a * chemotaxis_grad;
		// std::cout<<" I am coupled"<<std::endl;
	}

	else if (pCell->custom_data[link_anisotropy_and_bias_index] > 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 1
	{
		// NON-ECM linked way to signal 
		d_motility = (1-pCell->custom_data[chemotaxis_bias_index])*d_random + pCell->custom_data[chemotaxis_bias_index]*chemotaxis_grad;
		// std::cout<<" I am UNcoupled"<<std::endl;
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
	std::vector<double> d_perp = d_motility - dot_product_ext(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product_ext( d_motility , d_perp ); 
	double c_2 = dot_product_ext( d_motility, f ); 

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

	// std::cout<<"migration speed(l1236) = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	// double magnitude = norm( phenotype.motility.motility_vector);	

	// std::cout<<"Magnitutude of motility vector is "<< magnitude<<std::endl;

	// if(magnitude > 0.00000001)
	// {
	// 	std::cout<<"Cell is moving!!!!"<<std::endl;
	// }

	/****************************************END new migration direction update****************************************/


	// /*********************************************Begin speed update***************************************************/
	
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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	}

	// else //if (ECM_density >= rho_high)
	// {
	// 	pCell->phenotype.motility.migration_speed = 0.0;
	// }

	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// std::cout<<"cell adhesion = "<<pCell->phenotype.mechanics.cell_cell_adhesion_strength <<std::endl;
	// std::cout<<"cell repulsion = "<<pCell->phenotype.mechanics.cell_cell_repulsion_strength <<std::endl;
	// std::cout<<"cell persistence time ="<<pCell->phenotype.motility.persistence_time <<std::endl;
	// std::cout<<"cell transition rates = "<<phenotype.death.rates[apoptosis_index] <<std::endl;
	// std::cout<<"cell death rates = "<<phenotype.death.rates[necrosis_index] <<std::endl;

	// END New speed update 

	return; 
}

void fibroblast_ECM_informed_motility_update_w_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
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
	static int inflam_sig_index = microenvironment.find_density_index( "inflammatory_signal" ); 

	static int link_anisotropy_and_bias_index = pCell->custom_data.find_variable_index( "link_anisotropy_and_bias" );
	if (link_anisotropy_and_bias_index < 0) 
    {
        std::cout << "        static int link_anisotropy_and_bias_index = " <<link_anisotropy_and_bias_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
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
	std::vector<double> chemotaxis_grad = pCell->nearest_gradient(inflam_sig_index);

	// std::cout<<"D chemo"<<chemotaxis_grad<<std::endl;

	normalize( &chemotaxis_grad ); 

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)

	// New bias - bias such that the agents can more closely follow the gradient IF the written signals are stronger. 

	std::vector<double> d_motility;

	// d_motility = {0,0,0};
	
	if (pCell->custom_data[link_anisotropy_and_bias_index] < 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 0
	{
		d_motility = (1-a) * d_random + a * chemotaxis_grad;
		// std::cout<<" I am coupled"<<std::endl;
	}

	else if (pCell->custom_data[link_anisotropy_and_bias_index] > 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 1
	{
		// NON-ECM linked way to signal 
		d_motility = (1-pCell->custom_data[chemotaxis_bias_index])*d_random + pCell->custom_data[chemotaxis_bias_index]*chemotaxis_grad;
		// std::cout<<" I am UNcoupled"<<std::endl;
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
	std::vector<double> d_perp = d_motility - dot_product_ext(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product_ext( d_motility , d_perp ); 
	double c_2 = dot_product_ext( d_motility, f ); 

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

	// std::cout<<"migration speed(l1236) = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	// double magnitude = norm( phenotype.motility.motility_vector);	

	// std::cout<<"Magnitutude of motility vector is "<< magnitude<<std::endl;

	// if(magnitude > 0.00000001)
	// {
	// 	std::cout<<"Cell is moving!!!!"<<std::endl;
	// }

	/****************************************END new migration direction update****************************************/


	// /*********************************************Begin speed update***************************************************/
	
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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	}

	// else //if (ECM_density >= rho_high)
	// {
	// 	pCell->phenotype.motility.migration_speed = 0.0;
	// }

	// std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// std::cout<<"cell adhesion = "<<pCell->phenotype.mechanics.cell_cell_adhesion_strength <<std::endl;
	// std::cout<<"cell repulsion = "<<pCell->phenotype.mechanics.cell_cell_repulsion_strength <<std::endl;
	// std::cout<<"cell persistence time ="<<pCell->phenotype.motility.persistence_time <<std::endl;
	// std::cout<<"cell transition rates = "<<phenotype.death.rates[apoptosis_index] <<std::endl;
	// std::cout<<"cell death rates = "<<phenotype.death.rates[necrosis_index] <<std::endl;

	// END New speed update 

	return; 
}

void ECM_informed_motility_update_w_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout<<"cell speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;

	
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
	std::vector<double> d_perp = d_motility - dot_product_ext(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product_ext( d_motility , d_perp ); 
	double c_2 = dot_product_ext( d_motility, f ); 

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

		pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
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


void update_adhesion( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Find correct items
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	double anisotropy = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy;
	int adhesion_change_threshold_index = pCell->custom_data.find_variable_index( "adhesion_change_threshold");
	if (adhesion_change_threshold_index < 0) 
    {
        std::cout << "        static int adhesion_change_threshold_index = " <<adhesion_change_threshold_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
	
	if(anisotropy>pCell->custom_data[adhesion_change_threshold_index])
	{pCell->phenotype.mechanics.cell_cell_adhesion_strength = 5.0;}
	else
	{pCell->phenotype.mechanics.cell_cell_adhesion_strength = 10.0;}
}

void custom_cancer_cell_ECM_remodeling_and_adhesion_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	ECM_remodeling_function( pCell, phenotype, dt );

	if (parameters.bools("reduce_adhesion_on_groomed_matrix")==true)
	{
		update_adhesion( pCell, phenotype, dt );
	}

	return;
}

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
	if (Cell_ECM_production_rate_index < 0) 
    {
        std::cout << "        static int Cell_ECM_production_rate_index = " <<Cell_ECM_production_rate_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
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
		
		ddotf = dot_product_ext(ECM_orientation, norm_cell_motility);
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
