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

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>
#include <omp.h>
#include <fstream>

#include "./core/PhysiCell.h"
#include "./modules/PhysiCell_standard_modules.h" 
#include "./custom_modules/cell_ECM_interactions.h"

// custom user modules 

#include "./custom_modules/fibrosis.h" 
// #include "./custom_modules/AMIGOS-invasion.h" 
// #include "./custom_modules/ECM.h"
	
using namespace BioFVM;
using namespace PhysiCell;


// set number of threads for OpenMP (parallel computing)

// void ecm_update(void); // I think this from an old implentation. Noted and commented out 09.01.20


int main( int argc, char* argv[] )
{
	// std::cout<<"test"<<std::endl;
	// load and parse settings file(s)
 
	bool XML_status = false; 
	char copy_command [1024]; 
	if( argc > 1 )
	{
		XML_status = load_PhysiCell_config_file( argv[1] ); 
		sprintf( copy_command , "cp %s %s" , argv[1] , PhysiCell_settings.folder.c_str() ); 
	}
	else
	{
		XML_status = load_PhysiCell_config_file( "./config/PhysiCell_settings.xml" );
		sprintf( copy_command , "cp ./config/PhysiCell_settings.xml %s" , PhysiCell_settings.folder.c_str() ); 
	}
	if( !XML_status )
	{ exit(-1); }

	// copy config file to output directry 
	system( copy_command ); 

	// OpenMP setup
	omp_set_num_threads(PhysiCell_settings.omp_num_threads);
	
	// PNRG setup 
	// SeedRandom(0);
	if( parameters.ints("unit_test_setup") == 1) 
	{SeedRandom(0);}
	
	// time setup 
	std::string time_units = "min"; 
	
	double t_max = PhysiCell_settings.max_time;  // 1 days

	/* Microenvironment setup */ 
	setup_microenvironment();
	setup_extracellular_matrix(); // NEW LINE!!!!

	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	create_cell_types();
	setup_tissue();
	// ECM_setup(microenvironment.number_of_voxels());
 
	/* Users typically start modifying here. START USERMODS */ 
	
	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 
	
	char filename[1024];
	sprintf( filename , "%s/initial" , PhysiCell_settings.folder.c_str() ); 
	// save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , PhysiCell_globals.current_time ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = AMIGOS_invasion_coloring_function;
	
	sprintf( filename , "%s/initial.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot_custom( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, parameters.strings("visual_guideline_pattern") );
	
	sprintf( filename , "%s/legend.svg" , PhysiCell_settings.folder.c_str() );
	create_plot_legend( filename , cell_coloring_function ); 

	display_citations(); 
	
	run_biotransport( parameters.doubles("duration_of_uE_conditioning") ); 

	if(parameters.bools("freeze_uE_profile")==true)
	{
		alter_cell_uptake_secretion_saturation();
	}

	// if(parameters.ints("unit_test_setup")==1 || parameters.ints("discrete_ECM_remodeling") == 1)
	// {
	set_cell_motility_vectors(); // Required for instant writing and unit test. To make all simulations have similar initial conditions, requiring it for all simulations at this time 05.27.22
	// }
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file;

	//variables for March Project
	double reset_Cells_interval = 1960.0; // for a 1000 by 1000 um computational domain (use 980.0 for speed = 1.0. Doulbe it for speed equals 0.5. )
	bool enable_cell_resets = true;

	// main loop 
	if( PhysiCell_settings.enable_legacy_saves == true )
	{	
		sprintf( filename , "%s/simulation_report.txt" , PhysiCell_settings.folder.c_str() ); 
		
		report_file.open(filename); 	// create the data log file 
		report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;
	}
	
	try 
	{	
		while( PhysiCell_globals.current_time < t_max + 0.1*diffusion_dt )
		{
			// save data if it's time. 
            if(  fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_full_save_time ) < 0.01 * diffusion_dt )
            {
				
				display_simulation_status( std::cout ); 
				if( PhysiCell_settings.enable_legacy_saves == true )
				{	
					log_output( PhysiCell_globals.current_time , PhysiCell_globals.full_output_index, microenvironment, report_file);
				}
				
				if( PhysiCell_settings.enable_full_saves == true )
				{	
					sprintf( filename , "%s/output%08u" , PhysiCell_settings.folder.c_str(),  PhysiCell_globals.full_output_index ); 
					
					save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
				}
				if( parameters.bools("enable_ecm_outputs") == true)
				{
					sprintf( filename , "%s/output%08u_ECM.mat" , PhysiCell_settings.folder.c_str(), PhysiCell_globals.full_output_index);
					write_ECM_Data_matlab( filename );
				}
				PhysiCell_globals.full_output_index++; 
				PhysiCell_globals.next_full_save_time += PhysiCell_settings.full_save_interval;
			}
			
			// save SVG plot if it's time
			if( fabs( PhysiCell_globals.current_time - PhysiCell_globals.next_SVG_save_time  ) < 0.01 * diffusion_dt )
			{
				if( PhysiCell_settings.enable_SVG_saves == true )
				{	
					sprintf( filename , "%s/snapshot%08u.svg" , PhysiCell_settings.folder.c_str() , PhysiCell_globals.SVG_output_index ); 
					SVG_plot_custom( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, parameters.strings("visual_guideline_pattern"));
					
					PhysiCell_globals.SVG_output_index++; 
					PhysiCell_globals.next_SVG_save_time  += PhysiCell_settings.SVG_save_interval;
				}
			}

			// Uncomment to run march test
 			if( fabs( PhysiCell_globals.current_time - reset_Cells_interval  ) <  0.1 * diffusion_dt && parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)	
			{
				if (enable_cell_resets == true )
				{
					reset_cell_position();
					reset_Cells_interval += 1960.0; // for a 1000 by 1000 um computational domain (use 980.0 for speed = 1.0. Doulbe it for speed equals 0.5. )
				}
				
			}
		
			// update the microenvironment

			// if(parameters.bools("freeze_uE_profile")==true)
			
			microenvironment.simulate_diffusion_decay( diffusion_dt );

			// copy_ECM_data_to_BioFVM();
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells( PhysiCell_globals.current_time );
			
			//add ECM update here!
            
            // This changes the cell speed and bias as based on the ECM. It is the funciton that makes teh cells "see" the ECM and react to it with changes in their dynamics.
            // In this current LS18 implementation, that means that only follower cells will see the ECM.
            
            // May need something that specifics that only followers do this (fixed with test on cell type). Could perhaps put that into the custom stuff. Would need something similar for the ECM realignment. Would like to move this out of diffusion loop so we can specify how frequently it updates.
            
            // cell_update_from_ecm();
            
			PhysiCell_globals.current_time += diffusion_dt;
			
			if( PhysiCell_settings.enable_legacy_saves == true )
			{			
			log_output(PhysiCell_globals.current_time, PhysiCell_globals.full_output_index, microenvironment, report_file);
			report_file.close();
			}
		}
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	sprintf( filename , "%s/final" , PhysiCell_settings.folder.c_str() ); 
	save_PhysiCell_to_MultiCellDS_v2( filename , microenvironment , PhysiCell_globals.current_time ); 
	
	sprintf( filename , "%s/final.svg" , PhysiCell_settings.folder.c_str() ); 
	SVG_plot_custom( filename , microenvironment, 0.0 , PhysiCell_globals.current_time, cell_coloring_function, parameters.strings("visual_guideline_pattern") );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}