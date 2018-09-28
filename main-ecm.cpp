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

// custom user modules 

#include "./custom_modules/AMIGOS-invasion.h" 
#include "./custom_modules/ECM.cpp"
	
using namespace BioFVM;
using namespace PhysiCell;


// set number of threads for OpenMP (parallel computing)
int omp_num_threads = 12; // set this to # of CPU cores x 2 (for hyperthreading)
void ecm_update(void);

int main( int argc, char* argv[] )
{
	// OpenMP setup
	omp_set_num_threads(omp_num_threads);
	
	// PNRG setup 
	SeedRandom(); 
	
	// time setup 
	std::string time_units = "min"; 
	double t = 0.0; // current simulation time 
	
	double t_output_interval = 10; // output once per hour WHY ISN'T THE CONFIG FILE WORKING??
	double t_max = 60*24*2;  // 1 days 
	double t_next_output_time = t; 
	int output_index = 0; // used for creating unique output filenames 

	/* Microenvironment setup */ 
	
	setup_microenvironment();

	/* PhysiCell setup */ 
 	
	// set mechanics voxel size, and match the data structure to BioFVM
	double mechanics_voxel_size = 30; 
	Cell_Container* cell_container = create_cell_container_for_microenvironment( microenvironment, mechanics_voxel_size );
	
	/*tumor_radius for setup_tissue and ECM_setup*/
	double tumor_radius = 175.0;
	
	create_cell_types();
	setup_tissue(tumor_radius);
	ECM_setup(microenvironment.number_of_voxels(), tumor_radius);
	
	/* Users typically start modifying here. START USERMODS */ 
	
	/* Users typically stop modifying here. END USERMODS */ 
	
	// set MultiCellDS save options 

	set_save_biofvm_mesh_as_matlab( true ); 
	set_save_biofvm_data_as_matlab( true ); 
	set_save_biofvm_cell_data( true ); 
	set_save_biofvm_cell_data_as_custom_matlab( true );
	
	// save a simulation snapshot 

	save_PhysiCell_to_MultiCellDS_xml_pugi( "initial" , microenvironment , t ); 
	
	// save a quick SVG cross section through z = 0, after setting its 
	// length bar to 200 microns 

	PhysiCell_SVG_options.length_bar = 200; 

	// for simplicity, set a pathology coloring function 
	
	std::vector<std::string> (*cell_coloring_function)(Cell*) = AMIGOS_invasion_coloring_function;
	
	SVG_plot( "initial.svg" , microenvironment, 0.0 , t, cell_coloring_function );
	
// Is this to initialize the uE with the leader and follower signals?
    
	run_biotransport( 5.0 ); 
	
	// set the performance timers 

	BioFVM::RUNTIME_TIC();
	BioFVM::TIC();
	
	std::ofstream report_file ("simulation_report.txt"); 	// create the data log file 
	report_file<<"simulated time\tnum cells\tnum division\tnum death\twall time"<<std::endl;

	// main loop 
	
	try 
	{	
		while( t < t_max + 0.1*diffusion_dt )
		{
			// save data if it's time. 
            if(  fabs( t - t_next_output_time ) < 0.01 * diffusion_dt )
            {
                log_output(t, output_index, microenvironment, report_file);
                
                char filename[1024];
                sprintf( filename , "Output/output%08u" , output_index );
                
//                                sprintf( filename , "output%08u_ECM.mat" , output_index);
                
                save_PhysiCell_to_MultiCellDS_xml_pugi( filename , microenvironment , t );
                
                sprintf( filename , "Output/output%08u_ECM.mat" , output_index);
                
                
                write_ECM_Data_matlab( filename );
                
                sprintf( filename , "SVG/snapshot%08u.svg" , output_index );
                
                SVG_plot( filename , microenvironment, 0.0 , t, cell_coloring_function );
                
                
                output_index++;
                t_next_output_time += t_output_interval;
            }
			// update the microenvironment
			microenvironment.simulate_diffusion_decay( diffusion_dt );
			if( default_microenvironment_options.calculate_gradients )
			{ microenvironment.compute_all_gradient_vectors(); }
			
			// run PhysiCell 
			((Cell_Container *)microenvironment.agent_container)->update_all_cells(t);
			
			//add ECM update here!
            
            // This changes the cell speed and bias as based on the ECM. It is the funciton that makes teh cells "see" the ECM and react to it with changes in their dynamics.
            // In this current LS18 implementation, that means that only follower cells will see the ECM.
            
            // Need somethign that specifics that only followers do this. Maybe put that into the custom stuff ... Not sure how to do that. Will need somethign similar for the ECM realignment.
            
            cell_update_from_ecm();
			
			t += diffusion_dt; 
		}
		log_output(t, output_index, microenvironment, report_file);
		report_file.close();
	}
	catch( const std::exception& e )
	{ // reference to the base of a polymorphic object
		std::cout << e.what(); // information from length_error printed
	}
	
	// save a final simulation snapshot 
	
	save_PhysiCell_to_MultiCellDS_xml_pugi( "final" , microenvironment , t ); 
	SVG_plot( "final.svg" , microenvironment, 0.0 , t, cell_coloring_function );
	
	// timer 
	
	std::cout << std::endl << "Total simulation runtime: " << std::endl; 
	BioFVM::display_stopwatch_value( std::cout , BioFVM::runtime_stopwatch_value() ); 

	return 0; 
}

