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
#include "./ECM.h"
Cell_Definition leader_cell; 
Cell_Definition follower_cell; 

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
	
	cell_defaults.parameters.o2_proliferation_saturation = 38.0;
	cell_defaults.parameters.o2_reference = 38.0;
	
	
	// set default uptake and secretion 
	// oxygen 
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0; 
	cell_defaults.phenotype.secretion.uptake_rates[0] = 10; 
	cell_defaults.phenotype.secretion.saturation_densities[0] = 38; 

//    cell_defaults.phenotype.secretion.secretion_rates[1] = 0;
//    cell_defaults.phenotype.secretion.uptake_rates[1] = 0;
//    cell_defaults.phenotype.secretion.saturation_densities[1] = 1;
//
//    cell_defaults.phenotype.secretion.secretion_rates[2] = 0;
//    cell_defaults.phenotype.secretion.uptake_rates[2] = 0;
//    cell_defaults.phenotype.secretion.saturation_densities[2] = 1;


	// set the default cell type to no phenotype updates 
	
//    cell_defaults.functions.update_phenotype = switching_phenotype_model;
	
	cell_defaults.name = "cancer cell"; 
	cell_defaults.type = 0; 
	
	// set default motility parameters (even for when off)
	
	cell_defaults.phenotype.motility.is_motile = true;
	cell_defaults.phenotype.motility.persistence_time = 15.0;
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles("default_cell_speed");
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	cell_defaults.phenotype.motility.migration_bias = 0.90;
	
	// add custom data 
	
	cell_defaults.custom_data.add_variable( "hypoxic switch value" , "mmHg", 38 );
	
	// leader cells 
	
	leader_cell = cell_defaults;
	leader_cell.name = "leader cell"; 
	leader_cell.type = 1; 

//    For SIAM LS18 Motility presentation - eliminating leader/follower signal
    
	// 10% proliferation 
    leader_cell.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) *= 0.10;
    
//    For SIAM LS18 Motility presentation - eliminating leader/follower signal
	
	// turn on motility 
	leader_cell.phenotype.motility.is_motile = parameters.bools("leader_motility_mode"); 
	
    leader_cell.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("leader_adhesion");
	std::cout<<leader_cell.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
    
	leader_cell.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("leader_repulsion");
    std::cout<<leader_cell.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
//    leader_cell.phenotype.secretion.secretion_rates[1] = 50; // leader signal
    
//    For SIAM LS18 Motility presentation - eliminating leader/follower signal

    // modify ECM
    
//    leader_cell.functions.custom_cell_rule = ecm_update_from_cell; // this meas that ... only leaders have this, right? So, it is only accessed when leaders are updated, so this shoudl autoamtically make it so that followers can't modfiy ECM ...
	
	// set functions
	
	leader_cell.functions.update_migration_bias = chemotaxis_oxygen; 
	
    leader_cell.functions.update_phenotype = leader_cell_phenotype_model;
	
	// follower cells

	follower_cell = cell_defaults;
	follower_cell.name = "follower cell"; 
	follower_cell.type = 2;
    
//    follower_cell.functions.update_migration_bias = chemotaxis_oxygen;
    
    follower_cell.functions.update_phenotype = follower_cell_phenotype_model;
    std::cout<<follower_cell.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;

	follower_cell.phenotype.mechanics.cell_cell_adhesion_strength = parameters.doubles("follower_adhesion");
    std::cout<<follower_cell.phenotype.mechanics.cell_cell_adhesion_strength<<std::endl;
    std::cout<<follower_cell.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
	follower_cell.phenotype.mechanics.cell_cell_repulsion_strength = parameters.doubles("follower_repulsion");
    std::cout<<follower_cell.phenotype.mechanics.cell_cell_repulsion_strength<<std::endl;
	follower_cell.phenotype.motility.is_motile = parameters.bools("follower_motility_mode");
    
//    follower_cell.phenotype.secretion.secretion_rates[2] = 50; // follower signal
    
//    For SIAM LS18 Motility presentation - eliminating leader/follower signal/differencse in adhesion/etc
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters

	default_microenvironment_options.X_range = {-1500, 1500};
	default_microenvironment_options.Y_range = {-1500, 1500};
	default_microenvironment_options.simulate_2D = true; 
	
	// no gradients needed for this example 
	
	default_microenvironment_options.calculate_gradients = true; 

//    For SIAM LS18 Motility presentation - eliminating leader/follower signal
    
	// 50 micron length scale 
    microenvironment.add_density( "leader signal", "dimensionless", 1e5 , 1 );
    microenvironment.add_density( "follower signal", "dimensionless", 1e5 , 1 );

//    For SIAM LS18 Motility presentation - eliminating leader/follower signal
    
	// let BioFVM use oxygen as the default 
	
	default_microenvironment_options.use_oxygen_as_first_field = true; 
	
	// set Dirichlet conditions 
	
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	default_microenvironment_options.Dirichlet_condition_vector[0] = 38; // normoxic conditions
    
//    For SIAM LS18 Motility presentation - eliminating leader/follower signal
    
//    default_microenvironment_options.Dirichlet_condition_vector[1] = 0; // normoxic conditions
//    default_microenvironment_options.Dirichlet_condition_vector[2] = 0; // normoxic conditions

//    For SIAM LS18 Motility presentation - eliminating leader/follower signal
    
	initialize_microenvironment(); 
	
	// run to get a decent startin conditoin 
	
	// now, let's set the leader signla to 1, so we don't hvae early swiching 
	/*
	for( int i=0 ; i < microenvironment.number_of_voxels() ; i++ )
	{
		microenvironment.density_vector(i)[1] = 1.0; 
	}
	*/

	return; 
}	

void ECM_setup(double numvox)
{
	//Save initial parameters from XML file once so we don't have to call fxn repeatedly
	double initial_density = parameters.doubles("initial_density_ecm");
	double tumor_radius = parameters.doubles("tumor_radius");
	
    ecm.sync_to_BioFVM();
    ecm.ecm_data.resize(numvox);
//    std::cout<<"Hi! 1"<<std::endl;
//    double ECM_radius = default_microenvironment_options.X_range[1]
    
//    std::cout<<ecm.mesh.voxels.size()<<std::endl;
//    double temp = ecm.mesh.voxels[numvox-1].center.at(0);
//    std::cout<<temp<<std::endl;
	
	for (int i = 0; i<numvox-1; i++)
    {
//  /*      std::cout<<"Hi! 2"<<std::endl;
//          std::cout<<i<<std::endl;
        
//	This block of code is to randomly orient the ecm fibers using vector randomization from BioFVM_vector.cpp, line 262
//  Pick a random angle from 0 to 2pi and then set components equal to sin(theta) and cos(theta)
		double theta = 6.2831853071795864769252867665590 * uniform_random(); 
		ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
		ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
		ecm.ecm_data[i].ECM_orientation[2] = 0.0;
		
		double buffer_region = 20;
		
		ecm.ecm_data[i].density = initial_density;
		//if sqrt(x^2 + y^2) > (tumor rad + buf region)
		//get coord from ecm.mesh.voxels[i].center[0] -x ecm.mesh.voxels[i].center[1]
		
		//This block of code is used to create an initial density field that is different what the
		//default ECM constructor instantiates
		/*double x = ecm.mesh.voxels[i].center[0];
		double y = ecm.mesh.voxels[i].center[1];
		if((sqrt((x*x) + (y*y)) > (tumor_radius + 20)))
		{
			ecm.ecm_data[i].density = 1.0;
		}*/
		
	//This block of code is used for orienting the ecm fibers radially outward from the origin
        /*double epsilon = 1E-6;
        double ECM_radius = sqrt(ecm.mesh.voxels[i].center.at(0) * ecm.mesh.voxels[i].center.at(0)
                                 +ecm.mesh.voxels[i].center.at(1) * ecm.mesh.voxels[i].center.at(1)
                                 + ecm.mesh.voxels[i].center.at(2) * ecm.mesh.voxels[i].center.at(2));
//        std::cout<<"Hi! 3"<<std::endl;
        ecm.ecm_data[i].ECM_orientation[0] = ecm.mesh.voxels[i].center[0]/(ECM_radius + epsilon);
        ecm.ecm_data[i].ECM_orientation[1] = ecm.mesh.voxels[i].center[1]/(ECM_radius + epsilon);
        ecm.ecm_data[i].ECM_orientation[2] = ecm.mesh.voxels[i].center[2]/(ECM_radius + epsilon);
//        std::cout<<"Hi! 4"<<std::endl;
        normalize(ecm.ecm_data[i].ECM_orientation);
//        std::cout<<"Hi! 5"<<std::endl;*/
        
        
    }
    return;
}

void run_biotransport( double t_max )
{
	std::cout << "working on initial conditions .. " << std::endl; 
	double t = 0.0;
	
	// make sure the associated cell has the correct rate vectors 
	
	for( int i=0 ; i < (*all_cells).size() ; i++ )
	{
		Cell* pCell = (*all_cells)[i];
		
//			pCell->secretion_rates = &secretion_rates; 
//			pCell->uptake_rates = &uptake_rates; 
//			pCell->saturation_densities = &saturation_densities; 
			
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
	// place a cluster of tumor cells at the center 
	
//    pCell->custom_data[0] = NormalRandom( 1.0, 0.33 );
//    if( pCell->custom_data[0] < 0.0 )
//    { pCell->custom_data[0] = 0.0; }
//    if( pCell->custom_data[0] > 2.0 )
//    { pCell->custom_data[0] = .0; }
    
	//Get tumor radius from XML parameters
	double tumor_radius = parameters.doubles("tumor_radius");
	double cell_radius = cell_defaults.phenotype.geometry.radius; 
	double cell_spacing = 0.95 * 2.0 * cell_radius; 
	
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double x_outer = tumor_radius; 
	double y = 0.0;
	
	double leader_cell_fraction = 1.0;
	
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
			{ pCell = create_cell(follower_cell); }
				
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
//                    std::cout<<"created follower cell"<<std::endl;
					pCell->assign_position( -x , -y , 0.0 );
				}
			}
			x += cell_spacing; 
			
		}
		
		y += cell_spacing * sqrt(3.0)/2.0; 
		n++; 
	}
		
	return; 
}


void chemotaxis_oxygen( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int o2_index = microenvironment.find_density_index( "oxygen" ); 
	
	// if attached, biased motility towards director chemoattractant 
	// otherwise, biased motility towards cargo chemoattractant 
	
	phenotype.motility.is_motile = true; 
	phenotype.motility.migration_bias = 0.95;
	phenotype.motility.migration_bias_direction = pCell->nearest_gradient(o2_index);
    std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;
	
	return; 
}

void change_migration_bias_vector_ecm(Cell* pCell , Phenotype& phenotype , double dt )
{
    //change bias

    //    std::vector<double> migration bias direction is the 3-D vector giving the cell's preferred
    //    direction of motility for biased Brownian motion. If the user modies this vector, they must ensure
    //    it is a unit vector :
    //    jjmigration bias directionjj = 1: (14)
    //    54
    //    5. double migration bias (with a value in [0,1]) sets the degree to which cell motility is biased
    //    along migration_bias_direction. If 0, then motion is completely Brownian. If 1, it is completely
    //    deterministc along the bias direction.

    int ecm_index =  pCell->get_current_voxel_index();
    double a = ecm.ecm_data[ecm_index].anisotropy;
    
    std::vector<double> d = pCell->phenotype.motility.motility_vector; ////// Changed to motility_vector instead of bias - so the blending is of the actual velocity instead of just the direction it wants to go in
    //    d = norm(d);
    std::vector<double> f = ecm.ecm_data[ecm_index].ECM_orientation;
    double ddotf = 0.0;
//normalize( &( phenotype.motility.migration_bias_direction ) );
    normalize(d);
    normalize(f);

    //     pCell->phenotype.motility.migration_bias = a;

    //change bias direction
    for( int i=0; i < d.size() ; i++ )
    {
        //         double temp = d[i] * f[i];
        ddotf += d[i] * f[i];
    }

    if(ddotf < 0.0)
    {
        for( int i=0; i< f.size(); i++)
        {
            f[i] *= -1.0;
            //            ecm.ecm_data[ecm_index].ECM_orientation[i] = -1.0 * f[i];
        }
    }

    for( int i=0; i < d.size() ; i++ )
    {
        pCell->phenotype.motility.migration_bias_direction[i] = a * f[i]  + (1.0-a) * d[i];
    }

    normalize(pCell->phenotype.motility.migration_bias_direction);
    
    phenotype.motility.migration_bias = a;

    return;
}

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell , Phenotype& phenotype , double dt ) 
{	
	// o2-based birth and death, nothing more 
	update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
	return; 
}


void follower_cell_phenotype_model0( Cell* pCell , Phenotype& phenotype , double dt )
{
    // o2-based birth and death, nothing more
    tumor_cell_phenotype_with_oncoprotein(pCell,phenotype,dt);

    return;
}


/* old */ 


std::vector<std::string> ECM_anisotropy_coloring_function( Cell* pCell)
{
    std::vector< std::string > output( 4, "black" );
    char color [1024];
    int ecm_index = pCell->get_current_voxel_index();
//    std::cout<<ecm_index<<std::endl;
    double anisotropy = ecm.ecm_data[ecm_index].anisotropy;
    sprintf( color , "rgb(%d,%d,%d)" , int(anisotropy*255.0), int(anisotropy* 255.0), int(255-anisotropy*255));
//    std::cout<<color<<std::endl;
//    rgb(0,0,255)
    output[0] = color;
    output[2] = color;
    return output;
}

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

    if( pCell->type == 2 )
    {
        output[0] = "yellow";
        output[2] = "yellow";
        return output;
    }
	
	// live followers cells are red,
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

    phenotype.motility.is_motile = true;
    
//    if( pO2 > pCell->custom_data[hypoxic_i] )
//    {
//        // proliferate (don't overwrite)
////        phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) =
////            10.0 * pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
//        // turn off motility
//        phenotype.motility.is_motile = true ;
//    }
//    else
//    {
//        // don't proliferate,
////        phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
//        // turn on motility
//        phenotype.motility.is_motile = true;
//    }
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
    
    phenotype.motility.is_motile = true;
//    phenotype.motility.migration_bias = 1.0;
    
    //    if( pO2 > pCell->custom_data[hypoxic_i] )
    //    {
    //        // proliferate (don't overwrite)
    ////        phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) =
    ////            10.0 * pCell->parameters.pReference_live_phenotype->cycle.data.transition_rate(start_phase_index,end_phase_index);
    //        // turn off motility
    //        phenotype.motility.is_motile = true ;
    //    }
    //    else
    //    {
    //        // don't proliferate,
    ////        phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) *= 0.1;
    //        // turn on motility
    //        phenotype.motility.is_motile = true;
    //    }
    return;
}

// commented out by JPM 11.04.18 to eliminate chance of seg fault from currnetly non-exsistent fields.

//void switching_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt )
//{
//    static int hypoxic_i = pCell->custom_data.find_variable_index( "hypoxic switch value" );
//    static int oxygen_i = pCell->get_microenvironment()->find_density_index( "oxygen" );
//
//
//    int cycle_start_index = live.find_phase_index( PhysiCell_constants::live );
//    int cycle_end_index = live.find_phase_index( PhysiCell_constants::live );
//
//    double pO2 = (pCell->nearest_density_vector())[oxygen_i]; // PhysiCell_constants::oxygen_index];
//
//    // set death and birth
//    update_cell_and_death_parameters_O2_based(pCell,phenotype,dt);
//
//    // if a leader and now happy, switch to follower
//    if( pO2 > pCell->custom_data[hypoxic_i] && pCell->type == 1 )
//    {
//        pCell->convert_to_cell_definition( follower_cell );
//
//    }
//
//    // if a folower and now unhappy, switch to leader, if not supressed
//
//    if( pO2 < pCell->custom_data[hypoxic_i] && pCell->type == 2 )
//    {
//        if( (pCell->nearest_density_vector())[1] < 0.1 )
//        {
//         pCell->convert_to_cell_definition( leader_cell );
//        }
//    }
//
//    return;
//}

void ecm_update_from_cell(Cell* pCell , Phenotype& phenotype , double dt) // NOTE - not currently supporting ECM density increasing or anisotropy decreasing!!! 03.30.18
{
//    Cell* pCell = (*all_cells)[i];
    int ecm_index = pCell->get_current_voxel_index();
    
    // Type 2 are the followers. This prohibits them from modifying the ECM - the will be affected by its direction and ani, but can't change those properties.
    
//    if( pCell->type == 2 )
//    {
//        return;
//    }
    
    // Cell-ECM density interaction
    
    double density = ecm.ecm_data[ecm_index].density;
    double r = 1.0;
    
    ecm.ecm_data[ecm_index].density = density + r * dt  * (0.5 - density);
    
    // END Cell-ECM density interaction
    
    // Cell-ECM Fiber realingment

    std::vector<double> ECM_orientation = ecm.ecm_data[ecm_index].ECM_orientation;
    
    double motility_vector_norm = norm( phenotype.motility.motility_vector );
    if( motility_vector_norm < 1e-12 )
    { return; }
    else
        
        // FIX ME!!!!
        
    { phenotype.motility.motility_vector /= motility_vector_norm; }
//    std::vector<double> d = normalize(pCell->phenotype.motility.motility_vector);
    double anisotropy = ecm.ecm_data[ecm_index].anisotropy;
    double migration_speed = pCell->phenotype.motility.migration_speed;
//    std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;

    double r_0 = 1/1.0*migration_speed; // min-1 // NOTE!!! on 08.06.18 run - this wasn't multiplied by migration_speed!!! should be the same but worth noting!!!!

    double r_realignment = r_0 * (1-anisotropy);
    
    // START HERE! Double check math as well as if time step is working!

    /*for( int i=0; i < phenotype.motility.motility_vector.size() ; i++ )
    {
        double temp = ECM_orientation[i] + dt * (r_realignment * (ECM_orientation[i] - phenotype.motility.motility_vector[i]));
//        std::cout<<temp<<std::endl<<std::endl;
        ecm.ecm_data[ecm_index].ECM_orientation[i] = temp;
		
		/*Make sure each vector component of the ecm orientation is going the same direction as its corresponding component*
		 *of the cell's motility vector
		if(ecm.ecm_data[ecm_index].ECM_orientation[i] * phenotype.motility.motility_vector[i] < 0.0)
	       ecm.ecm_data[ecm_index].ECM_orientation[i] *= -1.0;
    }*/

	double ddotf;
	std::vector<double> temp;
	temp.resize(3,0.0);
	for(int i = 0; i < 3; i++)
	{
		temp[i] = ECM_orientation[i] * phenotype.motility.motility_vector[i];
	}
	ddotf = temp[1] + temp[2] + temp[3];
	
	if(ddotf < 0)
	{
		for(int i = 0; i < 3; i++)
		{
		   ECM_orientation[i] *= -1.0;
		}
	}
	
	std::vector<double> f_minus_d;
	f_minus_d.resize(3,0.0);
	for(int i = 0; i < 3; i++)
	{
		f_minus_d[i] = ECM_orientation[i] - phenotype.motility.motility_vector[i];
		ecm.ecm_data[ecm_index].ECM_orientation[i] -= dt * r_realignment * f_minus_d[i];
	}
	
	
    normalize(&(ecm.ecm_data[ecm_index].ECM_orientation));


    // End Cell-ECM Fiber realingment
    
    // Cell-ECM Anisotrophy Modification
    
    double r_a0 = 1.0/100.0; // min-1 - changes on same time scale as fiber realignment????
    double r_anisotropy = r_a0 * migration_speed; // What is a typical cell speed????
//    std::cout<<migration_speed<<std::endl;
    ecm.ecm_data[ecm_index].anisotropy = anisotropy + r_anisotropy * dt  * (1- anisotropy);
    
    // END Cell-ECM Anisotropy Modification
    
    return;
    
    
    
}
