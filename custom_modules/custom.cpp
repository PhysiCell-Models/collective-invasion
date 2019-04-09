/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
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

#include "./custom.h"

// declare cell definitions here 

Cell_Definition motile_cell; 

std::vector< std::vector<double> > ECM_fiber_alignment; 

void create_cell_types( void )
{
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 

	// housekeeping 
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	cell_defaults.type = 0; 
	cell_defaults.name = "cell"; 
	
	// set default cell cycle model 
	cell_defaults.functions.cycle_model = live; // let's keep it simple 
	
	// set default_cell_functions; 
	cell_defaults.functions.update_phenotype = NULL; // no need for fancy stuff here 
	cell_defaults.functions.update_migration_bias = ECM_motility_aligned; 
	
	// needed for a 2-D simulation: 
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	// make sure the defaults are self-consistent. 
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// first find indices for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 
	int nCycle_start = live.find_phase_index( PhysiCell_constants::live );
	int nCycle_end = live.find_phase_index( PhysiCell_constants::live );

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 
	
	// add custom data here, if any 
	cell_defaults.custom_data.add_variable( "max speed", "micron/min" , 
		parameters.doubles( "max_migration_speed") ); 
	
	// enable motility 
	cell_defaults.phenotype.motility.is_motile = true; 
	cell_defaults.phenotype.motility.persistence_time = parameters.doubles( "migration_persistence_time" ); // 15.0; 
	cell_defaults.phenotype.motility.migration_speed = parameters.doubles( "max_migration_speed" ); // 0.25 micron/minute 
	cell_defaults.phenotype.motility.migration_bias = 0.0;// completely random 
	
	// Set birth and death rates to zero 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	cell_defaults.phenotype.cycle.data.transition_rate(nCycle_start,nCycle_end) = 0.0; 
	 
	return; 
}

void setup_microenvironment( void )
{
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
	microenvironment.add_density( "ECM", "dimensionless", 0.0 , 0.0 ); 
	microenvironment.add_density( "ECM anisotropy" , "dimensionless" , 0.0 , 0.0 ); 
	
	// calculate gradients, just in case 
	default_microenvironment_options.calculate_gradients = true; 
	
	// set Dirichlet conditions 
	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector = { 38.0 , 0.5 , 1 };  // 5% o2 , half max ECM , isotropic  
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// initialize BioFVM 
	initialize_microenvironment(); 	
	
	// set up ECM profile
	int nECM = microenvironment.find_density_index( "ECM" ); 
	int nECM_A = microenvironment.find_density_index( "ECM anisotropy" ); 
	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		std::vector<double> position = microenvironment.mesh.voxels[n].center; 
		if( fabs( position[0] ) > 200 || fabs( position[1] ) > 200 )
		{
			microenvironment(n)[nECM] = 0.0; 
			microenvironment(n)[nECM_A] = 1.0; 
		}
	}
	
	// set up ECM alignment 
	std::vector<double> fiber_direction = { 1.0 , 0.0, 0.0 }; 
	ECM_fiber_alignment.resize( microenvironment.mesh.voxels.size() , fiber_direction );  
	
	for( int n = 0; n < microenvironment.mesh.voxels.size() ; n++ )
	{
		std::vector<double> position = microenvironment.mesh.voxels[n].center; 
		normalize( position ); 
		ECM_fiber_alignment[n] = { -position[1],position[0],0}; // position; 
	}
	
	return; 
}

void setup_tissue( void )
{
	Cell* pC;
	
	for( int n = 0 ; n < 200 ; n++ )
	{
		pC = create_cell(); 
		pC->assign_position( -450 + 900*UniformRandom() , -450 + 900*UniformRandom() , 0.0 );
	}
		
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	// start with flow cytometry coloring 
	
	std::vector<std::string> output( 4 , "red" ); 
	
	return output; 
}

void ECM_motility( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find location of variables 
	static int nECM = microenvironment.find_density_index( "ECM" ); 
	static int nMaxSpeed = pCell->custom_data.find_variable_index( "max speed" ); 
	
	// sample ECM 
	double ECM = pCell->nearest_density_vector()[nECM]; 
	phenotype.motility.migration_speed = pCell->custom_data[nMaxSpeed]*ECM*(1-ECM)*4.0; 
	
	return; 
}

double dot_product( const std::vector<double>& v , const std::vector<double>& w )
{
 double out = 0.0; 
 for( unsigned int i=0 ; i < v.size() ; i++ )
 { out += ( v[i] * w[i] ); }
 return out; 
}

void ECM_motility_aligned( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find location of variables 
	static int nECM = microenvironment.find_density_index( "ECM" ); 
	static int nMaxSpeed = pCell->custom_data.find_variable_index( "max speed" ); 
	static int nECM_A = microenvironment.find_density_index( "ECM anisotropy" ); 
	
	// sample ECM 
	double ECM = pCell->nearest_density_vector()[nECM]; 
	double a = pCell->nearest_density_vector()[nECM_A]; 
	
	int n = pCell->get_current_voxel_index();

	
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_mot = { cos(angle) , sin(angle) , 0.0 }; 
	
	// fiber direction
	std::vector<double> f = ECM_fiber_alignment[n]; 
	
	// part of d_mot that is perpendicular to f; 
	std::vector<double> d_perp = d_mot - dot_product(d_mot,f)*f; 
	normalize( d_perp ); 
	
	double c_1 = dot_product( d_mot , d_perp ); 
	double c_2 = dot_product( d_mot, f ); 
	
	double sensitivity = 1.0; 
	double theta = a*sensitivity; 
	
	phenotype.motility.migration_bias_direction = (1.0-theta)*c_1*d_perp + c_2*f;
	normalize( phenotype.motility.migration_bias_direction ); 
	phenotype.motility.migration_bias = 1.0; 

	phenotype.motility.migration_speed = pCell->custom_data[nMaxSpeed]*ECM*(1-ECM)*4.0; 
	return; 
}

void ECM_motility_aligned_faster( Cell* pCell, Phenotype& phenotype, double dt )
{
	// find location of variables 
	static int nECM = microenvironment.find_density_index( "ECM" ); 
	static int nMaxSpeed = pCell->custom_data.find_variable_index( "max speed" ); 
	static int nECM_A = microenvironment.find_density_index( "ECM anisotropy" ); 
	
	// sample ECM 
	double ECM = pCell->nearest_density_vector()[nECM]; 
	double a = pCell->nearest_density_vector()[nECM_A]; 
	
	int n = pCell->get_current_voxel_index();

	
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_mot = { cos(angle) , sin(angle) , 0.0 }; 

	std::vector<double>* pDmot = &( phenotype.motility.migration_bias_direction ); 
	std::vector<double>* pF = &( ECM_fiber_alignment[n] ); 
	
	// part of d_mot that is perpendicular to f; 
	double c2 = dot_product( *pDmot , *pF ); 

	std::vector<double> d_perp = *pF;
	d_perp *= c2; // (dMot.f)*f
	d_perp *= -1; // -(dMot.f)*f
	d_perp += *pDmot; // d_mot - dot(d_mot, f)*f; 
	normalize( d_perp ); 
	
	double c1 = dot_product( *pDmot , d_perp ); 
	
	double sensitivity = 1.0; 
	double theta = a*sensitivity; 
	
	*pDmot = d_perp; 
	*pDmot *= c2; 
	c1 *= (1-theta); // c1*(1-theta) 
	axpy( pDmot , c1 , (d_perp) ); 
	
	
	normalize( *pDmot ); 
	phenotype.motility.migration_bias = 1.0; 

	phenotype.motility.migration_speed = pCell->custom_data[nMaxSpeed]*ECM*(1-ECM)*4.0; 
	return; 
}
