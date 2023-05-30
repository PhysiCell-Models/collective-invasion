#include "./extracellular_matrix.h"

using namespace BioFVM; 
using namespace PhysiCell;

void test( Cell* pCell, Phenotype& phenotype, double dt )
{
	std::cout << "Hello, world! I am all out of bubble gum" << std::endl;
}

// void ECM_based_cell_motility_update( Cell* pCell, Phenotype& phenotype, double dt )
// {
    // /*********************************************Chemotaxsis update***************************************************/
	
	// // sample uE
	// static int inflam_sig_index = microenvironment.find_density_index( "inflammatory_signal" ); 
	
	// phenotype.motility.is_motile = true; 
	// phenotype.motility.migration_bias = parameters.doubles("inflam_sig_migration_bias_for_fibroblasts"); //0.95;
	// phenotype.motility.migration_bias_direction = pCell->nearest_gradient(inflam_sig_index);

	// normalize( &( phenotype.motility.migration_bias_direction ) );

	// /*********************************************Begin speed update***************************************************/

	// // sample ECM - only changes for decoupling **should** be here as nothign gets written to the ECM...
	// std::vector<double> cell_position = pCell->position;
	// int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	// double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 

	// // static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max speed" ); 
	// // static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis bias");
	// // static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM sensitivity");
	// // static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min ECM motility density");
	// // static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max ECM motility density");
	// // static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal ECM motility density");

	// // rwh: use underscores now that they are in the .xml as tags
	// static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max_speed" ); 
    // // std::cout << "        static int max_cell_speed_index = " <<max_cell_speed_index << std::endl;
	// static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis_bias");
	// static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
    // // std::cout << "        static int ECM_sensitivity_index = " <<ECM_sensitivity_index << std::endl;
	// static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min_ECM_motility_density");
    // if (min_ECM_mot_den_index < 0) 
    // {
    //     std::cout << "        static int min_ECM_mot_den_index = " <<min_ECM_mot_den_index << std::endl;
    //     std::exit(-1);  //rwh: should really do these for each
    // }
	// static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max_ECM_motility_density");
    // if (max_ECM_mot_den_index < 0) std::exit(-1);
	// static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal_ECM_motility_density");
    // if (ideal_ECM_mot_den_index  < 0) std::exit(-1);
	
	// // New speed update (06.18.19) - piece wise continous
	
	// double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	// double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	// double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];

	// // std::cout<<"ECM_density = "<<ECM_density<<std::endl;
	// pCell->phenotype.motility.migration_speed = 1.0;
	// if (ECM_density <= rho_low)
	// {
	// 	pCell->phenotype.motility.migration_speed = 0.0;

	// }

	// else if (rho_low < ECM_density && ECM_density <= rho_ideal)
	// {

	// 	// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
	// 	// Assuming that y_1 = 0 --> y = m (x - x_1)
	// 	// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
	// 	// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
	// 	// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
	// 	// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))
	// 	pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
	// 	// std::cout<<"max_cell_speed = "<<pCell->custom_data[max_cell_speed_index]<<std::endl;
	// 	// std::cout<<"rho_ideal = "<<rho_ideal<<std::endl;
	// 	// std::cout<<"rho_low = "<<rho_low<<std::endl;
	// 	// std::cout<<"rho_high = "<<rho_high<<std::endl;
	// 	// std::cout<<"speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// 	// pCell->phenotype.motility.migration_speed = 10.0;
	// }

	// else if (rho_ideal < ECM_density && ECM_density < rho_high )
	// {

	// 	// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
	// 	// Assuming that y_1 = 0 --> y = m (x - x_1)
	// 	// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
	// 	// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
	// 	// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
	// 	// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

	// 	pCell->phenotype.motility.migration_speed *= pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	// }

	// else //if (ECM_density >= rho_high)
	// {
	// 	pCell->phenotype.motility.migration_speed = 0.0;
	// }

	// if(phenotype.death.dead == true)
	// {
	// 	phenotype.motility.is_motile=false;
	// 	pCell->functions.update_phenotype = NULL;
	// 	pCell->functions.update_migration_bias = NULL;
	// 	std::cout<<"Cell is dead"<<std::endl;


	// }	
	// // std::cout<<"Volume= "<<phenotype.volume.total<<std::endl;
   	// // std::cout<<pCell->phenotype.motility.migration_speed<<std::endl;
	
	// return; 
// }

// void ECM_remodeling_function( Cell* pCell, Phenotype& phenotype, double dt )
// {
//     std::cout << "Hello, world! I am all out of bubble gum" << std::endl;
// }

