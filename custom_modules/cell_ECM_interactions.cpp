#include "./extracellular_matrix.h"

using namespace BioFVM; 
using namespace PhysiCell;

extern ECM ecm;

double dot_product_ext( const std::vector<double>& v , const std::vector<double>& w )
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

void ECM_based_cell_motility_update_with_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
    /*********************************************Chemotaxsis update***************************************************/
	
	// sample uE

	std::cout<<"In development - do not use. needs generalized. Exiting until it is fixed!\n";  
    std::exit(-1);  


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

	// static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max speed" ); 
	// static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis bias");
	// static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM sensitivity");
	// static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min ECM motility density");
	// static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max ECM motility density");
	// static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal ECM motility density");

	// rwh: use underscores now that they are in the .xml as tags
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


// uses cell motility vector for realigning ECM. 
void ECM_remodeling_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	// this is based in ecm_update_from_cell_motility_vector from PC_ECM_extension v.1.x
	
	//*********************************** REMODELING ***********************************//
	
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
    
	// std::cout<<"ECM density = "<<ecm.ecm_voxels[nearest_ecm_voxel_index].density<<std::endl;
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
    
	//****************************************** END REMODELING ******************************************************//

	/*********************************************Begin speed update***************************************************/
	
	// // New speed update (06.18.19) - piece wise continous

	// // rwh: use underscores now that they are in the .xml as tags
	// static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max_speed" ); 
	// if (max_cell_speed_index < 0) std::exit(-1);
    // // std::cout << "        static int max_cell_speed_index = " <<max_cell_speed_index << std::endl;
	// // static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis_bias");
	// // static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
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
	
	// double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	// double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	// double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];
	// // std::cout<<"rho_low = "<<rho_low<<std::endl;
	// // std::cout<<"rho_high = "<<rho_high<<std::endl;
	// // std::cout<<"rho_ideal = "<<rho_ideal<<std::endl;
	// // std::cout<<"ECM_density = "<<ECM_density<<std::endl;

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

	// 	// std::cout<<"l2438"<<std::endl;

	// 	pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
	// 	// std::cout<<"migration speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// 	// std::cout<<"max_cell_speed = "<<pCell->custom_data[max_cell_speed_index]<<std::endl;
	// 	// std::cout<<"speed = "<<pCell->phenotype.motility.migration_speed<<std::endl;
	// }

	// else if (rho_ideal < ECM_density && ECM_density < rho_high )
	// {

	// 	// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
	// 	// Assuming that y_1 = 0 --> y = m (x - x_1)
	// 	// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
	// 	// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
	// 	// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
	// 	// So finally: speed = max_speed * (1/(rho_ideal - rho_l) * (rho - rho_l))

	// 	pCell->phenotype.motility.migration_speed = pCell->custom_data[max_cell_speed_index] * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
	// }

	// else //if (ECM_density >= rho_high)
	// {
	// 	pCell->phenotype.motility.migration_speed = 0.0;
	// }

    // return;
}

