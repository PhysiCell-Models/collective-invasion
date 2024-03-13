#include "./extracellular_matrix.h"

using namespace BioFVM; 
using namespace PhysiCell;

extern ECM ecm;




void copy_ECM_data_to_BioFVM( void )
{


	// This enables the use of rules to change the cell behaviors as well as aspects of visualization in the Studio
	// this is ONE WAY!!! DO NOT MODIFY WITHING BIOFVM!!! OR USE CELLS TO CHANGE THIS!!!!!!!!!

	// found it greatly increase run time - could possibly improve this by calling at the phenotype time step - thats when rules are applied
	// look into this later - for now just going the route of signals/behaviors or even raw Hill functions
	// Note also that if you want to use the ECM to update a mechanics thing (like speed), rules might not be the right approach anyway (called to infrequently) - it would depend on the concept/intend - is it a phenotypic change or a physical reaction to environment?

    int number_of_voxels = ecm.ecm_mesh.voxels.size();
	
	static int ECM_anisotropy_index = BioFVM::microenvironment.find_density_index( "ECM_anisotropy" ); 
	if (ECM_anisotropy_index < 0) 
    {
        std::cout << "        static int ECM_anisotropy_index = " <<ECM_anisotropy_index << std::endl;
		std::cout << "        ADD ECM_anisotropy field to your simulation!!!" << std::endl;
		std::cout << "        This feature only works with ECM_anisotropy field added to your simulation. Set decay and diffusion constant to 0" << std::endl;
		std::cout << "        AND match ECM element size to the diffusion voxel size!!!" << std::endl;
		std::cout << "        Halting!!!!!!" << std::endl;
        std::exit(-1);   
    }
	static int ECM_density_index = BioFVM::microenvironment.find_density_index( "ECM_density" ); 
	if (ECM_density_index < 0) 
    {
        std::cout << "        static int ECM_density_index = " <<ECM_density_index << std::endl;
		std::cout << "        ADD ECM_density field to your simulation!!!" << std::endl;
		std::cout << "        This feature only works with ECM_density field added to your simulation. Set decay and diffusion constant to 0!!!" << std::endl;
		std::cout << "        AND match ECM element size to the diffusion voxel size!!!" << std::endl;
		std::cout << "        Halting!!!!!!" << std::endl;
        std::exit(-1);   
    }
    for( int n=0; n < number_of_voxels ; n++ )
    {
		BioFVM::microenvironment.density_vector(n)[ECM_anisotropy_index] = ecm.ecm_voxels[n].anisotropy;
		BioFVM::microenvironment.density_vector(n)[ECM_density_index] = ecm.ecm_voxels[n].density;
    }

    return;
}

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

void custom_update_motility_vector(Cell* pCell, Phenotype& phenotype, double dt_ )
{
	// Modified version of standard update_motility_vector function in PhysiCell. Search "update_motility_vector" to find original
	// Takes into account ECM anisotropy for anisotropy linked chemotactic migration
	// Also, outputs/modifies migration_bias_direction versus motility_vector (as in the original function)

	if( phenotype.motility.is_motile == false )
	{
		phenotype.motility.motility_vector.assign( 3, 0.0 ); 
		return; 
	}
	
	if( UniformRandom() < dt_ / phenotype.motility.persistence_time || phenotype.motility.persistence_time < dt_ )
	{
		static int link_anisotropy_and_bias_index = pCell->custom_data.find_variable_index( "link_anisotropy_and_bias" );
		if (link_anisotropy_and_bias_index < 0) 
		{
			std::cout << "        static int link_anisotropy_and_bias_index = " <<link_anisotropy_and_bias_index << std::endl;
			std::exit(-1);   
		}

		static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
		if (ECM_sensitivity_index < 0) 
		{
			std::cout << "        static int ECM_sensitivity_index = " <<ECM_sensitivity_index << std::endl;
			std::exit(-1);   
		}

		std::vector<double> randvec(3,0.0);
		if( phenotype.motility.restrict_to_2D == true )
		{ randvec = UniformOnUnitCircle(); }
		else
		{ randvec = UniformOnUnitSphere(); }
		
		// if the update_bias_vector function is set, use it  
		if( pCell->functions.update_migration_bias )
		{
			pCell->functions.update_migration_bias( pCell, phenotype, dt_ ); 
		}

		if (pCell->custom_data[link_anisotropy_and_bias_index] < 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 0
		{
			// sample ECM 
			std::vector<double> cell_position = pCell->position;
			double s = pCell->custom_data[ECM_sensitivity_index]; 
			int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
			double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 
			double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; 

			phenotype.motility.migration_bias_direction *= s; // motility = bias*bias_vector 
			phenotype.motility.migration_bias_direction *= a; // motility = bias*bias_vector 
			double one_minus_anistropy_s = 1.0 - a*s;
			axpy( &(phenotype.motility.migration_bias_direction), one_minus_anistropy_s, randvec ); // motility = bias*bias_vector + (1-bias)*randvec 
			// y = y + a*x 
			// void axpy( std::vector<double>* y, double& a , std::vector<double>& x );
			// std::cout<<"ECM_density = "<<ECM_density<<std::endl;
		}

		else if (pCell->custom_data[link_anisotropy_and_bias_index] > 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 1
		{
			phenotype.motility.migration_bias_direction *= phenotype.motility.migration_bias; // motility = bias*bias_vector 
		
			double one_minus_bias = 1.0 - phenotype.motility.migration_bias; 
			axpy( &(phenotype.motility.migration_bias_direction), one_minus_bias, randvec ); // motility = (1-bias)*randvec + bias*bias_vector

		}
		else
		{
			std::cout<<"Must specify reader chemotaxis modeling mode - see XML parameter \"link_anisotropy_and_bias\" Halting!!!!!!"<<std::endl;
			abort();
			return;
		}

		normalize( &(phenotype.motility.migration_bias_direction) ); 
	}
	return;
}

void ECM_to_cell_interactions( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Function calculates takes chemotaxis direction and combines with fiber following to produce combined random biased motility and 
	// fiber following vector. Then calculates agent speed based on ECM density. 

	Cell_Definition* pCD = find_cell_definition(pCell->type_name);	

	// sample ECM 
	std::vector<double> cell_position = pCell->position;

	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 
	double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; 
	std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;

	// Get cell level values
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis_bias");
	if (chemotaxis_bias_index < 0) 
    {
        std::cout << "        static int chemotaxis_bias_index = " <<chemotaxis_bias_index << std::endl;
        std::exit(-1);   
    }

	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
	if (ECM_sensitivity_index < 0) 
    {
        std::cout << "        static int ECM_sensitivity_index = " <<ECM_sensitivity_index << std::endl;
        std::exit(-1);   
    }
    
	static int link_anisotropy_and_bias_index = pCell->custom_data.find_variable_index( "link_anisotropy_and_bias" );
	if (link_anisotropy_and_bias_index < 0) 
    {
        std::cout << "        static int link_anisotropy_and_bias_index = " <<link_anisotropy_and_bias_index << std::endl;
        std::exit(-1);   
    }
	
	static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min_ECM_motility_density");
    if (min_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int min_ECM_mot_den_index = " <<min_ECM_mot_den_index << std::endl;
        std::exit(-1);   
    }
	static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max_ECM_motility_density");
    if (max_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int max_ECM_mot_den_index = " <<max_ECM_mot_den_index << std::endl;
        std::exit(-1);   
    }
	
	static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal_ECM_motility_density");
    if (ideal_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int ideal_ECM_mot_den_index = " <<ideal_ECM_mot_den_index << std::endl;
        std::exit(-1);   
    }

	static int migration_bias_norm_index = pCell->custom_data.find_variable_index( "migration_bias_norm");
    if (migration_bias_norm_index < 0) 
    {
        std::cout << "        static int migration_bias_norm_index = " <<migration_bias_norm_index << std::endl;
        std::exit(-1);   
    }

	static int rules_based_speed_multiplier_index = pCell->custom_data.find_variable_index( "rules_based_speed_multiplier");
	if (rules_based_speed_multiplier_index < 0) 
	{
		std::cout << "        static int rules_based_speed_multiplier_index = " <<rules_based_speed_multiplier_index << std::endl;
		std::exit(-1);  
	}

	// ******************************************************************************************************//
	// ****** Make linear combination of random biased migration (motility_vector) and ECM orientation following. ********//
	// ******************************************************************************************************//
	
	if (pCell->custom_data[ECM_sensitivity_index] < 0.0001)
	{
		// if ECM_sensitivity is very low, then the migration bias vector is used as the motility vector
		phenotype.motility.motility_vector = phenotype.motility.migration_bias_direction;
	}

	else
	{
		std::vector<double> d_motility = {0,0,0};

		d_motility = phenotype.motility.migration_bias_direction; 

		normalize( &d_motility ); 

		// to determine direction along f, find part of d_choice that is perpendicular to f; 
		
		// There are rare times where there istrouble with the accuracy of the approximation in one or both of the following lines - needs polished in the future. 
		std::vector<double> d_perp = d_motility - dot_product(d_motility,f)*f; 	
		normalize( &d_perp );


		double c_1 = dot_product( d_motility , d_perp ); 
		double c_2 = dot_product( d_motility, f ); 

		// calculate bias away from directed motitility - combination of sensitity to ECM and anisotropy
		double gamma = pCell->custom_data[ECM_sensitivity_index] * a; // at low values, directed motility vector is recoved. At high values, fiber direction vector is recovered.

		phenotype.motility.motility_vector = (1.0-gamma)*c_1*d_perp + c_2*f;
		
		if(parameters.bools("normalize_ECM_influenced_motility_vector") == true)
		{
			// if the vector is to be normalized, we, by definition, already know the magnitude will be 1.0
			pCell->custom_data[migration_bias_norm_index] = 1.0;
		}
		else  
		{
			pCell->custom_data[migration_bias_norm_index] = norm( phenotype.motility.motility_vector);
		}

		normalize( &(phenotype.motility.motility_vector) ); 
	}
	/****************************************END new migration direction update****************************************/

	/*********************************************Begin speed update***************************************************/
	
	// needed to reassign speed after update.
	
	double rho_low = pCell->custom_data[min_ECM_mot_den_index];
	double rho_high = pCell->custom_data[max_ECM_mot_den_index];
	double rho_ideal = pCell->custom_data[ideal_ECM_mot_den_index];

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
		// So finally: speed = normalized_value * max_speed * rules_modifier * (1/(rho_ideal - rho_l) * (rho - rho_l))

		// pCell->phenotype.motility.migration_speed = pCell->custom_data[migration_bias_norm_index] * get_single_base_behavior( pCD, "migration speed" ) * rules_modifier * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); 

		pCell->phenotype.motility.migration_speed = pCell->custom_data[migration_bias_norm_index];
		pCell->phenotype.motility.migration_speed *= get_single_base_behavior( pCD, "migration speed" ); 
		pCell->phenotype.motility.migration_speed *= pCell->custom_data[rules_based_speed_multiplier_index];
		pCell->phenotype.motility.migration_speed *= ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); 
	}

	else if (rho_ideal < ECM_density && ECM_density < rho_high )
	{

		// for base speed: y - y_1 = m (x - x_1) or y = m (x - x_1) + y_1
		// Assuming that y_1 = 0 --> y = m (x - x_1)
		// m = rise/run = (speed(rho_ideal) - speed(rho_l)/(rho_ideal - rho_l)). Same for rho_h
		// Assuming that speed(rho_ideal) = 1.0 and speed(rho_l (or rho_h)) = 0.0, m = 1/(rho_ideal - rho_l)
		// y = 1/(x_2 - x_1) * (x - x_1) --> speed_base = 1/(rho_ideal - rho_l) * (rho - rho_l)
		// So finally: speed = normalized_value * max_speed * rules_modifier * (1/(rho_ideal - rho_l) * (rho - rho_l))

		// pCell->phenotype.motility.migration_speed = pCell->custom_data[migration_bias_norm_index] * get_single_base_behavior( pCD, "migration speed" ) * rules_modifier * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); 
		
		pCell->phenotype.motility.migration_speed = pCell->custom_data[migration_bias_norm_index]; 
		pCell->phenotype.motility.migration_speed *= get_single_base_behavior( pCD, "migration speed" ); 
		pCell->phenotype.motility.migration_speed *= pCell->custom_data[rules_based_speed_multiplier_index];
		pCell->phenotype.motility.migration_speed *= ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); 
	}

	else //if (ECM_density >= rho_high)
	{
		pCell->phenotype.motility.migration_speed = 0.0;
	}

	phenotype.motility.motility_vector *= phenotype.motility.migration_speed;

	/*********************************************END speed update***************************************************/
	return; 
}

void ECM_remodeling_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// uses cell motility vector for realigning ECM (versus velocity vector)

	// Find attributes needed for updating ECM
	std::vector<double> cell_position = pCell->position;
	
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   

	static int Cell_ECM_target_density_index = pCell->custom_data.find_variable_index( "target_ECM_density");
	if (Cell_ECM_target_density_index < 0) 
    {
        std::cout << "        static int Cell_ECM_target_density_index = " <<Cell_ECM_target_density_index << std::endl;
        std::exit(-1);   
    }
	static int Cell_ECM_production_rate_index = pCell->custom_data.find_variable_index( "ECM_production_rate");
	if (Cell_ECM_production_rate_index < 0) 
    {
        std::cout << "        static int Cell_ECM_production_rate_index = " <<Cell_ECM_production_rate_index << std::endl;
        std::exit(-1);   
    }
	static int Cell_anistoropy_rate_of_increase_index = pCell->custom_data.find_variable_index( "Anisotropy_increase_rate");
	if (Cell_anistoropy_rate_of_increase_index < 0) 
    {
        std::cout << "        static int Cell_anistoropy_rate_of_increase_index = " <<Cell_anistoropy_rate_of_increase_index << std::endl;
        std::exit(-1);   
    }	
	static int Cell_fiber_realignment_rate_index = pCell->custom_data.find_variable_index( "Fiber_realignment_rate");
	if (Cell_fiber_realignment_rate_index < 0) 
    {
        std::cout << "        static int Cell_fiber_realignment_rate_index = " <<Cell_fiber_realignment_rate_index << std::endl;
        std::exit(-1);   
    }	

    // Cell-ECM density interaction

    double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density;
    double r = pCell->custom_data[Cell_ECM_production_rate_index];

    ecm.ecm_voxels[nearest_ecm_voxel_index].density = ECM_density + r * dt  * (pCell->custom_data[Cell_ECM_target_density_index] - ECM_density);


	// End Cell-ECM density interaction

	// Cell-ECM Fiber realignment and anisotropy remodeling - continous then instantaneous

	if( parameters.ints("discrete_ECM_remodeling") == 1)
	{	
		std::vector<double> ECM_orientation = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment; 

		double anisotropy = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy;
		double migration_speed = pCell->phenotype.motility.migration_speed;
		double r_0 = pCell->custom_data[Cell_fiber_realignment_rate_index]*migration_speed; 

		double r_realignment = r_0 * (1-anisotropy);
		double ddotf;
		std::vector<double> norm_cell_motility = phenotype.motility.motility_vector;

		normalize(&norm_cell_motility);

		ddotf = dot_product_ext(ECM_orientation, norm_cell_motility);
		ECM_orientation = sign_function(ddotf) * ECM_orientation; // flips the orientation vector so that it is aligned correctly with the moving cell for proper reoirentation later.
		std::vector<double> f_minus_d;
		f_minus_d.resize(3,0.0);

		for(int i = 0; i < 3; i++)
		{
			if (ddotf<0.0)
			{
				ECM_orientation = -1.0 * ECM_orientation;
			}
			f_minus_d[i] = ECM_orientation[i] - norm_cell_motility[i]; 
			ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment[i] -= dt * r_realignment * f_minus_d[i]; 
		}

		normalize(&(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment)); 

		// End Cell-ECM Fiber realingment

		// Cell-ECM Anisotropy Modification

		double r_a0 = pCell->custom_data[Cell_anistoropy_rate_of_increase_index] ; 

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
			// Cell-ECM Fiber realignment
			if (norm(phenotype.motility.motility_vector) == 0)
			{std::cout<<"Motility vector norm = 0"<<std::endl;}

			// Test for lack of realignment parameters (non-remodeling cells)
			if (pCell->custom_data[Cell_anistoropy_rate_of_increase_index] < 0.00001 ||pCell->custom_data[Cell_fiber_realignment_rate_index] < 0.00001 )
			{
				// std::cout<<"pass code"<<std::endl;
			}
			
			// if a remodeling cell, then realign
			else
			{

				ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment = phenotype.motility.motility_vector;

				normalize(&(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment));

				if (norm(ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment) == 0)
				{std::cout<<"ECM orientation vector norm = 0"<<std::endl;}

				ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy = 1;
			}
		}
	}

	else
	{
		std::cout<<"Must specify ECM remodeling mode - see XML parameter \"discrete_ECM_remodeling\" Halting!!!!!!"<<std::endl;
		abort();
		return;
	}
}

void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	// Replaces the standard_update_cell_velocity in fiber following/senstive agents. 
	// Assign this to update_cell_velocity function pointer to get fiber following and density based speed changes
	

	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	pCell->state.neighbors.clear(); 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		pCell->add_potentials(*neighbor);
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			pCell->add_potentials(*neighbor);
		}
	}

	// non-standard motility update - part of ECM package
	custom_update_motility_vector(pCell, phenotype, dt);

	// ECM following and speed update
	ECM_to_cell_interactions(pCell, phenotype, dt); 

	// standard update cell velocity - after this update proceeds as "conventional" PhysiCell
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}

void create_default_ECM_compatible_agent( void )
{
	initialize_default_cell_definition(); 
	// cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	cell_defaults.functions.update_velocity = custom_update_cell_velocity;

	cell_defaults.functions.custom_cell_rule = ECM_remodeling_function; 

	// cell_defaults.functions.custom_cell_rule = NULL; 

	cell_defaults.functions.update_phenotype = NULL; 

	// 	cell_defaults.functions.cycle_model = Ki67_advanced; 
	
	// cell_defaults.functions.volume_update_function = standard_volume_update_function;
	// cell_defaults.functions.update_migration_bias = NULL; 
	
	// cell_defaults.functions.update_phenotype = update_cell_and_death_parameters_O2_based; // NULL; 
	// cell_defaults.functions.custom_cell_rule = NULL; 
	
	// cell_defaults.functions.update_velocity = standard_update_cell_velocity;
	// cell_defaults.functions.add_cell_basement_membrane_interactions = NULL; 
	// cell_defaults.functions.calculate_distance_to_membrane = NULL; 
	
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
/* 
 if( ECM.TellRows() > 0 )
 {
  // find the k corresponding to z_slice
  
  
  
  Vector position; 
  *position(2) = z_slice; 
  

  // 25*pi* 5 microns^2 * length (in source) / voxelsize^3
  
  for( int j=0; j < ratio*ECM.TellCols() ; j++ )
  {
   // *position(1) = *Y_environment(j); 
   *position(1) = *Y_environment(0) - dy_stroma/2.0 + j*voxel_size + half_voxel_size; 
   
   for( int i=0; i < ratio*ECM.TellRows() ; i++ )
   {
    // *position(0) = *X_environment(i); 
    *position(0) = *X_environment(0) - dx_stroma/2.0 + i*voxel_size + half_voxel_size; 
	
    double E = evaluate_Matrix3( ECM, X_environment , Y_environment, Z_environment , position );	
	double BV = normalizer * evaluate_Matrix3( OxygenSourceHD, X_environment , Y_environment, Z_environment , position );
	if( isnan( BV ) )
	{ BV = 0.0; }

	vector<string> Colors;
	Colors = hematoxylin_and_eosin_stroma_coloring( E , BV );
	Write_SVG_rect( os , *position(0)-half_voxel_size-X_lower , *position(1)-half_voxel_size+top_margin-Y_lower, 
	voxel_size , voxel_size , 1 , Colors[0], Colors[0] );
   
   }
  }
 
 }
*/
	os << "  </g>" << std::endl; 
 
	// Now draw vessels

	/*
	 std::vector<std::string> VesselColors = hematoxylin_and_eosin_stroma_coloring( 0,1 );

	 os << " <g id=\"BloodVessels\">" << endl; 
	 extern vector<BloodVesselSegment*> BloodVesselSegments; 
	 Vector Offset; 
	 *Offset(0) = X_lower; 
	 *Offset(1) = Y_lower-top_margin;
	*/
 

 
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


	// plot intersecting BM points
	/* 
	 for( int i=0 ; i < BasementMembraneNodes.size() ; i++ )
	 {
		// vector<string> Colors = false_cell_coloring( pC ); 
		BasementMembraneNode* pBMN = BasementMembraneNodes[i]; 
		double thickness =0.1; 
		
		if( fabs( *(pBMN->Position)(2) - z_slice ) < thickness/2.0 ) 
		{
		 string bm_color ( "rgb(0,0,0)" );
		 double r = thickness/2.0; 
		 double z = fabs( *(pBMN->Position)(2) - z_slice) ; 

		 os << " <g id=\"BMN" << pBMN->ID << "\">" << std::endl; 
		 Write_SVG_circle( os,*(pBMN->Position)(0)-X_lower, *(pBMN->Position)(1)+top_margin-Y_lower, 10*thickness/2.0 , 0.5 , bm_color , bm_color ); 
		 os << " </g>" << std::endl;
		}
		// pC = pC->pNextCell;
	 }
	*/ 
	
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