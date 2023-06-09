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

void ECM_based_cell_motility_update_including_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt )
{
    /*********************************************Chemotaxsis update***************************************************/
	
	// sample uE

	// std::cout<<"In development - do not use. needs generalized. Exiting until it is fixed!\n";  
    // std::exit(-1);  
	// std::cout<<"Cell name "<< pCell->type_name<<std::endl;

	Cell_Definition* pCD = find_cell_definition(pCell->type_name);	
	// get_single_base_behavior( Cell_Definition* pCD , std::string behavior ); 
	// get_single_base_behavior( Cell* pC, std::string behavior) ; 
	// get_single_base_behavior( pCD , "migration speed" ); 
	// get_single_base_behavior( pCell, "migration speed" ); 

	// sample ECM 
	std::vector<double> cell_position = pCell->position;
	int nearest_ecm_voxel_index = ecm.ecm_mesh.nearest_voxel_index( cell_position );   
	double ECM_density = ecm.ecm_voxels[nearest_ecm_voxel_index].density; 
	// std::cout<<"ECM density "<<ECM_density<<std::endl;
	double a = ecm.ecm_voxels[nearest_ecm_voxel_index].anisotropy; 
	// std::cout<<"ECM anisotropy "<<a<<std::endl;
	std::vector<double> f = ecm.ecm_voxels[nearest_ecm_voxel_index].ecm_fiber_alignment;
	// std::cout<<"ECM fiber alignment "<<f[0]<<", "<<f[1]<<", "<<f[2]<<std::endl;

	// Get custom data indices
	// static int max_cell_speed_index = pCell->custom_data.find_variable_index( "max_speed" ); --> comes from cell XML definition now!!!
	// if (max_cell_speed_index < 0) 
    // {
    //     std::cout << "        static int max_cell_speed_index = " <<max_cell_speed_index << std::endl;
    //     std::exit(-1);  //rwh: should really do these for each
    // }	
    
	static int chemotaxis_bias_index = pCell->custom_data.find_variable_index( "chemotaxis_bias");
	if (chemotaxis_bias_index < 0) 
    {
        std::cout << "        static int chemotaxis_bias_index = " <<chemotaxis_bias_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }

	static int ECM_sensitivity_index = pCell->custom_data.find_variable_index( "ECM_sensitivity");
	if (ECM_sensitivity_index < 0) 
    {
        std::cout << "        static int ECM_sensitivity_index = " <<ECM_sensitivity_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
    
	static int link_anisotropy_and_bias_index = pCell->custom_data.find_variable_index( "link_anisotropy_and_bias" );
	if (link_anisotropy_and_bias_index < 0) 
    {
        std::cout << "        static int link_anisotropy_and_bias_index = " <<link_anisotropy_and_bias_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
	
	static int min_ECM_mot_den_index = pCell->custom_data.find_variable_index( "min_ECM_motility_density");
    if (min_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int min_ECM_mot_den_index = " <<min_ECM_mot_den_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
	static int max_ECM_mot_den_index = pCell->custom_data.find_variable_index( "max_ECM_motility_density");
    if (max_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int max_ECM_mot_den_index = " <<max_ECM_mot_den_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
	
	static int ideal_ECM_mot_den_index = pCell->custom_data.find_variable_index( "ideal_ECM_motility_density");
    if (ideal_ECM_mot_den_index < 0) 
    {
        std::cout << "        static int ideal_ECM_mot_den_index = " <<ideal_ECM_mot_den_index << std::endl;
        std::exit(-1);  //rwh: should really do these for each
    }
	
	/****************************************Begin migration direction update****************************************/


	// ****** Start with biased random migration. We do linear combination of random and chemotaxis vectors.********//

	// Select random direction (for random portion of motility vector) and begin building updated motility direction vector
	// (note - there is NO memeory of previous direction in this model - previous ECM-based motility used the current 
	// velocity vector to build off, not a random one - this could produce divergent behaviors between models)

	// See lab note book for more notes - MUST start with random vector. In the old method I defintiely used the previous motility vector in the method, but makes no sense here!

	// get random vector - cell's "intended" or chosen random direction
	double angle = UniformRandom() * 6.283185307179586;
	std::vector<double> d_random = { cos(angle) , sin(angle) , 0.0 };

	// std::cout<<"D random "<<d_random<<std::endl;

	// get chemotaxis vector - cell's "intended" or chosen chemotaxis direction - use inbuilt chemotaxis function from PhysiCell standard models
	chemotaxis_function(pCell, phenotype, dt); 
	// std::cout<<phenotype.motility.migration_bias_direction<<std::endl; // normalized in chemotaxis function call

	//combine cell chosen random direction and chemotaxis direction (like standard update_motlity function)

	// New bias - bias such that the agents can more closely follow the gradient IF the written signals are stronger. 

	std::vector<double> d_motility = {0,0,0};

	// std::vector<double> d_motility;

	// d_motility = {0,0,0};

	if (pCell->custom_data[link_anisotropy_and_bias_index] < 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 0
	{
		d_motility = (1-a) * d_random + a * phenotype.motility.migration_bias_direction;
		// std::cout<<"d_motility "<<d_motility<<std::endl;
	}
	else if (pCell->custom_data[link_anisotropy_and_bias_index] > 0.5) // no ints/bools in custom data - must be double - so use 0.5 instead of 1
	{
		// NON-ECM linked way to signal 
		d_motility = (1-pCell->custom_data[chemotaxis_bias_index])*d_random + pCell->custom_data[chemotaxis_bias_index]*phenotype.motility.migration_bias_direction;
		// std::cout<<"d_motility "<<d_motility<<std::endl;
		// std::cout<<"Migration bias "<<phenotype.motility.migration_bias<<std::endl;
		// std::cout<<"Migration bias direction "<<phenotype.motility.migration_bias_direction<<std::endl;
		// std::cout<<"Speed "<<phenotype.motility.migration_speed<<std::endl;
		// std::cout<<"xml speed"<< pCD->phenotype.motility.migration_speed<<std::endl;
		// std::cout<<"base speed= "<<get_single_base_behavior( pCell, "migration speed" )<<std::endl; 
	}
	else
	{
		std::cout<<"Must specify reader chemotaxis modeling mode - see XML parameter \"link_anisotropy_and_bias\" Halting!!!!!!"<<std::endl;
		abort();
		return;
	}

	normalize( &d_motility ); 

	// std::cout<<"D motility "<<d_motility<<std::endl;

	// ****** Finish with linear combination of random biased migration and ECM orientation following. ********//

	// to determine direction along f, find part of d_choice that is perpendicular to f; 
	std::vector<double> d_perp = d_motility - dot_product(d_motility,f)*f; 
	
	normalize( &d_perp ); 

	// std::cout<<"D perp"<<d_perp<<std::endl;

	// std::cout<<"Fiber "<<f<<std::endl;
	
	// find constants to span d_choice with d_perp and f
	double c_1 = dot_product( d_motility , d_perp ); 
	double c_2 = dot_product( d_motility, f ); 
	// std::cout<<"f = "<<f<<std::endl;
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

		pCell->phenotype.motility.migration_speed = get_single_base_behavior( pCD, "migration speed" ) * ( 1/(rho_ideal - rho_low) * (ECM_density - rho_low)); // magnitude of direction (from ~50 lines ago) * base speed * ECM density influence
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

		pCell->phenotype.motility.migration_speed = get_single_base_behavior( pCD, "migration speed" ) * ( 1/(rho_ideal - rho_high) * (ECM_density - rho_high)); // magnitude of direction (from ~60 lines ago) * base speed * ECM density influence
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

	/*********************************************END speed update***************************************************/
	
	return; 
}


