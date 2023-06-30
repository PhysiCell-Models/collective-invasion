#include "./extracellular_matrix.h"

extern ECM ecm;
extern BioFVM::Microenvironment microenvironment;

ECM_Cartesian_Mesh::ECM_Cartesian_Mesh()
{

Cartesian_Mesh();

}

ECM_Voxel::ECM_Voxel()
{

	anisotropy = 0;
	density = 0.5;
	ecm_fiber_alignment.assign( 3 , 0.0 ); 

}

ECM::ECM()
{
	ECM_Voxel template_ecm_voxel;

	// ECM_Cartesian_Mesh ecm_mesh;

	ecm_mesh.resize(1,1,1); 

	make_ecm_units();

	// ecm_voxels.push_back(template_ecm_voxel);
}

void ECM::make_ecm_units(void)
{
	// you want to ... grab the centers and ordering from ecm_mesh.voxels ... one by one. This will make it First - resize the mesh, then Second, resize the ECM. And you can make that happen in the constructor also, which might be the best idea???

	for (int i=0; i<ecm_mesh.voxels.size(); i++)
	{
		
		ecm_voxels.push_back(ecm_voxel);
		ecm_voxels[i].mesh_index = ecm_mesh.voxels[i].mesh_index;
	}


	initialize_ECM();

}

void ECM::resize_ecm_units_from_ecm_mesh(void)
{
	ecm_voxels.resize(0);

	for (int i=0; i<ecm_mesh.voxels.size(); i++)
	{	
		ecm_voxels.push_back(ecm_voxel);
		ecm_voxels[i].mesh_index = ecm_mesh.voxels[i].mesh_index;
		ecm_voxels[i].center = ecm_mesh.voxels[i].center;
		ecm_voxels[i].volume = ecm_mesh.voxels[i].volume;
	}

	initialize_ECM();
}

void ECM::initialize_ECM( void )
{
	
	if (ecm_mesh.voxels.size() != ecm_voxels.size())
	{
		std::cout<<"Resize ECM mesh to match ECM voxels before initializing ECM units to initial values"<<std::endl;
		std::cout<<" hit Enter to continue:"<<std::flush;
	 	std::cin.get();
	}

	// for( int n = 0; n < ecm_voxels.size() ; n++ )
	// {
	// 		double theta = 6.2831853071795864769252867665590 * uniform_random(); 
	// 		// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
	// 		// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
	// 		// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
	// 		ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0};

		// // For random 2-D initalization 
		// if(parameters.strings( "ECM_orientation_setup") == "random")
		// {
		// 	double theta = 6.2831853071795864769252867665590 * uniform_random(); 
		// 	// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
		// 	// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
		// 	// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
		// 	ecm_voxels[n].ecm_fiber_alignment = {cos(theta), sin(theta), 0.0};
		// }

		// // for starburst initialization 
		// else if(parameters.strings( "ECM_orientation_setup") == "starburst")
		// {
		// 	std::vector<double> position = microenvironment.mesh.voxels[n].center; 
		// 	normalize( &position ); 
		// 	ecm_voxels[n].ecm_fiber_alignment =  { position[0],position[1],0}; // oriented out (perpindeicular to concentric circles)
		// 	normalize(&ecm_voxels[n].ecm_fiber_alignment );
		// }

		// // for circular initialization 
		// else if(parameters.strings( "ECM_orientation_setup") == "circular")
		// {
		// 	std::vector<double> position = microenvironment.mesh.voxels[n].center; 
		// 	normalize( &position );
		// 	ecm_voxels[n].ecm_fiber_alignment  =  { position[1],-position[0],0}; // oriented in cirlce
		// 	normalize(&ecm_voxels[n].ecm_fiber_alignment );
		// }

		// else if(parameters.strings( "ECM_orientation_setup") == "horizontal")
		// {
		// 	ecm_voxels[n].ecm_fiber_alignment  =  { 1.0, 0.0, 0.0}; 
		// 	normalize(&ecm_voxels[n].ecm_fiber_alignment );
		// }

		// else if(parameters.strings( "ECM_orientation_setup") == "vertical")
		// {
		// 	ecm_voxels[n].ecm_fiber_alignment =  { 0.0, 1.0, 0.0}; 
		// 	normalize(&ecm_voxels[n].ecm_fiber_alignment );
		// }

		// else
		// {
		// 	std::cout<<"WARNING: NO ECM ORIENTATION SPECIFIED. FIX THIS!!!"<<std::endl;
		// 	std::cout<<"Halting program!!!"<<std::endl;
		// 	abort();
		// 	return;
		// }


		// if(parameters.ints("unit_test_setup")==1 && parameters.ints("march_unit_test_setup") == 0)
		// {

		// 	ecm_voxels[n].density = 0.5;
		// 	ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");

		// }

		// else if (parameters.ints("unit_test_setup") == 1 && parameters.ints("march_unit_test_setup") == 1)
		// {
			
		// 	ecm_voxels[n].density = 0.5;
		// 	ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");
			
		// }

		// else if(parameters.ints("unit_test_setup") == 0 && parameters.ints("march_unit_test_setup") == 0)
		// {
		// 	ecm_voxels[n].density = parameters.doubles("initial_ECM_density");
		// 	ecm_voxels[n].anisotropy = parameters.doubles("initial_anisotropy");
			
		// }

		// else
		// {
		// 	std::cout<<"ECM density and anisotropy not set correctly!!!! FIX!!!!!!!!!"<<std::endl;
		// 	std::cout<<"Halting!"<<std::endl;
		// 	abort();
		// 	return;
		// }
			// if(parameters.strings( "ECM_orientation_setup") == "random")
			// {
			// 	double theta = 6.2831853071795864769252867665590 * uniform_random(); 
			// 	// ecm.ecm_data[i].ECM_orientation[0] = cos(theta);
			// 	// ecm.ecm_data[i].ECM_orientation[1] = sin(theta);
			// 	// ecm.ecm_data[i].ECM_orientation[2] = 0.0;
			// 	ECM_fiber_alignment[n] = {cos(theta), sin(theta), 0.0};
			// }
		// }
		
	}

	

ECM_options::ECM_options()
{
 
}

// ECM ecm;