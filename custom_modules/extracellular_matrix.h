#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

class ECM_Cartesian_Mesh : public Cartesian_Mesh
{
	public:
	ECM_Cartesian_Mesh();


};

class ECM_Voxel : public Voxel
{
	public:

	ECM_Voxel();

	double anisotropy;
	double density;
	std::vector<double> ecm_fiber_alignment;

};

class ECM
{
	// so in theory, anythign that doens't need accessed outside of the class is private, to aid in debugging and not messing thigns up. Interesitng. 
	
	public:
	
	ECM_Voxel ecm_voxel;
	ECM_Cartesian_Mesh ecm_mesh;
	std::string spatial_units; 
	std::string name; 
	
	std::vector<ECM_Voxel> ecm_voxels;



	ECM();
	ECM(std::string name);

	void make_ecm_units (void);
	void resize_ecm_units_from_ecm_mesh(void); // destroys previous ECM vector. New ECM_units are put in place and any previous initial configuration is lost
	void initialize_ECM(void);

};

class ECM_options
{

	private:
 
	public: 
	ECM* pECM;
	std::string name; 

	

	std::string spatial_units; 
	double dx;
	double dy; 
	double dz; 

	// Currently only used to specify FIBER orientation, but could be expanded in the future
	std::string ecm_initial_configuration;

	ECM_options(); // needs defined!!!
};

// extern ECM_options default_ecm_options; 
// extern ECM ecm;
