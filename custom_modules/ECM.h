#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

namespace PhysiCell{
	
class ECM_DATA
{
	private:
	public:
	   double density, anisotropy;
	   std::vector<double> ECM_orientation;
	   ECM_DATA();
};//close ECM_DATA

class ECM
{
 private:
 public: 
	std::vector<ECM_DATA> ecm_data;
	Cartesian_Mesh mesh;
	Microenvironment* pMicroenvironment; 
	ECM();
	
	void sync_to_BioFVM( void );
	
	
};//close ECM

extern ECM ecm;
extern ECM_DATA ecm_data;

//void ECM_setup(double numvox);

void cell_update_from_ecm(void);
void change_speed_ecm(Cell* pCell);
void change_migration_bias_vector_ecm(Cell* pCell);
void change_bias_ecm(Cell* pCell);
void ecm_update_from_cell(double numvox, double dt);
void change_ecm_density(int ecm_index, double dt);

void write_ECM_Data_matlab( std::string filename );

};//close namespace
