#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

void copy_ECM_data_to_BioFVM( Cell* pCell, Phenotype& phenotype, double dt );

double dot_product_ext( const std::vector<double>& v , const std::vector<double>& w );

double sign_function (double number);

// uses cell motility vector for realigning ECM. 
void ECM_remodeling_function( Cell* pCell, Phenotype& phenotype, double dt );

void ECM_based_speed_update( Cell* pCell, Phenotype& phenotype, double dt );

void custom_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt);

// uses custom chemotaxis function - will upgrade to use the standard one later.
void ECM_and_chemotaxis_based_cell_migration_update( Cell* pCell, Phenotype& phenotype, double dt );

void combined_ECM_remodeling_and_speed_update( Cell* pCell, Phenotype& phenotype, double dt);