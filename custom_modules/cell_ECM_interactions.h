#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

double dot_product_ext( const std::vector<double>& v , const std::vector<double>& w );

double sign_function (double number);

// uses cell motility vector for realigning ECM. 
void ECM_remodeling_function( Cell* pCell, Phenotype& phenotype, double dt );

// uses custom chemotaxis function - will upgrade to use the standard one later.
void ECM_based_cell_motility_update_with_chemotaxis( Cell* pCell, Phenotype& phenotype, double dt );