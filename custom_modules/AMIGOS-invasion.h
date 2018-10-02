/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version 1.2.2) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 2017 (in review).                       #
#     preprint DOI: 10.1101/088773                                            #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#    llelized diffusive transport solver for 3-D biological simulations,      #
#    Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730   #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2017, Paul Macklin and the PhysiCell Project             #
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

#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

extern Cell_Definition leader_cell; 
extern Cell_Definition follower_cell; 

// overall rules 

void tumor_cell_phenotype_with_oncoprotein( Cell* pCell , Phenotype& phenotype , double dt ) ;// done 
void chemotaxis_oxygen( Cell* pCell , Phenotype& phenotype , double dt ); // done 
//void change_speed_ecm(Cell* pCell);

// follower cell rules 

void follower_cell_phenotype_model0( Cell* pCell , Phenotype& phenotype , double dt ); 

// leader cell rules

void leader_cell_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt ); 

void follower_cell_phenotype_model( Cell* pCell , Phenotype& phenotype , double dt ); 

void leader_cell_motility_model0( Cell* pCell , Phenotype& phenotype , double dt ); 

void switching_phenotype_model( Cell* pCell, Phenotype& phenotype, double dt ); 

class Options 
{
 public: 
	int model = 0; 	
};



// custom cell phenotype function to scale immunostimulatory factor with hypoxia 

// set the tumor cell properties, then call the function 
// to set up the tumor cells 
void create_cell_types( void ); // done 
void ECM_setup(double numvox, double radius);
void setup_tissue(double radius); // done 

// set up the microenvironment to include the immunostimulatory factor 
void setup_microenvironment( void );  // done 

std::vector<std::string> AMIGOS_invasion_coloring_function( Cell* );
std::vector<std::string> ECM_anisotropy_coloring_function( Cell* );
void ecm_update_from_cell(Cell* pCell , Phenotype& phenotype , double dt); // NOTE - not currently supporting ECM density increasing or anisotropy decreasing!!! 03.30.18
void change_migration_bias_vector_ecm(Cell* pCell , Phenotype& phenotype , double dt);
void change_speed_ecm(Cell* pCell , Phenotype& phenotype , double dt);
void run_biotransport( double t_max );
