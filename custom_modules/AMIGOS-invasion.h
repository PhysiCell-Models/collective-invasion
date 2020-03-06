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

void rightward_deterministic_cell_march (Cell* pCell , Phenotype& phenotype , double dt );

void reset_cell_position(void);

class Options 
{
 public: 
	int model = 0; 	
};



// custom cell phenotype function to scale immunostimulatory factor with hypoxia 

// set the tumor cell properties, then call the function 
// to set up the tumor cells 
void create_cell_types( void ); // done 
void ECM_setup(double numvox);
void setup_tissue(void); // done 

// set up the microenvironment to include the immunostimulatory factor 
void setup_microenvironment( void );  // done 

std::vector<std::string> AMIGOS_invasion_coloring_function( Cell* );
std::vector<std::string> ECM_anisotropy_coloring_function( Cell* );
void ecm_update_from_cell(Cell* pCell , Phenotype& phenotype , double dt); // Not currently supporting anisotropy decreasing!! 06.17.19
void ECM_informed_motility_update_w_chemotaxis ( Cell* pCell, Phenotype& phenotype, double dt );
void ECM_informed_motility_update_model_w_memory ( Cell* pCell, Phenotype& phenotype, double dt ); // Uses previous migration bias direction
void ECM_informed_motility_update_w_chemotaxis_w_variable_speed( Cell* pCell, Phenotype& phenotype, double dt ); // Sets speed based on fiber alignment. If chemo senstivitty set to 1, then cells will chemotax with speed varying strictly with alignemnt.
// void ECM_informed_motility_update_model_3 ( Cell* pCell, Phenotype& phenotype, double dt ); // uses previous velocity vector
// void ECM_informed_motility_update_model_4 ( Cell* pCell, Phenotype& phenotype, double dt ); // uses previous migration bias direction AND previous anisotropy
void change_migration_bias_vector_ecm(Cell* pCell , Phenotype& phenotype , double dt);
void run_biotransport( double t_max );
void alter_cell_uptake_secretion_saturation ( void );
void set_cell_motility_vectors(void); // Runs update_migration_bias for each cell present in a simulation
void write_ECM_Data_matlab( std::string filename );
double sign_function (double number);

// New hookean spring mechanics functions

void attach_cells( Cell* pCell_1, Cell* pCell_2 ); // used in immune_cell_attempt_attachment - copying directly. Tests to see if they are already attached using "state" What is state?
void dettach_cells( Cell* pCell_1 , Cell* pCell_2 ); // used in immune_cell_rule - copying directly

void add_elastic_velocity( Cell* pActingOn, Cell* pAttachedTo , double elastic_constant ); // calculation for velocity - copying directly
void extra_elastic_attachment_mechanics( Cell* pCell, Phenotype& phenotype, double dt ); // calls the calculation - copying directly

Cell* leader_follower_cell_check_neighbors_for_attachment( Cell* pCell , double dt ); // Used to get nearby cell pointers - cells the cell of interest could attach too - copied directly form Immun example
bool check_for_detachment( Cell* pCell , Cell* pAttachedTo , double dt ); // used to see if cells should detach - if true, dettach_cells is called. It is simple fthan than leader_follower_cell_check_neighbors_for_attachment b/c we already know a cell is attached if it is getting called. 
bool leader_follower_cell_attempt_attachment( Cell* pCell, Cell* pTarget , double dt ); // used in  immune_cell_check_neighbors_for_attachment - uses a probability based on oncoprotein amount

void leader_cell_rule( Cell* pCell, Phenotype& phenotype, double dt ); // wrapts them all together as a custom cell rule
void follower_cell_rule( Cell* pCell, Phenotype& phenotype, double dt ); // wrap them all together as a custom cell rule

