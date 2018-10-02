#include "./ECM.h"

namespace PhysiCell{
	
ECM ecm;
//ECM_DATA ecm_data;

ECM_DATA::ECM_DATA()
{
	density = 0.5;//arbitrary value for now
	ECM_orientation.resize(3,0.0);//more arbitrary values for now
	anisotropy = 0.0;//same as above
	return;
}

ECM::ECM()
{
	pMicroenvironment = NULL; 
	return;
}

void ECM::sync_to_BioFVM(void)
{
	if( pMicroenvironment == NULL )
	{ pMicroenvironment = get_default_microenvironment(); }
	
	int Xnodes = pMicroenvironment->mesh.x_coordinates.size(); 
	int Ynodes = pMicroenvironment->mesh.y_coordinates.size();  
	int Znodes = pMicroenvironment->mesh.z_coordinates.size();
	
	mesh.resize( default_microenvironment_options.X_range[0] , default_microenvironment_options.X_range[1] ,
		default_microenvironment_options.Y_range[0] , default_microenvironment_options.Y_range[1] ,
		default_microenvironment_options.Z_range[0] , default_microenvironment_options.Z_range[1] ,
		Xnodes, Ynodes, Znodes ); 
	mesh.units = default_microenvironment_options.spatial_units; 
	
    // So you don't have to call the default constructor for ECM_DATA to get the vector of ECM_data structures into the ECM object? That would be cool.
    
	//ECM_DATA.resize( mesh.voxels.size() );
	
	return;
}

//void ECM_setup(double numvox)
//{
//    ecm.sync_to_BioFVM();
//    ecm.ecm_data.resize(numvox);
//    
//    return; 
//}

//void cell_update_from_ecm( void ) // Called in main-ecm.cpp - mgiht want to change that.
//{
//    for(int i = 0; i < (*all_cells).size(); i++)
//    {
//        Cell* pCell = (*all_cells)[i];
//
//        // Ensures that leader cells are not influenced by ECM
//
//        if( pCell->type == 1) // This gets rid of changes in cell speed due to ECM for LEADERS - makes them ignore speed changes.
//        {
////            std::cout<<"Migration bias for type 1 "<<pCell->phenotype.motility.migration_bias << std::endl;
//
//
//            return;
//            std::cout<<"didn't return!"<<std::endl;
//        }
//
//        // Modifies the follower cells
////        std::cout<<"Did modifying cell!"<<pCell->type << std::endl;
////        change_bias_ecm(pCell);
//
//        std::cout<<"Migration bias for type 2 "<<pCell->phenotype.motility.migration_bias << std::endl;
//        change_speed_ecm(pCell); // As long as density is set to 0.5, this will have no impact on cell speed.
////        change_migration_bias_vector_ecm(pCell); MOVED TO follower_cell.functions.update_migration_bias = change_migration_bias_vector_ecm; in AMIGOS-INVASION
//    }
//    return;
//}
//
//void change_bias_ecm(Cell* pCell)
//{
//    int ecm_index =  pCell->get_current_voxel_index();
//    pCell->phenotype.motility.migration_bias = ecm.ecm_data[ecm_index].anisotropy;
//
//    std::cout<<"HEY! Hey! hey!!!"<<std::endl;
//
//    return;
//}
    
//void change_speed_ecm(Cell* pCell)
//{
//    double vmax = 1.0;
//    int ecm_index =  pCell->get_current_voxel_index();
//    double density = ecm.ecm_data[ecm_index].density;
//
//    if(pCell->phenotype.motility.is_motile == true)
//    {
//        // Needs min and max stuff
//        pCell->phenotype.motility.migration_speed = (-(4.0)*pow((density-0.5),2.0) + 1.0)*vmax;
//        std::cout<<"Speed "<<<<pCell->phenotype.motility.migration_speed<<std::endl;
//    }
//
//    return;
//}

//void change_migration_bias_vector_ecm(Cell* pCell)
//{
//     //change bias
//
////    std::vector<double> migration bias direction is the 3-D vector giving the cell's preferred
////    direction of motility for biased Brownian motion. If the user modies this vector, they must ensure
////    it is a unit vector :
////    jjmigration bias directionjj = 1: (14)
////    54
////    5. double migration bias (with a value in [0,1]) sets the degree to which cell motility is biased
////    along migration_bias_direction. If 0, then motion is completely Brownian. If 1, it is completely
////    deterministc along the bias direction.
//
//     int ecm_index =  pCell->get_current_voxel_index();
//     double a = ecm.ecm_data[ecm_index].anisotropy;
//     std::vector<double> d = pCell->phenotype.motility.migration_bias_direction;
//     std::vector<double> f = ecm.ecm_data[ecm_index].ECM_orientation;
//     double ddotf = 0.0;
//
//     normalize(d);
//     normalize(f);
//
////     pCell->phenotype.motility.migration_bias = a;
//
//     //change bias direction
//     for( int i=0; i < d.size() ; i++ )
//     {
////         double temp = d[i] * f[i];
//         ddotf += d[i] * f[i];
//     }
//
//     if(ddotf < 0.0)
//     {
//        for( int i=0; i< f.size(); i++)
//        {
//            f[i] *= -1.0;
////            ecm.ecm_data[ecm_index].ECM_orientation[i] = -1.0 * f[i];
//        }
//     }
//
//     for( int i=0; i < d.size() ; i++ )
//     {
//        pCell->phenotype.motility.migration_bias_direction[i] = a * f[i]  + (1.0-a) * d[i];
//     }
//
//    normalize(pCell->phenotype.motility.migration_bias_direction);
//
//     return;
//}

//void change_migration_bias_vector_ecm(Cell* pCell) MOVED TO AMIGOS-invasion file
//{
//    //change bias
//
//    //    std::vector<double> migration bias direction is the 3-D vector giving the cell's preferred
//    //    direction of motility for biased Brownian motion. If the user modies this vector, they must ensure
//    //    it is a unit vector :
//    //    jjmigration bias directionjj = 1: (14)
//    //    54
//    //    5. double migration bias (with a value in [0,1]) sets the degree to which cell motility is biased
//    //    along migration_bias_direction. If 0, then motion is completely Brownian. If 1, it is completely
//    //    deterministc along the bias direction.
//
//    int ecm_index =  pCell->get_current_voxel_index();
//    double a = ecm.ecm_data[ecm_index].anisotropy;
//    std::vector<double> d = normalize(pCell->phenotype.motility.motility_vector); ////// Changed to motility_vector instead of bias - so the blending is of the actual velocity instead of just the direction it wants to go in
////    d = norm(d);
//    std::vector<double> f = ecm.ecm_data[ecm_index].ECM_orientation;
//    double ddotf = 0.0;
//
//    normalize(d);
//    normalize(f);
//
//    //     pCell->phenotype.motility.migration_bias = a;
//
//    //change bias direction
//    for( int i=0; i < d.size() ; i++ )
//    {
//        //         double temp = d[i] * f[i];
//        ddotf += d[i] * f[i];
//    }
//
//    if(ddotf < 0.0)
//    {
//        for( int i=0; i< f.size(); i++)
//        {
//            f[i] *= -1.0;
//            //            ecm.ecm_data[ecm_index].ECM_orientation[i] = -1.0 * f[i];
//        }
//    }
//
//    for( int i=0; i < d.size() ; i++ )
//    {
//        pCell->phenotype.motility.migration_bias_direction[i] = a * f[i]  + (1.0-a) * d[i];
//    }
//
//    normalize(pCell->phenotype.motility.migration_bias_direction);
//
//    return;
//}
    
void write_ECM_Data_matlab( std::string filename )

{

    int number_of_data_entries = microenvironment.number_of_voxels();

    int size_of_each_datum = 8;

    FILE* fp = write_matlab_header( size_of_each_datum, number_of_data_entries,  filename, "ECM_Data" );  // Note - the size of datum needs to correspond exaectly to the lines of output or there is an error upon importing.

    for( int i=0; i < number_of_data_entries ; i++ )

    {

	    fwrite( (char*) &( ecm.mesh.voxels[i].center[0] ) , sizeof(double) , 1 , fp ); // 1

        fwrite( (char*) &( ecm.mesh.voxels[i].center[1] ) , sizeof(double) , 1 , fp ); // 2

        fwrite( (char*) &( ecm.mesh.voxels[i].center[2] ) , sizeof(double) , 1 , fp ); //3
		
		fwrite( (char*) &( ecm.ecm_data[i].anisotropy), sizeof(double) , 1 , fp ); // 4
	
        fwrite( (char*) &( ecm.ecm_data[i].density), sizeof(double) , 1 , fp ); // 5

        fwrite( (char*) &( ecm.ecm_data[i].ECM_orientation[0]), sizeof(double) , 1 , fp ); // 6

        fwrite( (char*) &( ecm.ecm_data[i].ECM_orientation[1]), sizeof(double) , 1 , fp ); // 7

        fwrite( (char*) &( ecm.ecm_data[i].ECM_orientation[2]), sizeof(double) , 1 , fp ); // 8


        // current voxel index of cell

    }



    fclose( fp );



    return;

}

};
