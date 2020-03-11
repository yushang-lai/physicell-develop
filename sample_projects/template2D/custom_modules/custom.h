#include <list>
#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

// any additional cell types (beyond cell_defaults)

extern Cell_Definition CSC;
extern Cell_Definition DCC; 

// temp for testing 
extern int diff_count; 
extern int div_count; 

 
// custom cell phenotype functions could go here 

void differentiation_function( Cell* pCell , Phenotype& phenotype, double dt ); 

// setup functions to help us along 

void create_cell_types( void );
void setup_tissue( void ); 

// set up the BioFVM microenvironment 
void setup_microenvironment( void ); 

// custom pathology coloring function 

std::vector<std::string> my_coloring_function( Cell* );

// radiation therapy functions 

double calculate_radiation_death_probability( Cell* pCell, Phenotype& phenotype, double dt ); 
// void check_for_radiation_therapy( void ); 
// void check_for_radiation_therapy1( std::list<double> myList); 
void apply_radiation_model( void ); 
/* void check_for_radiation_therapy1( std::list<int> myList)
 */

void update_transition_rate_and_differetion_rate_parameters_factorC_based( Cell* pCell, Phenotype& phenotype, double dt );
void update_transition_rate_and_differetion_rate_parameters_factorC_based_DCC( Cell* pCell, Phenotype& phenotype, double dt );

void show_pressure(void);