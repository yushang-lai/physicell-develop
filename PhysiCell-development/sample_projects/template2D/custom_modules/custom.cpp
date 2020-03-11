#include "./custom.h"

// declare cell definitions here 

Cell_Definition CSC;
Cell_Definition DCC; 

// temp for testing 
int diff_count = 0; 
int div_count = 0; 

void differentiation_function( Cell* pCell, Phenotype& phenotype, double dt )
{
	// am I at the start of a cycle?
	if( phenotype.cycle.data.elapsed_time_in_phase < 0.1 )
	{
//		#pragma omp critical 
//		{ div_count++; } 

		if(PhysiCell_globals.current_time > 6)
		{
			if( UniformRandom() <= parameters.doubles( "differentiation_probability" ) )
			{
	//			#pragma omp critical 
	//			{ diff_count++; } 
				pCell->convert_to_cell_definition( DCC ); 
			}
		}
		
	}

	return; 
} 

void dummy_function( Cell* pCell , Phenotype& phenotype , double dt )
{
	// std::cout << phenotype.cycle.data.transition_rate(0,0) << std::endl; 
}

void create_cell_types( void )
{
	// use the same random seed so that future experiments have the 
	// same initial histogram of oncoprotein, even if threading means 
	// that future division and other events are still not identical 
	// for all runs 
	
	SeedRandom( parameters.ints("random_seed") ); // or specify a seed here 
	
	// housekeeping 
	
	initialize_default_cell_definition();
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment ); 
	
	// Name the default cell type 
	
	cell_defaults.type = 0; 
	cell_defaults.name = "tumor cell"; 
	
	// set default cell cycle model 

	live.phase_link(0,0).fixed_duration = parameters.bools("fixed_cycle_duration"); 
	cell_defaults.functions.cycle_model = live; 
	
	// set default_cell_functions; 
	
	cell_defaults.functions.update_phenotype = NULL; 
	
	// needed for a 2-D simulation: 
	
	/* grab code from heterogeneity */ 
	
	cell_defaults.functions.set_orientation = up_orientation; 
	cell_defaults.phenotype.geometry.polarity = 1.0;
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	

				

	// make sure the defaults are self-consistent. 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 

	// set the rate terms in the default phenotype 

	// first find index for a few key variables. 
	int apoptosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Apoptosis" );
	int necrosis_model_index = cell_defaults.phenotype.death.find_death_model_index( "Necrosis" );
	int oxygen_substrate_index = microenvironment.find_density_index( "oxygen" ); 

	int G0G1_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::G0G1_phase );
	int S_index = flow_cytometry_separated_cycle_model.find_phase_index( PhysiCell_constants::S_phase );
	int live_index = live.find_phase_index( PhysiCell_constants::live ); 

	// initially no necrosis 
	cell_defaults.phenotype.death.rates[necrosis_model_index] = 0.0; 
	
	// initially no apoptosis 
	cell_defaults.phenotype.death.rates[apoptosis_model_index] = 0.0; 

	// set oxygen uptake / secretion parameters for the default cell type 
	cell_defaults.phenotype.secretion.uptake_rates[oxygen_substrate_index] = 10; 
	cell_defaults.phenotype.secretion.secretion_rates[oxygen_substrate_index] = 0; 
	cell_defaults.phenotype.secretion.saturation_densities[oxygen_substrate_index] = 38; 
	
	cell_defaults.custom_data.add_variable( "pressure" , "dimensionless", 0.0 ); 
	cell_defaults.custom_data.add_variable( "volume" , "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "transrate" , "dimensionless", 0.0 );

	// add custom data here, if any 
	

	// Now, let's define another cell type. 
	// It's best to just copy the default and modify it. 
	
	// define CSC type 
	
	CSC = cell_defaults; 
	CSC.type = 1;
	// set birth rate 
	CSC.phenotype.cycle.data.transition_rate(live_index,live_index) = parameters.doubles( "CSC_birth_rate" );  
	// set death rate 
	CSC.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "CSC_apoptosis_rate" ); 
	// set the phenotype function to non-NULL
	CSC.functions.update_phenotype = update_transition_rate_and_differetion_rate_parameters_factorC_based; 

	CSC.phenotype.motility.is_motile = true; 
	CSC.phenotype.motility.migration_speed = 2; 
	CSC.phenotype.motility.migration_bias = 0;
	
	// define DCC type  
	
	DCC = cell_defaults; 
	DCC.type = 2; 	
	// set birth rate 
	DCC.phenotype.cycle.data.transition_rate(live_index,live_index) = parameters.doubles( "DCC_birth_rate" );  
	// set death rate 
	DCC.phenotype.death.rates[apoptosis_model_index] = parameters.doubles( "DCC_apoptosis_rate" ); 

	DCC.phenotype.motility.is_motile = true; 
	DCC.phenotype.motility.migration_speed = 1; 
	DCC.phenotype.motility.migration_bias = 0;

	int factorC_substrate_index = microenvironment.find_density_index( "factorC" );
	DCC.phenotype.secretion.uptake_rates[ factorC_substrate_index ] = 0.0;
	DCC.phenotype.secretion.secretion_rates[ factorC_substrate_index ] = 0.1; //0.0000001;
	DCC.phenotype.secretion.saturation_densities[ factorC_substrate_index ] = 1;
	DCC.functions.update_phenotype = update_transition_rate_and_differetion_rate_parameters_factorC_based_DCC;
	// DCC.phenotype.secretion.saturation_densities[ factorC_substrate_index ] = 1.0;
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters 
	
/* now this is in XML 
	default_microenvironment_options.X_range = {-1000, 1000}; 
	default_microenvironment_options.Y_range = {-1000, 1000}; 
	default_microenvironment_options.simulate_2D = true; 
*/
	
	// make sure to override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	
/* now this is in XML 	
	// no gradients need for this example 

	default_microenvironment_options.calculate_gradients = false; 
	
	// set Dirichlet conditions 

	default_microenvironment_options.outer_Dirichlet_conditions = true;
	
	// if there are more substrates, resize accordingly 
	std::vector<double> bc_vector( 1 , 38.0 ); // 5% o2
	default_microenvironment_options.Dirichlet_condition_vector = bc_vector;
	
	// set initial conditions 
	default_microenvironment_options.initial_condition_vector = { 38.0 }; 
*/
	
	// put any custom code to set non-homogeneous initial conditions or 
	// extra Dirichlet nodes here. 
	
	// initialize BioFVM 
	
	initialize_microenvironment(); 	
	
	return; 
}

void setup_tissue( void )
{
	// create some cells near the origin
	
	Cell* pC;
	
	double Lx = default_microenvironment_options.X_range[1] - default_microenvironment_options.X_range[0]; 
	double Ly = default_microenvironment_options.Y_range[1] - default_microenvironment_options.Y_range[0]; 
	
	for( int n = 0 ; n < parameters.ints( "number_of_stem_cells" ) ; n++ )
	{
		double x = default_microenvironment_options.X_range[0] + 0.1*Lx + 0.8*UniformRandom()*Lx; 
		double y = default_microenvironment_options.Y_range[0] + 0.1*Ly + 0.8*UniformRandom()*Ly; 
		
		pC = create_cell( CSC ); 
		pC->assign_position( x,y, 0.0 );
		
	}
	
	for( int n = 0 ; n < parameters.ints( "number_of_differentiated_cells" ) ; n++ )
	{
		double x = default_microenvironment_options.X_range[0] + 0.1*Lx + 0.8*UniformRandom()*Lx; 
		double y = default_microenvironment_options.Y_range[0] + 0.1*Ly + 0.8*UniformRandom()*Ly; 
		
		pC = create_cell( DCC ); 
		pC->assign_position( x,y, 0.0 );
		
	}
	
	return; 
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
	
	std::vector<std::string> output = { "black", "black", "black", "black" } ; // simple_cell_coloring(pCell); 
		
	if( pCell->phenotype.death.dead == false && pCell->type == 1 )
	{
		 output[0] = "red"; 
		 output[2] = "red"; 
		 return output; 
	}
	
	if( pCell->phenotype.death.dead == false && pCell->type == 2 )
	{
		 output[0] = "cyan"; 
		 output[2] = "cyan"; 
		 return output; 
	}
	
	if( pCell->phenotype.death.dead == true )
	{
		 output[0] = "black"; 
		 output[2] = "black"; 
		 return output; 
	}
	return output; 
}

double calculate_radiation_death_probability( Cell* pCell, Phenotype& phenotype, double dt )
{
	static double alpha_CSC = 0.2;
	static double  beta_CSC = 0.05;
	static double alpha_DCC = 0.2;
	static double  beta_DCC = 0.05;
	static double  dose     = 1;
	static double  dedifferention_rate   = 0;
	if(pCell->type == 1)
	{
		double death_prob = 1.0-exp(-alpha_CSC*dose-beta_CSC*dose*dose);
		// std::cout <<
		// std::cout << "death probability is " << death_prob << " ... " << std::endl << std::endl ;
		return death_prob;
	} 

	// hard-coded for now. put more biology here later. 
	if(pCell->type == 2)
	{
		double death_prob = 1.0-exp(-alpha_DCC*dose-beta_DCC*dose*dose);
		// std::cout << "death probability is " << death_prob << " ... " << std::endl << std::endl  ;
		return death_prob;
	} 
	double death_prob = 1.0-exp(-alpha_DCC*dose-beta_DCC*dose*dose);
	death_prob = 1.0-exp(-alpha_DCC*dose-beta_DCC*dose*dose);
	return death_prob;
}

// void check_for_radiation_therapy( void )
// {
// 	static bool therapy_triggered = false;

// 	if( PhysiCell_globals.current_time >= parameters.doubles( "therapy_time" )-0.005 && therapy_triggered == false )
// 	{
// 		std::cout << "Triggering therapy at time " << PhysiCell_globals.current_time << " ... " << std::endl ; 
// 		apply_radiation_model(); 
// 		therapy_triggered = true; 
// 	}
// 	return; 
// }

/* void check_for_radiation_therapy1( std::list<int> myList)
{	

	std::list<int>:: iterator it;
	for (it = myList.begin();it!=myList.end(); ++it)
			{
				if(PhysiCell_globals.current_time == *it)
				{
					apply_radiation_model(); 
					std::cout << "Triggering therapy at time " << PhysiCell_globals.current_time << " ... " << std::endl ;
				}
			}
} */

// void check_for_radiation_therapy1( std::list<double> myList )
// {
// 	std::list<double> myList({360,720,770});
// 	std::list<double> :: iterator it;
// 	double upper_bd = PhysiCell_globals.current_time+0.005;
// 	double lower_bd = PhysiCell_globals.current_time-0.005;
// 	for(it=myList.begin(); it!=myList.end();++it)
// 	/* std::cout << *it<< ","<< std::endl;; */
// 	if(upper_bd > *it && lower_bd<*it)
// 		{
// 			apply_radiation_model(); 
// 			std::cout << "Triggering therapy at time " << PhysiCell_globals.current_time << " ... " << std::endl ;
// 		}
// 	return;	 
// }
// {	

// 	std::list<int>:: iterator it;
// 	for (it = myList.begin();it!=myList.end(); ++it)
// 			{
// 				if(PhysiCell_globals.current_time == *it)
// 				{
// 					apply_radiation_model(); 
// 					std::cout << "Triggering therapy at time " << PhysiCell_globals.current_time << " ... " << std::endl ;
// 				}
// 			}
// } */


void apply_radiation_model( void )
{
	static int apoptosis_model_index = CSC.phenotype.death.find_death_model_index( "apoptosis" );
	
	int number_of_cells = (*all_cells).size(); 
	std::cout << "therapy tested against " << number_of_cells << " cells at time " << PhysiCell_globals.current_time << "!" << std::endl << std::endl; 
	
	#pragma omp parallel for
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		// isolate a cell 
		Cell* pC = (*all_cells)[i]; 
		
		if(parameters.bools("thearpy_dedifferentiation"))
		{
			if(UniformRandom()<0.8 && pC->type==2)
			{
				pC->convert_to_cell_definition(CSC);
			}
		}
		// get its death probability 
		double death_prob = calculate_radiation_death_probability( pC, pC->phenotype, 0.0 ); 

		// if it's alive, apply death based on the probability 
		if( pC->phenotype.death.dead == false && UniformRandom() <= death_prob )
		{
			pC->start_death(apoptosis_model_index);
			pC->functions.update_phenotype = NULL; 
			pC->phenotype.flagged_for_removal =true;
			pC->phenotype.death.dead = true; 
			pC->functions.custom_cell_rule = NULL;
		}
	}
	
	return; 
}

void update_transition_rate_and_differetion_rate_parameters_factorC_based( Cell* pCell, Phenotype& phenotype, double dt )
{
	if(parameters.bools("disable_feedback_division_rate")==false)
	{
		if( phenotype.death.dead == true )
		{ return; }

		static int start_phase_index; // Q_phase_index; 
		static int end_phase_index; // K_phase_index;
	
		static int factorC_substrate_index = pCell->get_microenvironment()->find_density_index( "factorC" );
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		// necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		double factorc = (pCell->nearest_density_vector())[factorC_substrate_index]; // PhysiCell_constants::oxygen_index]; 
		int n = pCell->phenotype.cycle.current_phase_index();
		
		double multiplier = 1.0;
		static double L = 1;
		static int nn = 1;
		multiplier = 1.0;
		multiplier = 1/(1+L*pow(factorc,nn));
		// if(multiplier<0.5)
		// {
	    // 	std::cout << "Warning: multiplier is:" << multiplier <<".\n"<< std::endl;
		// }
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier * parameters.doubles("CSC_birth_rate");
		static int transrate_i = pCell->custom_data.find_variable_index( "transrate" );
		pCell->custom_data[transrate_i] = pCell->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
	}
		// Update necrosis rate
	static int pressure_i = pCell->custom_data.find_variable_index( "pressure" ); 
	static int volume_i = pCell->custom_data.find_variable_index( "volume" );  


		// am I at the start of a cycle?
	if(parameters.bools("contact_inhibition"))
	{	
		static int start_phase_index; // Q_phase_index; 
		static int end_phase_index; // K_phase_index;
		static double pressure_threshold = 1;
		double multiplier_pressure = 1.0;
		multiplier_pressure = (1.0 - pCell->state.simple_pressure)/pressure_threshold;
		static double upper_bound = 1;
		static double lower_bound = 0;
		if(multiplier_pressure>upper_bound)
		{
			multiplier_pressure = 1.0;
		}
		else if(multiplier_pressure<lower_bound)
		{
			multiplier_pressure = 0.0;
		}
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		// necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier_pressure*parameters.doubles("CSC_birth_rate");
	}

	if(parameters.bools("feedback_Number_on_Transition"))
	{	
		static int start_phase_index; 
		static int end_phase_index; 
		static double feedback_multiplier = 1;
		feedback_multiplier = 1/(1+0.0005*(*all_cells).size());
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) *= feedback_multiplier;
	}


	if( phenotype.cycle.data.elapsed_time_in_phase < 0.1 )
	{
//		#pragma omp critical 
//		{ div_count++; } 

		if(PhysiCell_globals.current_time > 6)
		{	
			double multiplier_prob = parameters.doubles( "differentiation_probability" );
			if(parameters.bools("feedback_Number_on_diff_prob"))
			{
				std::cout << "Triggering therapy at time "<< std::endl;
				multiplier_prob = parameters.doubles( "differentiation_probability" )/(1+0.005*(*all_cells).size());
			}
			if( UniformRandom() <= multiplier_prob )
			{
	//			#pragma omp critical 
	//			{ diff_count++; } 
				pCell->convert_to_cell_definition( DCC ); 
			}
		}
		
	}
		// now, update the necrosis rate 
		// pCell->phenotype.death.rates[necrosis_index] = multiplier * parameters.doubles("differentiation_probability"); 
	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	// necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
	end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	pCell->custom_data[pressure_i] = pCell->state.simple_pressure;
	pCell->custom_data[volume_i] = pCell->phenotype.volume.total;
	static int transrate_i = pCell->custom_data.find_variable_index( "transrate" );
	pCell->custom_data[transrate_i] = pCell->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
	return; 

}

void show_pressure(void)
{
	int number_of_cells = (*all_cells).size(); 
	#pragma omp parallel for
	for( int i=0; i < (*all_cells).size(); i++ )
	{
		// isolate a cell 
		Cell* pC = (*all_cells)[i]; 
		std::cout<<"pressure is:" << pC->state.simple_pressure << std::endl << std::endl;
	}
}

void update_transition_rate_and_differetion_rate_parameters_factorC_based_DCC( Cell* pCell, Phenotype& phenotype, double dt )
{
	if(parameters.bools("disable_feedback_division_rate")==false)
	{
		if( phenotype.death.dead == true )
		{ return; }

		static int start_phase_index; // Q_phase_index; 
		static int end_phase_index; // K_phase_index;
	
		static int factorC_substrate_index = pCell->get_microenvironment()->find_density_index( "factorC" );
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		// necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		double factorc = (pCell->nearest_density_vector())[factorC_substrate_index]; // PhysiCell_constants::oxygen_index]; 
		int n = pCell->phenotype.cycle.current_phase_index();
		
		double multiplier = 1.0;
		static double L = 1;
		static int nn = 1;
		multiplier = 1.0;
		multiplier = 1/(1+L*pow(factorc,nn));

		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier * parameters.doubles("CSC_birth_rate");
		static int transrate_i = pCell->custom_data.find_variable_index( "transrate" );
		pCell->custom_data[transrate_i] = pCell->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
	}
	static int pressure_i = pCell->custom_data.find_variable_index( "pressure" ); 
	static int volume_i = pCell->custom_data.find_variable_index( "volume" ); 
	pCell->custom_data[pressure_i] = pCell->state.simple_pressure;
	pCell->custom_data[volume_i] = pCell->phenotype.volume.total;
		// Update necrosis rate

		// am I at the start of a cycle?
		// now, update the necrosis rate 
		// pCell->phenotype.death.rates[necrosis_index] = multiplier * parameters.doubles("differentiation_probability"); 
	// std::cout << "Warning: multiplier is:" << parameters.bools("disable_contact_inhibition") <<".\n"<< std::endl;
	if(parameters.bools("contact_inhibition"))
	{	
		static double pressure_threshold = 1;
		double multiplier_pressure = 1.0;
		multiplier_pressure = (1.0 - pCell->state.simple_pressure)/pressure_threshold;
		static double upper_bound = 1;
		static double lower_bound = 0;
		if(multiplier_pressure>upper_bound)
		{
			multiplier_pressure = 1.0;
		}
		else if(multiplier_pressure<lower_bound)
		{
			multiplier_pressure = 0.0;
		}
		static int start_phase_index; // Q_phase_index; 
		static int end_phase_index; // K_phase_index;
		// std::cout << "Warning: multiplier from curerrent pressure is:" << multiplier_pressure <<".\n"<< std::endl;
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		// necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) = multiplier_pressure*parameters.doubles("DCC_birth_rate");
	}

	if(parameters.bools("feedback_Number_on_Transition"))
	{	
		static int start_phase_index; 
		static int end_phase_index; 
		static double feedback_multiplier = 1;
		feedback_multiplier = 1/(1+0.0005*(*all_cells).size());
		start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live ); 
		end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index) *= feedback_multiplier;
	}

	static int start_phase_index; // Q_phase_index; 
	static int end_phase_index; // K_phase_index;
	start_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
		// necrosis_index = phenotype.death.find_death_model_index( PhysiCell_constants::necrosis_death_model ); 
	end_phase_index = phenotype.cycle.model().find_phase_index( PhysiCell_constants::live );
	pCell->custom_data[pressure_i] = pCell->state.simple_pressure;
	pCell->custom_data[volume_i] = pCell->phenotype.volume.total;
	static int transrate_i = pCell->custom_data.find_variable_index( "transrate" );
	pCell->custom_data[transrate_i] = pCell->phenotype.cycle.data.transition_rate(start_phase_index,end_phase_index);
	return; 

}
