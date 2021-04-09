/*
    Code for the article “Invasive spread in meta-food-webs depends on landscape structure, fertilization and species characteristics”

    Copyright (C) 2020 Johanna Häussler

    Code is an extension of the code from the article “he biggest losers: Habitat isolation deconstructs
	complex food webs from top to bottom” by Remo Ryser, Johanna Häussler & Markus Stark Proceedings B (2019).

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    version 1.0
    last edit: 18.02.2020 by JH

**** Read the README.md for instructions on how to run the code ****

*/


// Include important libraries
#include "invasive_spread.h"

// Basic simulation variables:
int S;                                                                          // total number of species
int S_c;                                                                        // number consumers (animals, predators)
int S_b;                                                                        // number basal species (plants)
int N=2;                                                                        // number of nutrients
                                                         
char* TEND = getenv("TEND_BASH");                                                    
double tend = atof(TEND);  														// time at which simulation runs stop											
char* TEVAL = getenv("TEVAL_BASH");                                             
double teval = atof(TEVAL);														// length of time interval from which to determine mean(B) and CV(B)

double Delta_t=0.01;                                                            // step size in ODE integration (internal step size can be less or larger)
																			    
const double eps_rel= 1e-10;                                                    // relative error tolerance for ODE solver
const double eps_abs = 1e-10;                                                   // absolute error tolerance for ODE solver
const double EXTINCT = 1e-20;                                                   // extinction threshold for biomass densities
const double GLOBAL_EXTINCT = EXTINCT;                                          // global extinction threshold for species biomass densities
const double MIN_MIGRATION = 1e-17;		                                        // minimum migration threshold (should be higher than EXTINCT => otherwise the solver can get quite slow) 

// Parameters that determine the network topology
double zeta_c=6;                                                                // log_10 of the mean of the consumer (animal) body masses
double sigma_c=3;                                                               // width of the distribution of consumer (animal) body masses
double cutoff_c=1e5;                                                            // half relative cutoff for distribution of consumer body masses
double zeta_b=5;                                                                // log_10 of the mean of the basal (plant) body masses
double sigma_b=3;                                                               // width of the distribution of basal (plant) body masses
double cutoff_b=1e5;                                                            // half relative cutoff for distribution of basal body masses

double m_p_min = 0;                                                             // minimal log10 body masses in case of uniform distributions for plant species
double m_a_min = 2;                                                             // minimal log10 body masses in case of uniform distributions for animal (consumer) species
double m_p_max = 6;                                                             // maximal log10 body masses in case of uniform distributions for plant species
double m_a_max = 12;                                                            // maximal log10 body masses in case of uniform distributions for animal (consumer) species

double cutoff=0.01;                                                             // cutoff of the Ricker curve for setting a link between predator and prey
double R_opt=100;                                                               // optimal predator-prey body-mass ratio
double g=2;                                                                     // width of the Ricker curve (higher g-value -> narrower curve)

double f_herbiv = 0.0;                                                          // fraction of species that are strict herbivores
double f_pred = 0.0;                                                            // fraction of species that are strict predators

// Parameters of the functional response
double a_0=15;                                                                  // scaling factor for the attack rate
double a_c;                                                                     // exponent for predator body-mass scaling of the attack rate (value is drawn from distribution)
double mean_a_c = 0.42;                                                         // (mean for the above distribution)
double sigma_a_c = 0.05;                                                        // (standard deviation for the above distribution)
double a_p;                                                                     // exponent for prey body-mass scaling of the attack rate (value is drawn from distribution)
double mean_a_p = 0.19;                                                         // (mean for the above distribution)
double sigma_a_p = 0.04;                                                        // (standard deviation for the above distribution)
double a_plant = 3500;                                                          // constant attack rate coefficient for plants

double h_0=0.4;                                                                 // scaling factor for the handling time
double h_c;                                                                     // exponent for predator body-mass scaling of the handling time (value is drawn from distribution)
double mean_h_c = -0.48;                                                        // (mean for the above distribution)
double sigma_h_c = 0.03;                                                        // (standard deviation for the above distribution)
double h_p;                                                                     // exponent for prey body-mass scaling of the handling time (value is drawn from distribution)
double mean_h_p = -0.66;                                                        // (mean for the above distribution)
double sigma_h_p = 0.02;                                                        // (standard deviation for the above distribution)

double hill;                                                                    // Hill coefficient (value is drawn from distribution)
double mu_hill=1.5;                                                             // (mean for the above distribution)
double sigma_hill=0.2;                                                          // (standard deviation for the above distribution)
double C_interference;                                                          // interference competition (value is drawn from distribution)
double mu_Cint=0.8;                                                             // (mean for the above distribution)
double sigma_Cint=0.2;                                                          // (standard deviation for the above distribution)

double x_resp_a = 0.141;                                                        // intercept of animal respiration rates
double x_resp_p = 0.138;                                                        // intercept of producer respiration rates
double e_p=0.545;                                                               // assimilation efficiency for plant resources
double e_a=0.906;                                                               // assimilation efficiency for animal resources

double B_b=10;                                                                  // intercept initial basal biomass densities
double B_c=10;                                                                  // intercept initial consumer biomass densities

// Parameters of the nutrient model
double C1=1;                                                                    // content of first nutrient in plants
double C2=0.5;                                                                  // content of second nutrient in plants
double D_N=0.25;                                                                // nutrient turnover rate
double K_min=0.1;                                                               // minimal nutrient uptake half saturation density
double K_rel=0.1;                                                               // width of interval for nutrient uptake half saturation densities
double mu_S; 	                                                                // mean nutrient supply concentration (level of nutrient supply)
double sigma_S=2;                                                               // standard deviation of nutrient supply concentration

// Additional parameter
const double mass_crit = 0;                                                     // body mass threshold at which scaling for extinction boundary changes

// Default parameters for loops and flags
double cv_flag;                                                                 // indicates potentially pathological coefficients of variation
int EXT_FLAG;                                                                   // flag for extinction threshold
int VAR_COEFF = 1;                                                              // 1: calculate the coefficient of variation; 
																				// 0: the simulations run a bit faster => does not work properly => mean biomass not calculated! 
int CANNIBALISM = 0;                                                            // 0: remove all cannibalistic links if necessary (just to be on the safe side)
int KILLSWITCH = 1;                                                             // 0: initialize each species on each patch; 1: initialize arbitrary beta-diversity => remove a fraction of species from some patches during initialization
																		
// Basic simulation variables for spatial dynamics:
int Z; 						                                                    // number of patches (=number of subpopulations)
int D;	 					                                                    // dimension of the biotic system (D=S*Z)
int DN; 					                                                    // dimension of the abiotic system (D=N*Z)

const double exp_suc = 0.33;                                                    // exponent for scaling extinction_boundary

// Parameters for plant migration dynamics
double emigr_a_S_b;                                                             // max. emigration rate
double mean_emigr_a_S_b=0.1;                                                    // (mean for the above distribution)
double sigma_emigr_a_S_b=0.03;                                                  // (standard deviation for the above distribution)
double emigr_b_S_b;                                                             // shape of curve for emigration

// Parameters for animal migration dynamics
double emigr_a_S_c;                                                             // max. emigration rate
double mean_emigr_a_S_c=0.1;                                                    // (mean for the above distribution)
double sigma_emigr_a_S_c=0.03;                                                  // (standard deviation for the above distribution)
double emigr_b_S_c;                                                             // shape of curve for emigration

// Disperal parameters for RGG
double D_0=1.0;                                                                 // scaling factor for dispersal distances 
double eps = 0.05;                                                              // scaling factor of body mass for maximal dispersal distance d_max_i=D_0*m_i^eps
double theta = 1;                                                               // shape of curve for dispersal success

// Parameters and variables imported from bash script
char* web = getenv("WEB");                                                      // web id read in from bash
int webi = atoi(web);                                                           // converting string to integer	 
char* landscape = getenv("LANDSCAPE");                                          // read in patch number variable from bash
int landscapei = atoi(landscape);                                         		// converting string to integer	 
char* directory = getenv("DIR"); 
char* output_dir = getenv("OUTPUT_DIR"); 
char* invasion_dir = getenv("INVASION_DIR"); 
char* bash_seed = getenv("SEED");                                               // read in seed from bash
int seed = atoi(bash_seed);                                                     // different replicates need different random numbers!
char* bash_remspp = getenv("REMSPP");                                           // read in invasive species from bash 
int remspp = atoi(bash_remspp); 	                                            // converting string to integer	 
char* bash_invasion = getenv("INVASION_BASH");                                       // read in INVASION flag from bash => if TRUE S+1 species else S species 
int INVASION = atoi(bash_invasion); 	                                        // converting string to integer	 
char* bash_nutsupply = getenv("NUTSUPPLY");                                     // read in nutrient supply concentration from bash 
double nutsupply = atof(bash_nutsupply); 										// converting string to double 	

char* timeseries = getenv("TIMESERIES_BASH");                                   // read in patch number variable from bash
int TIMESERIES = atoi(timeseries);  
char* joint = getenv("JOINT_BASH");                                             // read in patch number variable from bash
int JOINT = atoi(joint);  
char* isolated = getenv("ISOLATED_BASH");                                        // read in patch number variable from bash
int ISOLATED = atoi(isolated); 

// Inialization of timeseries
char* timeseries_dir = getenv("TS_DIR"); 
char* timeseries_file = getenv("TS_FILE"); 

double counting_migration=0; 
double counting_migration2=0; 
double counting_migration3=0; 
double t_output=0; 

//***********************************
// Main function
//***********************************
/*
The main function first initializes the random number generator (GSL library).
Input Species number can be received from already determined food webs or is set in the code.
There is a difference between conumser species (Sc) and plant/autotroph species (Sb).
The RGG is determined by the number of habitable patches/communities.
Call of the next function (web_calc) which starts the simulation.
*/
int main(int argc, char* argv[])
{
    int i;

    gsl_rng_default_seed = seed;                                                // seed the rng (random number generator)
    gsl_rng *r=gsl_rng_alloc(gsl_rng_default);                                  // initialize the random number generator with the variable r

    web_calc(r);					                                            // initialise model

    gsl_rng_free(r);

    return(0);
}

//****************************************************************************************************************************************************************
// Generate a food web and RGG (if not imported), set food web and spatial parameters, simulate feeding and dispersal dynamics, and save the output 
//****************************************************************************************************************************************************************
/*
Assign local memory for maxtrices and vectors used in the simulations.
Call each function which is necessary to run simulations:
1. pdef_structure - Generate a food web structure from body masses for the simulation
2. set_parameters - Save food web parameters in vector params to use later in simulations of population dynamics
3. Generate_RGG_structure - Set up the spatial environment for simulations including locations of patches and dispersal restrictions for each species
4. set_spatial_parameters - Add spatial parameters to the params vector
5. initialise_biomass_densities - Initialize biomass densities for each species on each patch (optional to start with different subsets of species per patch)
6. solve_ode - solver is applied for the integration of ordinary differential equations and mean biomasses, coefficient of variation and variabilities are calculated continously
7. output - Files are preapred and written to output files for posthoc analysis
Free memory of matrices and vectors
*/
double *web_calc(gsl_rng *r)
{
	
    if(INVASION)
		getInputInvasionSpp();													// import number of basal and consumer species including (S_b,S_c) (with invasive species) 
    else
        getInputSpp();															// import number of basal and consumer species (S_b,S_c) (without invasive species) 

    getInputZ();															// import number of patches (Z)

	S = S_b+S_c;                                                 			// total number of species (S) 
	mu_S = nutsupply; 												    	// nutrient supply concentration
			
    D = S*Z;								                                // dimension of the biotic system
    DN = N*Z;								                                // dimension of the abiotic system
  
    // initialize vectors and matrices for calculations
    gsl_matrix *Ap=gsl_matrix_calloc(S,S);                                  // adjacency matrix 
    gsl_matrix *SW= gsl_matrix_calloc(S*Z,Z);		     	                // success-weighted dispersal matrix
    gsl_matrix *CovM= gsl_matrix_calloc((S+1)*Z, Z);		                // covariance matrix: submarix for each species S (ZxZ) + one for all species (ZxZ)
	
	gsl_vector *mass=gsl_vector_calloc(S);                                  // mean body masses of the species
	gsl_vector *InvSppvec=gsl_vector_calloc(S);                             // indicator for invasive species 
    
    gsl_vector *D_Max = gsl_vector_calloc(S);                                   // maximal dispersal distance for each species
    gsl_vector *Loc= gsl_vector_calloc(2*Z);		     		            // XY patch locations (X1,Y1,X2,Y2,...XZ,YZ)
    gsl_vector *con=gsl_vector_calloc(S);                                   // spatial connectance of each species

    double *B=(double *) calloc(D+DN,sizeof(double));          		        // biomass densities
    double *iniB=(double *) calloc(D+DN,sizeof(double));       		        // initial biomass denisties
    double *meanB=(double *) calloc(D+DN,sizeof(double));      		        // mean biomass densities
    double *meanB_tot=(double *) calloc(S+N,sizeof(double));   		        // mean total biomass densities
    double *VarCoeffB=(double *) calloc(D+DN,sizeof(double));  		        // coefficients of variation
    double *VarCoeffB_tot=(double *) calloc(S+N,sizeof(double)); 	        // coefficients of variation
    double *net_growth_Temp=(double *) calloc(3*D,sizeof(double));		    // values for the net growth rate
    double *Biomass_Temp= (double *) calloc(3*(D+DN),sizeof(double));       // biomass densities of three time steps
    double *RGG_params=(double *) calloc(10,sizeof(double));                // RGG-parameters (spatial connectance, mean distance, etc.)
    int paramslength = 2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+8+D+DN+S;          // length of 'params' array 
    double *params=(double *) calloc(paramslength,sizeof(double));          // 'params' array with all parameters for the local population and spatial dynamics (required input for Cvode_rates)
    
    // start the simulation         
    printf("Timeseries %d\n", TIMESERIES); 
    printf("Joint %d\n", JOINT); 
    printf("Isolated %d\n", ISOLATED); 
    printf("Invasion %d\n", INVASION); 
	
	if(INVASION)
		printf("1 Start simulation for web %d with S+1 = %d species (with invasive species %d).\nMean nutrient supply concentration is %.7g\n",webi,S,remspp,mu_S);
	else 
		printf("1 Start simulation for web %d with S = %d species (without species %d).\nMean nutrient supply concentration is %.7g\n",webi,S,remspp,mu_S);

    pdef_structure(r,Ap,mass);                                        			// import body masses and generate food web structure 
    printf("2 Foodweb %s generated with S = %d; S_b = %d; S_c = %d; N = %d\n",web,S,S_b,S_c,N);

    gsl_rng *ran=gsl_rng_alloc(gsl_rng_default);                                  // initialize the random number generator with the variable r
    set_parameters(ran,Ap,mass,params);                                           // set random parameters for the simulation run
    printf("3 Parameters set\n");
	
    Generate_RGG_structure(mass,D_Max,Loc,SW,con,params,RGG_params);                // immport XY locations for landscapes
    printf("4 Landscape %d generated with Z = %d\n", landscapei, Z);
	
    set_spatial_parameters(SW,Loc,params);                                    // set random spatial parameters for the simulation run
    printf("5 Spatial parameters set\n");

    initialise_biomass_densities(B, iniB, params, r);                   		// initialize biomass densities or if INVASION import t_end biomass densities from simulation run without invasive species
    printf("6 Biomass initialised\n");

    solve_ode(B,meanB,meanB_tot,VarCoeffB,VarCoeffB_tot,params,CovM,Biomass_Temp,net_growth_Temp);     // apply solver for integration of the ODE
    printf("7 ODE solved\n");
    printf("%.2g percent of all migrations were above the threshold\n", counting_migration2/counting_migration3*100);
    printf("%.2g percent of all migrations were below the threshold\n", counting_migration/counting_migration3*100);
    printf("%.g total number of migrations\n", counting_migration3);
     
    output(SW,mass,con,iniB,meanB,meanB_tot,VarCoeffB,VarCoeffB_tot,params,paramslength,Biomass_Temp,RGG_params,net_growth_Temp,InvSppvec,B);   // generate output files for evaluation
    printf("8 Output generated\n");

//  *********** free memories ***********
    gsl_matrix_free(Ap);
    gsl_matrix_free(SW);
    gsl_matrix_free(CovM);
    gsl_vector_free(mass);
    gsl_vector_free(InvSppvec);
    gsl_vector_free(con);
    gsl_vector_free(Loc);
    gsl_vector_free(D_Max);
    gsl_rng_free(ran);
    free(params);
    free(B);
    free(iniB);
    free(meanB);
    free(meanB_tot);
    free(net_growth_Temp);
    free(VarCoeffB);
    free(VarCoeffB_tot);
    free(Biomass_Temp);

    return 0;
}

//********************************************************************
// Generate (or import) the body-mass based food web structure  
//********************************************************************
/*
Body masses are either imported or drawn from a normal or a uniform distribution within certain boundaries.
Feeding links are generated following allometric rules.
*/
static void pdef_structure(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass)
{
	int i,j,flag=0;
    double Wkeit,sigma_i,zeta_act;
    double temp1,temp2,m_opt,R;
    double m_crit;
          
    while(flag==0)
    {
        flag=1;

        gsl_matrix_set_zero(Ap);
		
		if(INVASION)
			getInputInvasionBodymass(mass); 
		else 	
			getInputMass(mass);     										// import body masses of S species (without invasive species)  => determined from a uniform distribution
                                                    
//  ********** fill the adjacency matrix with Ricker feeding links (interaction matrix)  *************
        for(i=0; i<S_c; i++)
        {
            temp1=gsl_vector_get(mass,S_b+i);

            for(j=0; j<S; j++)
            {
                temp2=gsl_vector_get(mass,j);
                R = temp1/temp2;
                if(pow((R / R_opt) * exp(1 - (R / R_opt)), g) >= cutoff)
                    gsl_matrix_set(Ap,S_b+i,j,1);
            }

            if(!CANNIBALISM)
                gsl_matrix_set(Ap,S_b+i,S_b+i,0);                               // remove canibalistic links if necessary

            gsl_vector_view tempp=gsl_matrix_row(Ap,S_b+i);                     // reject networks with consumers or predators without prey
            flag=flag*(1-gsl_vector_isnull(&tempp.vector));
        }

        for(i=0; i<S_b; i++)
        {
            gsl_vector_view tempp = gsl_matrix_column(Ap,i);                    // reject networks with uncontrolled basal species
            flag = flag*(1-gsl_vector_isnull(&tempp.vector));
       }

    }
    
    return;
}

//********************************************************************************************************
// Set and write all parameters required for the local population dynamics to the array 'params'
//********************************************************************************************************
/*
Set random parameters for population dynamics and save them in a vector called params.
Those are for instance attack rates, handling times, respiration rates, assimilation efficiencies, etc.
*/
static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, double params[])
{
    int i,j,k;
    double temp1,temp2,temp3,temp4,temp5,R;

    gsl_matrix *A=gsl_matrix_calloc(S,S);                                       // attack rates
    gsl_matrix *H=gsl_matrix_calloc(S,S);                                       // handling times

    gsl_matrix_memcpy(A,Ap);

    hill = get_random_parameter(mu_hill,sigma_hill,1,2,r);                      // hill coefficient
    C_interference = get_random_parameter(mu_Cint,sigma_Cint,mu_Cint - 3*sigma_Cint, mu_Cint + 3*sigma_Cint,r);     // interference competition

    a_c = get_random_parameter(mean_a_c, sigma_a_c, mean_a_c - 3*sigma_a_c, mean_a_c + 3*sigma_a_c,r);
    a_p = get_random_parameter(mean_a_p, sigma_a_p, mean_a_p - 3*sigma_a_p, mean_a_p + 3*sigma_a_p,r);

    h_c = get_random_parameter(mean_h_c, sigma_h_c, mean_h_c - 3*sigma_h_c, mean_h_c + 3*sigma_h_c,r);
    h_p = get_random_parameter(mean_h_p, sigma_h_p, mean_h_p - 3*sigma_h_p, mean_h_p + 3*sigma_h_p,r);

    for(i=0; i<S; i++)
    {
        if(i >= S_b)       														// the following lines are only for consumer species
        {
            gsl_vector_view tempp=gsl_matrix_row(A,i);
            temp1=gsl_blas_dasum(&tempp.vector);                                // temp1 stores the number of prey species of predator i
            gsl_vector_scale(&tempp.vector,1/temp1);                            // reduce attack rates for generalists
            for(j=0; j<S; j++)
            {
                if(j < S_b)
                {
                    temp1=pow(gsl_vector_get(mass,i),a_p);
                    temp2 = a_plant;
                }
                else
                {
                    temp1 = pow(gsl_vector_get(mass,i),a_c);
                    temp2 = a_0 * pow(gsl_vector_get(mass,j),a_c);
                }
                
                R = gsl_vector_get(mass,i) / gsl_vector_get(mass,j);            // predator-prey body mass ratio for Ricker curve
                temp3 = temp2 * temp1 * pow((R / R_opt) * exp(1 - (R / R_opt)),g);
                gsl_matrix_set(A,i,j,temp3*gsl_matrix_get(A,i,j));          // attack rates
                
                temp4 = h_0 * pow(gsl_vector_get(mass,i),h_c) * pow(gsl_vector_get(mass,j),h_p);
                gsl_matrix_set(H,i,j,temp4);                                    // handling times
                
            }
        }

        for(j=0; j<S; j++)
        {
            *(params+i*S+j)=gsl_matrix_get(A,i,j);                              // save attack rates in 'params'
            *(params+S*S+i*S+j)=gsl_matrix_get(H,i,j);                          // save handling times in 'params'
        }

        if(i < S_b)
        {
            *(params+2*S*S+i)=x_resp_p*pow(gsl_vector_get(mass,i),-0.25);       // plant respiration rates
            *(params+2*S*S+2*S+i)=1;                                            // identifyer for basal species
            *(params+2*S*S+5*S+i) = e_p;                                        // assimilation efficiency for plant resources
        }
        else
        {
            *(params+2*S*S+i)=x_resp_a*pow(gsl_vector_get(mass,i),-0.305);      // animal respiration rates
            *(params+2*S*S+5*S+i) = e_a;                                        // assimilation efficiency for animal resources
        }

        *(params+2*S*S+S+i)=C_interference;                                     // interference competition
        *(params+2*S*S+3*S+i)=hill;                                             // Hill coefficient
        *(params+2*S*S+4*S+i)=gsl_vector_get(mass,i);                           // body masses
    }
	
	if(!INVASION){
    for(i=0; i<S_b; i++)                                                        // the following lines are only for basal species
    {
        for(j=0; j<N; j++){
            temp1=K_min+K_rel*gsl_rng_uniform(r);           // nutrient uptake half saturation densities
            *(params+2*S*S+6*S+i*N+j)=temp1;           // nutrient uptake half saturation densities
		}
    }
	}
	
	if(INVASION){
    for(i=0; i<S_b; i++)                                                        // the following lines are only for basal species
    {
        for(j=0; j<N; j++){
            if(i==remspp) temp1=0;
            else temp1=K_min+K_rel*gsl_rng_uniform(r);           // nutrient uptake half saturation densities
            *(params+2*S*S+6*S+i*N+j)=temp1;           // nutrient uptake half saturation densities
		}
    }
    
	for(j=0; j<N; j++){
		temp1=K_min+K_rel*gsl_rng_uniform(r);						// nutrient uptake half saturation densities invasive species
		*(params+2*S*S+6*S+remspp*N+j)=temp1; 
	}
	}
	


    for(i=0; i<N; i++)
    {
        temp1=mu_S; 
        if(temp1<=0)
            i--;
        else
            *(params+2*S*S+6*S+S_b*N+i)=temp1;                                  // nutrient supply concentrations
    }


    for(i=0; i<S_b; i++)
        *(params+2*S*S+6*S+S_b*N+N+i)=pow(gsl_vector_get(mass,i),-0.25);        // max. nutrient uptake rates

    *(params+2*S*S+6*S+S_b*N+N+S_b) = C1;
    *(params+2*S*S+6*S+S_b*N+N+S_b+1) = C2;


    for(i=0; i<S; i++)
    {
        if(gsl_vector_get(mass,i) < mass_crit)
            *(params+2*S*S+6*S+S_b*N+N+S_b+2+i) = EXTINCT*pow(gsl_vector_get(mass,i),exp_suc);      // extinction boundary for small animals
        else
            *(params+2*S*S+6*S+S_b*N+N+S_b+2+i) = EXTINCT*gsl_vector_get(mass,i);                   // extinction boundary for large animals

    }
   
    gsl_matrix_free(A);
    gsl_matrix_free(H);

}

//***********************************************************************************************
// Set and write all parameters required for the spatial dynamics to the array 'params' 
//***********************************************************************************************
/*
Set random parameters for spatial dynamics and save them in the pararms vector.
Those are for instance dispersal success of each species, maximal emigration rates, location of the patches, etc.
*/
static void set_spatial_parameters(gsl_matrix *SW, gsl_vector *Loc, double params[])
{

    int i=0;
    for(int j=0; j<S*Z; j++)
    {
        for(int k=0; k<Z; k++)
        {
            *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+i) = gsl_matrix_get(SW, j, k);   //success-weighted dispersal matrix
            i++;
        }
    }

    for(int j=0; j<Z; j++)
    {
        *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+0+2*j) = gsl_vector_get(Loc, 0+j*2);  // x-coordinates
        *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+1+2*j) = gsl_vector_get(Loc, 1+j*2);  // y-coordinates
    }
    
    emigr_a_S_c = 0.1; // max. emigration rate consumers
    emigr_b_S_c = 10; // shape of emgiration function consumers

    emigr_a_S_b = 0.1; // max. emigration rate plants
    emigr_b_S_b = 10; // shape of emgiration function plants

    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z) = emigr_a_S_c;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+1) = emigr_b_S_c;

    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+2) = emigr_a_S_b;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+3) = emigr_b_S_b;

    // parameters for the RGG
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+4) = eps;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+5) = D_0;
    *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+6) = theta;

    // parameters for net growth rate 
    for(int j=0; j<D; j++)
        *(params+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+7+j) = 0;
    
    return;
}

//******************************************************************************
// Initialize biomass densities and nutrient concentrations 
//******************************************************************************
/*
1.Initial biomass densities for each patch are drawn from a uniform distribution and written into the biomass array B[] and saved in iniB[].
2. Killswitch: If TRUE, a scenario where an initial ß-diversity is introduced which is implemented in the way that on each patch biomasses of a random 
fraction of species are set to zero.
*/
static void initialise_biomass_densities(double B[], double iniB[], double params[], gsl_rng *r)
{

    int i,j;
    double temp1;
    double *pparams=(double *)params;                                   // pointer to first element of params...

	if(INVASION)
		getInputInvasionBiomass(B); 
	else 
	{
		// Initialize biomass densities for S species on Z patches 
		for(i=0; i<S; i++)
		{
			for(j=0; j<Z; j++)
			{

					if(i<S_b)
						B[j*S+i]=B_b*gsl_rng_uniform_pos(r);                        // basal species
					if(i>=S_b)
						B[j*S+i]=B_c*gsl_rng_uniform_pos(r);                        // consumer species
			}
		}

		// Initialize N nutrients on Z patches
		for(i=0; i<N; i++)
		{
			for(j=0; j<Z; j++)
			{
					temp1=(double)*(params+2*S*S+5*S+2*S_b*N+i);					// define location in params where they are to append 
					B[D+j*N+i] = 0.5*temp1+0.5*temp1*gsl_rng_uniform_pos(r);        // nutrients

			}
		}

		int test_iniB_b = 1; 
		if(!INVASION && KILLSWITCH) // If KILLSWITCH TRUE, generate an arbitrary beta-diversity => set biomass densities for a random fraction of species to 0
		{
			while(test_iniB_b>0){
			
			for(j=0; j<Z; j++)
			{

				int no_species_killed = S*0.4; 										// set number of species to be killed on each patch to 40 % 
				
				double temp3[no_species_killed],temp2[S];                           // initialise arrays for gsl_ran_choose

				for (i = 0; i < S; i++)
					temp2[i] = (double) i;                                          // fill array with species numbers		
				
				gsl_ran_choose(r, temp3, no_species_killed, temp2, S, sizeof (double));    // randomly select species
					
				for(i=0; i<S; i++)
				{
					for(int m=0; m<no_species_killed; m++)
					{
						if(i==temp3[m]){
							B[j*S+i]=0;                                        // initialise biomass with 0
						}
					}
				}
			}
		
			for(j=0; j<Z; j++)													  // testing to make sure that on each patch is at least one basal species initialized
			{	
				double sum_B_basal = 0; 
				for(i=0; i<S_b; i++){
					sum_B_basal = sum_B_basal + B[j*S+i];                        // basal species
				}
				if(sum_B_basal == 0) {
					test_iniB_b = 2;   
					printf("All basal species are extinct %.9g on patch %d\n", sum_B_basal, j);    	
				}
			}
		
			if(!test_iniB_b)
				printf("test_iniB_b %d: on one patch there are no basal species \n", test_iniB_b);    	
			else 
				test_iniB_b = 0; 
		}
		}
	}

    for(i=0; i<D+DN; i++) 								
		iniB[i]=B[i]; 													  // saving initial biomass values
        	
    return;
}


//***************************************************************************************************
// Efficient environment for solving the ODE, including reducing the dimension if possible
//***************************************************************************************************
/*
1. Assign local memory for maxtrices and vectors used in the solving environment & get fixed parameters from params matrix.
2. Initialization of a file in which the timeseries is saved as .txt if TIMESERIES = 1.
3. Initialization of the solving environment from the sundials library.
4. Start of the while loop which ends when t = tend
    a) Always use CV_NORMAL, i.e. the solver uses a fixed step length => otherwise dispersal per time step doesn't make sense;
    b) Optimization of the minimum step length during simulations;
    c) Test whether local biomasses are below an extinction threshold & global biomass of a species is below a threshold - if yes then bioamss is set to zero;
    d) If TIMESERIES = 1, biomasses are printed to the timeseries file;
    e) If simulation time (t) has reached a certain simulation time (tend-2*teval) a first evaluation starts
        - For teval steps of the simulation the biomasses are summed up for each species and each nutrient on each patch
        - After those time steps the mean biomass of each species and each nutrient on each patch is calculated
        - During the last teval time steps the coefficient of variation & variability of the timeseries are calculated following the approach by Wang&Loreau (2014).
5. Free memory of matrices and vectors.
*/
static void solve_ode(double B[], double meanB[], double meanB_tot[], double CV[], double CV_tot[], double params[], gsl_matrix *CovM, double Biomass_Temp[], double net_growth_Temp[])
{
	
	printf("Start the solver...\n"); 
	
    int i,j,m,flag,ERR_FLAG,tempBm;
    double t = 0;                                                       		// starting time
    double *pparams=(double *)params;          									// pointer to first element of params...

    double *meanB2=(double *) calloc(D+DN,sizeof(double));             			// mean biomass densities (for comparison)
    gsl_vector_view Mvec2 = gsl_vector_view_array(meanB2,D+DN);       			// mean biomass denisties (including resources), for comparison

    gsl_vector_view Bvec = gsl_vector_view_array(B,D+N*Z);              		// biomass densities (including resources)
    gsl_vector_view Mvec = gsl_vector_view_array(meanB,D+N*Z);          		// mean biomass denisties (including resources)
    gsl_vector_view MTvec = gsl_vector_view_array(meanB_tot,S);         		// mean total biomass densities
    gsl_vector_view netvec = gsl_vector_view_array(net_growth_Temp,3*D);        // net growth rates temp 1&2 and their Standard deviation from 2*D onwards
    gsl_vector_view MTNvec = gsl_vector_view_array_with_stride(meanB_tot,S,N);  // mean total nutrient densities
    gsl_vector_view Bvec_temp = gsl_vector_view_array(Biomass_Temp,3*(D+N*Z));  // biomass densities (including resources)

    gsl_vector *Mvec_all = gsl_vector_calloc(Z);								// aggregate mean biomasses of all species per patch
    gsl_vector *Nvec_all = gsl_vector_calloc(Z);								// aggregate mean nutrient level per patch

    gsl_vector_view net_growthrate = gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+7,D); // net growth rate
    gsl_vector *netgrowthrate = &net_growthrate.vector;

    gsl_vector_view net_temp1 = gsl_vector_subvector(&netvec.vector,0,D);       // net growth rates temporally summed
    gsl_vector_view net_temp2 = gsl_vector_subvector(&netvec.vector,D,D);       // net growth rates temporally summed

    // If TRUE initialize timeseries file 
    FILE *timeseries;
    if(TIMESERIES)
        Prepare_timeseries_file(timeseries);

    double tout = 0.0001; //0.01; //Delta_t;                                             // initial step size

    // ODE Solver environment
    N_Vector y = N_VNew_Serial(D+DN);
    N_VSetArrayPointer(B, y);

    void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
    CVodeSetUserData(cvode_mem,params);
    CVodeInit(cvode_mem, Cvode_rates, t, y);
    CVodeSStolerances(cvode_mem, eps_rel, eps_abs);
    CVDense(cvode_mem, D+DN);

// ********* integrate ODE to reach an attractor; optional: output of time series ********
    ERR_FLAG = 0;
    while(t < tend)
    {
        flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);                    // Solver with fixed step length

        if(flag != CV_SUCCESS)
        {
            ERR_FLAG = 1;
            break;
        }
 
        else if(tout < Delta_t)                                                 // set minimum step length of tout
            tout *= 10;
        else
            tout += Delta_t;

        // Set biomass to zero below a certain threshold
        EXT_FLAG = 0;
        for (i=0; i<D+DN; i++)
        {
            if (B[i] < EXTINCT && B[i] != 0)                     		        // determine whether a species just fell below the extinction threshold
                EXT_FLAG = 1;
            }
        if (EXT_FLAG)                                          			        // if an extinction occured:
        {
            Extinct_Species(B,S,Z);									            // set the biomass of a whole species to 0 if the total biomass is below a certain threshold
            Extinction(B,D+DN);                                     	        // set the biomass of the species to 0
            tout = t+0.01;                                       		        // set the next output time to a reasonable value
            CVodeReInit(cvode_mem, t, y);                     // and re-initialise the solver
            CVDense(cvode_mem, D+DN);

        }

        // If TRUE write timeseries to a file
        if(TIMESERIES)
            Write_timeseries_to_file(B,t,timeseries);

        else
        // Evaluation of timeseries => calculate the coefficient of variation (only if the step size is constant)
        {

            if(fabs(t - (tend - 2*teval)) < 0.1*Delta_t)
            // 	Copy the current biomasses to a vector for evaluation (tend-20.000)
            {
                gsl_vector_view Bvec_Temp1 = gsl_vector_subvector(&Bvec_temp.vector,0,D+N*Z);
                gsl_vector_memcpy(&Bvec_Temp1.vector,&Bvec.vector);
				t_output = t; 
            }

            if(t >= tend - (1+VAR_COEFF)*teval && t < tend - teval)
            // Sum up the biomasses (&net growth rates) for a certain time period
            {
                gsl_vector_add(&Mvec.vector,&Bvec.vector);
                gsl_vector_add(&net_temp1.vector,netgrowthrate);
                t_output = t; 
            }

            if(fabs(t - (tend - teval)) < 0.1*Delta_t)
            // Calculate mean biomasses
            {
                gsl_vector_scale(&Mvec.vector,Delta_t/teval);					// scale the biomasses to one time step
                gsl_vector_scale(&net_temp1.vector,Delta_t/teval);				// scale the net growth rate to one time step

                gsl_vector_view M1vec = gsl_vector_subvector(&Mvec.vector,0,S); // mean biomass densities in first patch
                gsl_vector_view M2vec = gsl_vector_subvector(&Mvec.vector,S,S); // mean biomass densities in second patch

                gsl_vector_memcpy(&MTvec.vector,&M2vec.vector);
                gsl_vector_add(&MTvec.vector,&M1vec.vector);                    // total biomasses

                gsl_vector_view M1Nvec = gsl_vector_subvector_with_stride(&Mvec.vector,0,D,N);   // mean nutrient densities in first patch
                gsl_vector_view M2Nvec = gsl_vector_subvector_with_stride(&Mvec.vector,D,N,N);   // mean nutrient densities in second patch

                gsl_vector_memcpy(&MTNvec.vector,&M2Nvec.vector);
                gsl_vector_add(&MTNvec.vector,&M1Nvec.vector);                  // total nutrient densities

                // 	Copy the current biomasses to a vector for evaluation (tend-10.000)
                gsl_vector_view Bvec_Temp2 = gsl_vector_subvector(&Bvec_temp.vector,D+N*Z,D+N*Z);
                gsl_vector_memcpy(&Bvec_Temp2.vector,&Bvec.vector);
				t_output = t; 

                for(m=0; m<Z; m++)
                {

                    gsl_vector_view Mvec_sp = gsl_vector_subvector(&Mvec.vector,m*(S+N),S); // mean biomass densities of a single patch without resource
                    gsl_vector_view Mvec_nu = gsl_vector_subvector(&Mvec.vector,S+m*(S+N),N); // mean biomass densities of a single patch without resource

                    tempBm = 0;
                    for(i=0; i<S; i++)
                    {
                        tempBm += gsl_vector_get(&Mvec_sp.vector, i);
                    }
                    gsl_vector_set(Mvec_all, m, tempBm);				        // total biomass of one patch without resource

                    tempBm = 0;
                    for(i=0; i<N; i++)
                    {
                        tempBm += gsl_vector_get(&Mvec_nu.vector, i);
                    }
                    gsl_vector_set(Nvec_all, m, tempBm);				        // total nutrient level of each patch

                }

                gsl_vector_scale(Mvec_all, 1./S);                               // Mean species biomass of every patch
                gsl_vector_scale(Nvec_all, 1./N);                               // Mean nutrient level of every patch

            }

            if(t >= tend - teval)
            //  Calcualte mean biomass a second time
            {
                gsl_vector_add(&Mvec2.vector,&Bvec.vector);                     // compute mean biomass densities second time
                gsl_vector_add(&net_temp2.vector,netgrowthrate);		        // copy the current net growth rates to a vector for evaluation (tend-10.000)
            }

            if(fabs(t - tend) < 0.1*Delta_t)
            {
                gsl_vector_scale(&net_temp2.vector,Delta_t/teval);				// scale the net growth rate to one time step
				t_output = t; 
            }
        }
    }
	
    // 	Copy the current biomasses to a vector for evaluation (tend)
    gsl_vector_view Bvec_Temp3 = gsl_vector_subvector(&Bvec_temp.vector,2*(D+N*Z),D+N*Z); // if solver stops earlier then this will not be at the set tend!
    gsl_vector_memcpy(&Bvec_Temp3.vector,&Bvec.vector);

    gsl_vector_scale(&Mvec2.vector,Delta_t/teval);                              // scale the mean biomass vector to one time step

    double norm1 = gsl_blas_dasum(&Mvec.vector);
    double norm2 = gsl_blas_dasum(&Mvec2.vector);

    cv_flag = norm1 - norm2;

    // ****  Free memories *****
    CVodeFree(&cvode_mem);
    free(meanB2);
    gsl_vector_free(Mvec_all);
    gsl_vector_free(Nvec_all);

    return;
}


//*********************************************************************************************************
// Set biomass density for species i on patch j to 0 if it is below the extinction threshold
//*********************************************************************************************************
/*
Function checks whether the biomass of each species on every patch is below a certain threshold level.
If that is the case the biomass will be set to zero.
*/
static void Extinction(double B[], int Num)
{
    int i;
    for (i=0; i<Num; i++)
    {
        if (B[i] < EXTINCT)
            B[i] = 0;
    }

    return;
}


//***************************************************************************************
// Set biomass density of species i to 0 if species is globally extinct 
//***************************************************************************************
/*
Function checks whether the biomass of a whole species is below a certain threshold level.
If that is the case the biomass of the whole species will be set to zero.
*/
static void Extinct_Species(double B[], int Num, int Pat)
{
    int i,m;
    for (i=0; i<Num; i++)
    {
        double B_Species = 0;                                                   // Count of biomass level for each species

        for (m=0; m<Pat; m++)
        {
            B_Species = B_Species + B[i+Num*m];                                 // Sum up biomass for each species
        }

        if (B_Species < GLOBAL_EXTINCT)                                         // If the sum of biomass is below the extinction threshold, the biomass of the whole species will be set to zero
        {
            for (m=0; m<Pat; m++)
            {
                B[i+Num*m] = 0;
            }
        }
    }

    return;
}

//*******************************************************************
// This function is required by the cvode-solver
//*******************************************************************
/*
Call of the function dynamics
That function calculates alle the local dynamics for each time step
*/
static int Cvode_rates(double t, N_Vector y, N_Vector ydot, void *params)
{
    dynamics(N_VGetArrayPointer(y),N_VGetArrayPointer(ydot),params);
    return(0);
}

//**********************************************************
// Define population dynamics equations
//**********************************************************
/*
1. Define vectors and matrices for calculating population dynamics
2. Calcualation of net growth rates in the function net_rates
3. Immigrating biomasses are calculated for each species on each patch from emigrating biomasses and the dispersal success
4. Minimum Immigration biomass threshold discards very low immigrating biomass levels
5. Add up growth/loss terms and immigrating biomass on each patch
6. Check whether there are some negative biomass levels, if yes => set to zero
7. Free memory of vectors
*/
static void dynamics(double B[], double Bdot[], void *params)
{
    int i,j;
    double *pparams=(double *)params;                                           // pointer to first element of params

    gsl_vector_view TB_vec = gsl_vector_view_array(B,D+DN);                     // vector with all biomass densities and nutrient concentration
    gsl_vector *TBvec = &TB_vec.vector;

    gsl_vector_view Bdot_vec=gsl_vector_view_array(Bdot,D);                     // vector for time derivatives of biomass densities
    gsl_vector *Bdotvec=&Bdot_vec.vector;

    gsl_vector_view Ndot_vec=gsl_vector_view_array(Bdot+D,DN);                  // vector for time derivatives of nutrient concentrations
    gsl_vector *Ndotvec=&Ndot_vec.vector;

    gsl_vector *Ivec=gsl_vector_calloc(D);                                      // ingestion: nutrient uptake, consumption, predation
    gsl_vector *Dvec=gsl_vector_calloc(D);                                      // death (predation losses)
    gsl_vector *Xvec=gsl_vector_calloc(D);                                      // respiration

    gsl_vector *N_in=gsl_vector_calloc(DN);                                     // influx of nutrients
    gsl_vector *N_out=gsl_vector_calloc(DN);                                    // outflux of nutrients
    gsl_vector *N_up=gsl_vector_calloc(DN);                                     // uptake of nutrients by plants

    gsl_vector *Imvec=gsl_vector_calloc(D);                                     // immigration
    gsl_vector *Emvec=gsl_vector_calloc(D);                                     // emigration

    gsl_matrix_view SW_mat=gsl_matrix_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2+S,D,Z);  // dispersal succes weighted matrix
    gsl_matrix *SW =&SW_mat.matrix;

    net_rates(TBvec,Ivec,Dvec,Xvec,Emvec,N_in,N_out,N_up,params);               // calculate net rates

//  ************* calculate dispersal success *************
       
    for(j=0; j<Z; j++)
    {
        for(i=0; i<S; i++)
        {
            gsl_vector_view Im_vec_i = gsl_vector_subvector_with_stride(Imvec, i, S, Z);
            gsl_vector *Imvec_i=&Im_vec_i.vector;

            gsl_vector_view Em_vec_i = gsl_vector_subvector_with_stride(Emvec, i, S, Z);
            gsl_vector *Emvec_i=&Em_vec_i.vector;

            gsl_matrix_view SW_i_mat = gsl_matrix_submatrix(SW,i*Z,0,Z,Z);      // dispersal success of each species 
            gsl_matrix *SW_i = &SW_i_mat.matrix;
				
            // Calcualtion of Immigration
            gsl_blas_dgemv(CblasNoTrans,1,SW_i,Emvec_i,0,Imvec_i);
        }
    }

// Threshold limit for minimum migrating biomasses
    for(i=0; i<D; i++)
    {
        if(gsl_vector_get(Imvec,i)<MIN_MIGRATION && gsl_vector_get(Imvec,i) != 0){
            gsl_vector_set(Imvec,i,0);
            counting_migration=counting_migration+1; 
           }
        else counting_migration2=counting_migration2+1;
        counting_migration3=counting_migration3+1; 
    }
	
//  ************* assmble dynamical equations (Bdot) *************
    gsl_vector_memcpy(Ndotvec,N_in);
    gsl_vector_sub(Ndotvec,N_out);
    gsl_vector_sub(Ndotvec,N_up);

    gsl_vector_memcpy(Bdotvec,Ivec);
    gsl_vector_sub(Bdotvec,Dvec);
    gsl_vector_sub(Bdotvec,Xvec);
    
    if(!JOINT) 															// If JOINT scenario: 100% dispersal success => patches behave as one (no emigration & no immigrationss) 
		gsl_vector_sub(Bdotvec,Emvec); 								   	// just setting the success to 1 makes solver get stuck early in simulation
    
//  ************* control of negative biomass densities *************
    for(i=0; i<D; i++)
    {
        if (B[i] < EXTINCT)
            Bdot[i] = 0;
    }
    
    if(!ISOLATED && !JOINT)												// If ISOLATED scenario: 0% dispersal scucess 
		gsl_vector_add(Bdotvec,Imvec); 									// => only dispersal loss due to emigrating biomass (no immigration) 
	
//  ********** free memory ************
    gsl_vector_free(Ivec);
    gsl_vector_free(Dvec);
    gsl_vector_free(Xvec);
    gsl_vector_free(Emvec);
    gsl_vector_free(Imvec);
    gsl_vector_free(N_in);
    gsl_vector_free(N_out);
    gsl_vector_free(N_up);

    return;
}

//********************************************************************************
// Calculate biomass population rates between species and resources
//********************************************************************************
/*
For the calulations all biomass levels of resources, plant and consumer species are needed (TBvec),
an ingestion term (Ivec), death term (Dvec), respiration term, emigration vector , in- and outflow of nutrients are defined
in the function dynamics and are calculated within this function

1. Get all constant parameters and allocate it in a vector
2. Define vectors and matrices for calculation of population dynamics
3. For all habitat patches population dynamics are calculated
    a) If a species is attributed a Type III functional response, it is calculated
    b) For all other species functional response Type II is calculated
    c) Rate of change for the nutrients is calculated according to a chemostat model
    d) Add up all the growth and loss terms for each species
    e) Calculate an emigration rate depending on the net growth rate
    f) Get emigrating biomasses
4. Free memory of vectors and matrix
*/
static void net_rates(gsl_vector *TBvec, gsl_vector *Ivec, gsl_vector *Dvec, gsl_vector *Xvec, gsl_vector *Emvec, gsl_vector *N_in, gsl_vector *N_out, gsl_vector *N_up, void *params)
{

    int i,p,k;
    double *pparams=(double *)params;                                           // pointer to the first element of params

    gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);                   // a_ij*f_ij - attack rates
    gsl_matrix *Amat=&A_mat.matrix;

    gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);               // H_ij - handling times
    gsl_matrix *Hmat=&H_mat.matrix;

    gsl_vector_view R_vec=gsl_vector_view_array(pparams+2*S*S,S);               // respiration rates
    gsl_vector *Rvec=&R_vec.vector;

    gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);           // interference vector
    gsl_vector *Cvec=&C_vec.vector;

    gsl_vector_view P_vec=gsl_vector_view_array(pparams+2*S*S+2*S,S);           // basal vector (1 for basal species, 0 otherwise)
    gsl_vector *Pvec=&P_vec.vector;

    gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);           // Hill coefficients
    gsl_vector *Hvec=&H_vec.vector;

    gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);           // body masses
    gsl_vector *Mvec=&M_vec.vector;

    gsl_vector_view E_vec=gsl_vector_view_array(pparams+2*S*S+5*S,S);           // assimilation efficiencies
    gsl_vector *Evec=&E_vec.vector;

    gsl_matrix_view K_mat=gsl_matrix_view_array(pparams+2*S*S+6*S,S_b,N);       // nutrient uptake half saturation densities
    gsl_matrix *Kmat=&K_mat.matrix;                                             // (P rows, N columns)

    gsl_vector_view S_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N,N);     // nutrient supply concentrations
    gsl_vector *Svec=&S_vec.vector;

    gsl_vector_view U_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N,S_b); // max. nutrient uptake rates
    gsl_vector *Uvec=&U_vec.vector;

    gsl_vector_view Cc_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b,N); // nutrient contents in plants
    gsl_vector *Ccvec=&Cc_vec.vector;

    gsl_vector_view net_growthrate = gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b+2+S+S*Z*Z+2*Z+7,D); // net growth rate
    gsl_vector *netgrowthrate = &net_growthrate.vector;

    gsl_vector *Prey=gsl_vector_calloc(S);                                      // prey biomasses

    gsl_matrix *mmat=gsl_matrix_calloc(S,S);
    gsl_vector *rvec=gsl_vector_calloc(S);
    gsl_vector *svec=gsl_vector_calloc(S);
    gsl_vector *nvec=gsl_vector_calloc(N);
    gsl_vector *mvec=gsl_vector_calloc(N);
    gsl_vector *pvec=gsl_vector_calloc(S_b);
    gsl_vector *netvec=gsl_vector_calloc(S);								    // net growth rate without migration

    for(p=0; p<Z; p++)
    {

        gsl_vector_view B_vec_p	= gsl_vector_subvector(TBvec, p*S, S); 		    // all biomasses (plants + animals) on patch p
        gsl_vector *Bvec_p = &B_vec_p.vector;

        gsl_vector_view PB_vec_p	= gsl_vector_subvector(TBvec, p*S, S_b);    // basal biomasses (plants) on patch p
        gsl_vector *PBvec_p = &PB_vec_p.vector;

        gsl_vector_view NB_vec_p	= gsl_vector_subvector(TBvec, D+p*N, N);    // nutrient biomasses on patch p
        gsl_vector *NBvec_p = &NB_vec_p.vector;

        gsl_vector_view D_vec_p	= gsl_vector_subvector(Dvec, p*S, S); 		    // death rates on patch p
        gsl_vector *Dvec_p = &D_vec_p.vector;

        gsl_vector_view X_vec_p	= gsl_vector_subvector(Xvec, p*S, S); 		    // respiration rates on patch p
        gsl_vector *Xvec_p = &X_vec_p.vector;

        gsl_vector_view I_vec_p	= gsl_vector_subvector(Ivec, p*S, S); 		    // ingestion: nutrient uptake, consumption, predation on patch p
        gsl_vector *Ivec_p = &I_vec_p.vector;

        gsl_vector_view Em_vec_p = gsl_vector_subvector(Emvec, p*S, S); 	    // emigratio rates on patch p
        gsl_vector *Emvec_p = &Em_vec_p.vector;

        gsl_vector_view N_in_p = gsl_vector_subvector(N_in, p*N, N); 	        // nutrient influx on patch p
        gsl_vector *Nin_p = &N_in_p.vector;

        gsl_vector_view N_out_p = gsl_vector_subvector(N_out, p*N, N); 	        // nutrient outflux on patch p
        gsl_vector *Nout_p = &N_out_p.vector;

        gsl_vector_view N_up_p	= gsl_vector_subvector(N_up, p*N, N); 	        // nutrient uptake on patch p
        gsl_vector *Nup_p = &N_up_p.vector;

        gsl_vector_view net_growthrate_p = gsl_vector_subvector(netgrowthrate,p*S,S); // species-specific subvector for each patch
        gsl_vector *netgrowth_p = &net_growthrate_p.vector;


//  ************ in case of Type 3 functional response ***********

            gsl_vector_memcpy(Prey,Bvec_p);

            for(i=0; i<S; i++)
            {
                if(gsl_vector_get(Prey,i)<0)
                    gsl_vector_set(Prey,i,0);                                   // paranoia (this doesn't mean that you can safely remove these lines!)

                gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));  // filling prey vector (prey biomasses to the power of hill coefficient)
            }

//  ************* calcualate functional responses *************

            gsl_matrix_memcpy(mmat,Hmat);
            gsl_matrix_mul_elements(mmat,Amat);

            gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);                     // handling time term: r_i=sum_j a_ij*h_ij*Prey_j
            gsl_vector_memcpy(svec,Bvec_p);
            gsl_vector_mul(svec,Cvec);                                          // predator interference: s_i=C_i*B_i
            gsl_vector_add(rvec,svec);
            gsl_vector_add_constant(rvec,1);
            gsl_vector_mul(rvec,Mvec);                                          // rvec: denominator of functional response

            gsl_vector_memcpy(svec,Evec);                                       // assimilation efficiency (prey-specific)
            gsl_vector_mul(svec,Prey);
            gsl_blas_dgemv(CblasNoTrans,1,Amat,svec,0,Ivec_p);                  // numerator of prey intake term: I_i = sum_j e_j*A_ij*Prey_j

            gsl_vector_div(Ivec_p,rvec);                                        // divide by denominator...
            gsl_vector_mul(Ivec_p,Bvec_p);                                      // multiply with target species' biomass

            gsl_vector_memcpy(Xvec_p,Rvec);                                     // compute the total respiration rate:
            gsl_vector_mul(Xvec_p,Bvec_p);                                      // per unit biomass respiration times biomass density

            gsl_vector_memcpy(svec,Bvec_p);                                     // predator biomass
            gsl_vector_div(svec,rvec);                                          // ... divide by denominator of the functional response
            gsl_blas_dgemv(CblasTrans,1,Amat,svec,0,Dvec_p);                    // per capita predation losses:
            gsl_vector_mul(Dvec_p,Prey);                                        // multiply with target species' biomass

            //  ************* terms for basal species *************
            gsl_vector_view u_vec=gsl_vector_subvector(rvec,0,S_b);             // basal species growth rates
            gsl_vector *uvec=&u_vec.vector;

            for(i=0; i<S_b; i++)
            {
                gsl_vector_view K_vec=gsl_matrix_row(Kmat,i);
                gsl_vector *Kvec=&K_vec.vector;                                 // Kvec is a vector of length N (Kvec_j=Kmat_ij)

                gsl_vector_memcpy(nvec,NBvec_p);
                gsl_vector_add(nvec,Kvec);
                gsl_vector_memcpy(mvec,NBvec_p);
                gsl_vector_div(mvec,nvec);                                      // m_j=NB_j/(Kvec_j+NB_j)

                gsl_vector_set(uvec,i,gsl_vector_min(mvec));                    // remember: uvec is just a vector_view of rvec!
            }

            gsl_vector_mul(uvec,Uvec);
            gsl_vector_mul(uvec,PBvec_p);                                       // uvec_i=PB_i*mas_i^-0.25*Min_j(NB_j/(NB_j+Kmat_ij))

            gsl_vector_mul(rvec,Pvec);                                          // growth function only for plant species
            gsl_vector_add(Ivec_p,rvec);


//  ********** put the rates for the nutrient dynamics together ********** // nutrient dynamics on non-habitat patches ???
            gsl_vector_memcpy(Nin_p,Svec);
            gsl_vector_scale(Nin_p,D_N);                                        // nutrient influx: D_N * S_i

            gsl_vector_memcpy(Nout_p,NBvec_p);
            gsl_vector_scale(Nout_p,D_N);                                       // nutrient outflux: D_N * N_i

            gsl_vector_memcpy(Nup_p,Ccvec);
            gsl_vector_scale(Nup_p,gsl_blas_dasum(uvec));                       // nutrient uptake: C_i * sum_j (r_j*G_j*P_j)


//  ********** calculate net growth rates for basal and consumer species **********
        gsl_vector_memcpy(netvec,Ivec_p);                                       // copy biomass gain by ingestion into netvec
        gsl_vector_sub(netvec,Dvec_p);                                          // substract biomass lost by death
        gsl_vector_sub(netvec,Xvec_p);                                          // substract biomass lost by respiration
        gsl_vector_div(netvec,Bvec_p);                                          // divide by biomass to get net growth rate

        for (i=0;i<S;i++)
        {
            if (gsl_vector_get(netvec, i)!= gsl_vector_get(netvec, i)) // Replace potential NaNs with 0
                gsl_vector_set(netvec, i, 0);
        }

//  ************* dispersal rates *************
        for(i=0; i<S_b; i++)		                                            // basal (plant) species -- passive dispersal (random) TEST DIFFERENT FUNCTIONS FOR PLANT DISPERSAL (include different dispersal directions??)
            gsl_vector_set(Emvec_p,i,(emigr_a_S_b/(1+(exp(emigr_b_S_b*(-gsl_vector_get(Rvec, i) - gsl_vector_get(netvec, i)))))));

        for(i=S_b; i<S; i++)		                                            // consumer (animal) species -- active dispersal (body mass dependent)
            gsl_vector_set(Emvec_p,i,(emigr_a_S_c/(1+(exp(-emigr_b_S_c*(-gsl_vector_get(Rvec, i) - gsl_vector_get(netvec, i)))))));

//  ************* cp net growth rates *************
        gsl_vector_memcpy(netgrowth_p,netvec);

        // Final emigration rate
        gsl_vector_mul(Emvec_p,Bvec_p);

    } // end loop patches


//  ********** free memory ************
    gsl_vector_free(Prey);
    gsl_vector_free(svec);
    gsl_vector_free(rvec);
    gsl_vector_free(nvec);
    gsl_vector_free(mvec);
    gsl_vector_free(netvec);
    gsl_vector_free(pvec);
    gsl_matrix_free(mmat);

    return;
}


/* The following contains functions to generate the spatial topology structures.
 * The network structures are returned in the form of weighted success matrices SW (S*Z,Z),
 * i.e., SW is non-zero if and only if patch i is linked with patch j, and the
 * actual value equals the success weighted dispersal rate of this link.
 *
 * Spatial topology generating functions:
 *  - Generate_RGG_structure - random geometric graph (RGG) build
 *
 * Requirements for spatial topology generating functions:
 *  1. pointer to an instance of the random number generator; gsl_rng *r
 *  2. double array of length S for the body masses of the species; gsl_vector *mass
 *  3. pointer to gsl_vector of length 2*Z to be filled  with the location of the patches; gsl_vector *Loc (returned)
 *  4. pointer to gsl_matrix of size S*Z,Z to be filled with success weighted dispersal rates; gsl_vector *SW (returned)
 *
 * The following also contains functions which are used within the above functions.
*/

//***************************************************************************************
// Basic RGG network structure
//***************************************************************************************
static void Generate_RGG_structure(gsl_vector *mass, gsl_vector *D_Max, gsl_vector *Loc, gsl_matrix *SW,  gsl_vector *con, double params[], double RGG_params[])
{
	if(INVASION)
		getInputInvasionDisp(D_Max); 
	else if (!INVASION) 
		getInputDisp(D_Max); 
	
	int i,j,k,m,n; 	

    double X_1,X_2,Y_1,Y_2,DistX,DistY,d_mn,d_spec,sum,temp5,D_max_scale,D_max_max;

	D_max_max=0.5;                                                              // maximal dispersal distance for the largest species
    D_max_scale=D_max_max/pow(1e12,eps);                                        // scaling factor acording to D_max_max for D_max for maximal body size of 1e12 for consumers

    gsl_matrix *D_spez = gsl_matrix_calloc(D,Z);                                // matrix containing the specific distances d_spec=d_mn / d_max
    gsl_matrix *SU = gsl_matrix_calloc(D,Z);                                    // dispersal success
    gsl_matrix *W = gsl_matrix_calloc(D,Z);                                     // dispersal weights
    gsl_matrix *DP = gsl_matrix_calloc(Z,Z);                                    // distance between patches (P_Dist)
    
    for (i=0; i<S; i++)
        *(params+2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+8+D+DN+i) = gsl_vector_get(D_Max,i);  // save D_Max in params 
	
	getInputLocX(Loc); 	 														// import Z X locations 
	getInputLocY(Loc); 															// import Z Y locations
														
    calc_patch_distances(DP, Loc, RGG_params); 	                                // calculate distances d_mn betweeen patches with pythagoras

    for (m=0; m<Z; m++)
    {
        for (n=m+1; n<Z; n++)
        {
            // species-specific distances => can species i connect two patches
            for (i=0; i<S; i++)
            {
                // Create submatrices
                gsl_matrix_view D_spez_i_mat = gsl_matrix_submatrix(D_spez,i*Z,0,Z,Z);
                gsl_matrix *D_spez_i = &D_spez_i_mat.matrix;
                gsl_matrix_view SU_i_mat = gsl_matrix_submatrix(SU,i*Z,0,Z,Z);
                gsl_matrix *SU_i = &SU_i_mat.matrix;
                gsl_matrix_view W_i_mat = gsl_matrix_submatrix(W,i*Z,0,Z,Z);
                gsl_matrix *W_i = &W_i_mat.matrix;

                d_spec = gsl_matrix_get(DP, m,n)/gsl_vector_get(D_Max,i);       // specific distance d_spec = d_mn/d_max

                if(d_spec <=1)                                                  // d_max_i >= distance between patch m and n d_mn
                {
                    gsl_matrix_set(D_spez_i, m, n, d_spec);
                    gsl_matrix_set(D_spez_i, n, m, d_spec);

                    // calculate dispersal success s_i,nm = (1-d_i,nm)^theta and save it in matrix SU_S*Z,Z
                    gsl_matrix_set(SU_i, m, n, pow((1-d_spec),theta));          // just setting the success to 1 makes solver get stuck early in simulation
                    gsl_matrix_set(SU_i, n, m, pow((1-d_spec),theta));
                }                                                               // else 0 (initial value)
            }
        }
    }

    for(i=0; i<S; i++)
    {
        gsl_matrix_view D_spez_i_mat = gsl_matrix_submatrix(D_spez,i*Z,0,Z,Z);
        gsl_matrix *D_spez_i = &D_spez_i_mat.matrix;
        gsl_matrix_view W_i_mat = gsl_matrix_submatrix(W,i*Z,0,Z,Z);
        gsl_matrix *W_i = &W_i_mat.matrix;
        gsl_matrix_view SU_i_mat = gsl_matrix_submatrix(SU,i*Z,0,Z,Z);
        gsl_matrix *SU_i = &SU_i_mat.matrix;
        gsl_matrix_view SW_i_mat = gsl_matrix_submatrix(SW,i*Z,0,Z,Z);
        gsl_matrix *SW_i = &SW_i_mat.matrix;

        for(j=0; j<Z; j++)                                                      // iteration over columns
        {
            sum=0;                                                              // counts from 0 for every row
            for(k=0; k<Z; k++)                                                  // iteration over number of patch es
            {
                if(gsl_matrix_get(D_spez_i, j, k)!=0)
                {
                    sum += 1-gsl_matrix_get(D_spez_i, j, k);
                    gsl_matrix_set(W_i,j,k, (1-gsl_matrix_get(D_spez_i, j, k)));
                }
            }
            gsl_vector_view vrow = gsl_matrix_row(W_i, j);
            gsl_vector *vrow_j = &vrow.vector;

            if(sum!=0)
                gsl_vector_scale(vrow_j, 1/sum);
        
         }
        // calculate success-weighted dispersal matrix elementwise SU*W
        // SW is the matrix we pass to the function (it is our weighted-success-matrix (the end results we want to attain))
        gsl_matrix_memcpy(SW_i, W_i);                                           // copy contents of W_i to SW_i
        gsl_matrix_mul_elements(SW_i, SU_i);                                    // multiply all elements of W (in SW) and SU to get SW (end result)
        
	}
    

    calc_spatial_connectance(SW,con,RGG_params);                                // calculate the connectance within the RGG

    //*** free memory *****
    gsl_matrix_free(D_spez);
    gsl_matrix_free(SU);
    gsl_matrix_free(DP);
    gsl_matrix_free(W);
}

//***************************************************************************************
// Calculate the distances between patches, saved in gsl_matrix P_Dist
//***************************************************************************************
/*
A Function which returns a distance matrix (P_Dist) Z*Z, reflecting the distances between each patch.
Each location of the patch is saved in the vector Loc (Z*2).
RGG_params is an empty vector which is filled with some properties of the RGG such as Mean patch distance & mean nearest neighbor distance

1. Get the locations of each patch (X,Y)
2. Calculate distances according to Pythagoras
3. Calcualte mean distance and nearest neighbor distance and save them in RGG_params
4. Free memory
*/
static void calc_patch_distances(gsl_matrix *P_Dist, gsl_vector *Loc, double RGG_params[])
{

    int m,n;
    double X_1,X_2,Y_1,Y_2,DistX,DistY, d_mn, sum_nn_dist;

    gsl_vector_view X_Vec = gsl_vector_subvector_with_stride(Loc,0,2,Z);        // stride = 2 (alternating XY values)
    gsl_vector *XVec=&X_Vec.vector;
    gsl_vector_view Y_Vec = gsl_vector_subvector_with_stride(Loc,1,2,Z);
    gsl_vector *YVec=&Y_Vec.vector;

    // calculate distances d_mn betweeen patches with pythagoras
    for (m=0; m<Z; m++)
    {

        X_1 = gsl_vector_get(XVec,m);                                           // take one X-value for calculations
        Y_1 = gsl_vector_get(YVec,m);                                           // take one Y-value for calculations

        for (n=m+1; n<Z; n++)
        {
            X_2 = gsl_vector_get(XVec,n);                                       // choose another X-Value to be substracted by X1
            DistX = X_1-X_2;                                                    // Substraction

            Y_2 = gsl_vector_get(YVec,n);                                       // choose another X-Value to be substracted by X1
            DistY = Y_1-Y_2;                                                    // Substraction

            // Pythagoras for distances between two patches
            d_mn = sqrt(pow(DistY,2)+pow(DistX,2));
            gsl_matrix_set(P_Dist, m, n, d_mn);
            gsl_matrix_set(P_Dist, n, m, d_mn);
        }
    }

    // calcualte mean distance and nearest neighbor distance
    sum_nn_dist = 0;
    gsl_vector *temp = gsl_vector_calloc(Z);
    gsl_vector *temp2 = gsl_vector_calloc(Z*Z);

    for(int i=0; i<Z; i++)
    {

        gsl_vector_view P_Dist_Row = gsl_matrix_row(P_Dist,i);                  // create matrix view for each row i (length Z)
        gsl_vector *PDist_Row=&P_Dist_Row.vector;

        for(int j=0; j<Z; j++)
            gsl_vector_set(temp2,i*Z+j,gsl_vector_get(PDist_Row,j));            // save all distances in temp 2

        gsl_vector_set(PDist_Row,i,5);                                          // set distance to patch itself to X to remove 0 as min dist
        sum_nn_dist += gsl_vector_min(PDist_Row);                               // sum up of minimum distances for each patch
        gsl_vector_set(temp,i,gsl_vector_min(PDist_Row));
        gsl_vector_set(PDist_Row,i,0);                                          // set distance to patch itself back to 0

    }

    double mean_dist = calc_mean(temp2,0,Z*Z,(Z*Z)-Z);
    double sd_nn_dist = calc_sd(temp,0,Z,Z);
    double sd_dist = calc_sd(temp2,0,Z*Z,(Z*Z)-Z);

    RGG_params[0] = mean_dist;                                                  // mean distance
    RGG_params[8]= sd_dist;                                                     // standard deviation of distance

    RGG_params[1]= sum_nn_dist/Z;                                               // mean nearest neighbor distance
    RGG_params[9]= sd_nn_dist;                                                  // standard deviation nearest neighbor distance

    // *** free memory ***
    gsl_vector_free(temp);
    gsl_vector_free(temp2);
}

//***************************************************************************************
// Calculate spatial connectance, save in RGG_parms[]
//***************************************************************************************
/*
A Function which returns the spatial connectance of the patches. The input needed consits of the success weighted matrix (SW)
which provides information of all realized links between patches. The realized conectance for each species is saved in the vector con
and for evaluation RGG quantities are saved in RGG_params such as mean spatial connactence

1. Realized links of each species
2. Spatial connectance of each species
3. Overall mean spatial connectance and SD are calculated
*/
static void calc_spatial_connectance(gsl_matrix *SW, gsl_vector *con, double RGG_params[])
{

    int i, j, k;
    double count2, mean_all, mean_b, mean_c, sd_all, sd_b, sd_c;

    for(i=0; i<S; i++)
    {
        gsl_matrix_view SW_i_mat = gsl_matrix_submatrix(SW,i*Z,0,Z,Z);          // submatrix of success weighted matrix
        gsl_matrix *SW_i = &SW_i_mat.matrix;

        count2=0;
        for(j=0; j<Z; j++)
        {
            for(k=0; k<Z; k++)
            {
                if(gsl_matrix_get(SW_i,j,k)>0) count2 ++;                       // number of realized links
            }
        }

        gsl_vector_set(con,i,count2/(pow(Z,2)-Z));                              // calculate spatial connectance on species-level con = no.realized.links/no.possible.links (Z^2-Z)
    }

    mean_all = calc_mean(con,0,S,S);
    mean_b = calc_mean(con,0,S_b,S_b);
    mean_c = calc_mean(con,S_b,S,S_c);

    sd_all = calc_sd(con,0,S,S);
    sd_b = calc_sd(con,0,S_b,S_b);
    sd_c = calc_sd(con,S_b,S,S_c);

    RGG_params[2] = mean_all;                                                   // mean spatial connectance
    RGG_params[3] = sd_all;                                                     // sd spatial connectance

    RGG_params[4] = mean_b;                                                     // mean spatial connectance plants
    RGG_params[5] = sd_b;                                                       // sd spatial connectance plants

    RGG_params[6] = mean_c;                                                     // mean spatial connectance consumers
    RGG_params[7] = sd_c;                                                       // sd spatial connectance consumers

}

//****************************************************************************************************************************
// Write output to files: web.csv; mass.csv; global.csv; params.csv (optional).
//****************************************************************************************************************************
/* Output files containing all relevant simulation parameters and variables are created to evaluate the simulation.
1. file 1 (web_*.out) contains the adjacency matrix with interaction strengths (S*S).
2. file 2 (mass_*.out) contains relevant patch- and species-specific quantities for posthoc analysis such as mean biomasses, body masses, variability or mean net growth rate.
3. file 3 (global_*.out) contains important global (i.e. simulation-specific) parameters including patch number, mean patch distance, maximal emigration rates or the shape of the emigration function.
4. Option to call function to write all simulation parameters stored in params[] to an output file (params_*.out).
*/
static void output(gsl_matrix *SW, gsl_vector *mass,  gsl_vector *con, double iniB[], double meanB[], double meanB_tot[], double CV[], double CV_tot[], double params[], int paramslength, double Biomass_Temp[], double RGG_params[], double net_growth_Temp[], gsl_vector *InvSppvec, double B[])
{
    int i,j,k,m;

	if(INVASION) 
		getInputInvasionInvSpp(InvSppvec);		
	else 
		gsl_vector_set_zero(InvSppvec); 
					
    FILE *file1,*file2,*file3;

    char str[99];
    char web_out[99] = "web_";
    char mass_out[99] = "masses_";
    char global_out[99] = "global_";
    
    // Directory of the files
    char resdir[99] = "res_";
	strcpy(resdir,output_dir); 
	
    strcpy(web_out,resdir);
    strcpy(mass_out,resdir);
    strcpy(global_out,resdir);

    // Naming of the files
    strcat(web_out,"web_");
    strcat(mass_out,"mass_");
    strcat(global_out,"global_");
        
    strcat(web_out,web);
    strcat(mass_out,web);
    strcat(global_out,web);
    
    strcat(web_out,"_");
    strcat(mass_out,"_");
    strcat(global_out,"_");
   
    strcat(web_out,landscape);
    strcat(mass_out,landscape);
    strcat(global_out,landscape);
    
    strcat(web_out,"_");
    strcat(mass_out,"_");
    strcat(global_out,"_");

    strcat(web_out,bash_remspp);
    strcat(mass_out,bash_remspp);
    strcat(global_out,bash_remspp);
    
    strcat(web_out,"_");
    strcat(mass_out,"_");
    strcat(global_out,"_");
    
    strcat(web_out,bash_nutsupply);
    strcat(mass_out,bash_nutsupply);
    strcat(global_out,bash_nutsupply);
    
    strcat(web_out,"_");
    strcat(mass_out,"_");
    strcat(global_out,"_");
    
    strcat(web_out,bash_invasion);
    strcat(mass_out,bash_invasion);
    strcat(global_out,bash_invasion);
    
    // Data type of the files
    strcat(web_out,".csv");
    strcat(mass_out,".csv");
    strcat(global_out,".csv");

    file1 = fopen(web_out,"w");
    file2 = fopen(mass_out,"w");
    file3 = fopen(global_out,"w");

	//web.csv (SxS matrix) 
    for(i=0; i<S; i++)
    {
        for(j=0; j<S; j++)
            fprintf(file1,"%g ",params[i*S+j]);
        fprintf(file1,"\n");
    }

    double Surv=0;                                                              // global survival
    double con_count=0;                                                         // number of extant consumer species
    double bas_count=0;                                                         // number of extant basal species

    for(int i=0; i<S; i++)
    {
        double B_Species = 0;

        for (int m=0; m<Z; m++)
        {
            B_Species = B_Species + B[i+S*m];
        }
        if(B_Species > GLOBAL_EXTINCT)
            Surv++;
        if(i >= S_b && B_Species > GLOBAL_EXTINCT)
            con_count++;
        if(i < S_b && B_Species > GLOBAL_EXTINCT)
            bas_count++;
    }

    gsl_vector_view P_vec=gsl_vector_view_array(params+2*S*S+2*S,S);            // basal vector (1 for basal species, 0 otherwise)
    gsl_vector *Pvec=&P_vec.vector;
	
    // global.csv (one line with important global simulation variables)
    fprintf(file3,"%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","web","landscape","number.of.spp","number.of.plants","number.of.consumers",
		"rem.spp", "number.of.patch","rng.seed","max.emigr.rate.plants","shape.emigr.function.plants","max.emigr.rate.consumers","shape.emigr.function.consumers","D_0","theta","eps",
		"mean.patch.dist","sd.patch.dist","mean.nn.dist","sd.nn.dist","mean.con.rgg","sd.con.rgg","mean.con.rgg.plants","sd.con.rgg.plants","mean.con.rgg.consumers","sd.con.rgg.consumers",
		"nutrient.supply","extant.species","extant.plants","extant.consumers","INVASION","JOINT","ISOLATED","number.of.successfull.migrations","number.of.totalmigrations","tend_stopped");
    fprintf(file3,"%d,%d,%d,%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%.9g,%.9g,%.9g\n",
			webi,landscapei,S,S_b,S_c,remspp,Z,seed,emigr_a_S_b,emigr_b_S_b,emigr_a_S_c,emigr_b_S_c,D_0,theta,eps,RGG_params[0],RGG_params[8],RGG_params[1],RGG_params[9],RGG_params[2],
			RGG_params[3],RGG_params[4],RGG_params[5],RGG_params[6],RGG_params[7],mu_S, Surv, bas_count, con_count, INVASION, JOINT, ISOLATED, counting_migration2, counting_migration3,t_output);
   
    // mass.csv (table of length Z*(S+N) with patch- and species-specific biomass densities resp. nutrient concentrations)
    fprintf(file2, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n","patch","species","bodymass", "if.basal.spp", "if.inv.spp","init.biomass", "mean.biomass","biomass.variance",
		"tot.mean.biomass","tot.biomass.variance","Biomassess_tend-20.000","Biomassess_tend-10.000","Biomassess_tend","spatial.connectance","Mean.net.growth1","Mean.net.growth2",
		"net.growth.sd"); 
    
    for(j=0; j<Z; j++)
    {
        for(i=0; i<S; i++){
			int ii; 
			if(!INVASION){
				if(i<remspp)
					ii = i;
			else 
				ii = i+1; 
			}
			if(INVASION) ii = i; 
		fprintf(file2,"%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g\n", j, ii, gsl_vector_get(mass,i),gsl_vector_get(Pvec,i),gsl_vector_get(InvSppvec,i),iniB[j*S+i],meanB[j*S+i],CV[j*S+i],meanB_tot[i],CV_tot[i],Biomass_Temp[i+j*S],Biomass_Temp[D+DN+i+j*S],Biomass_Temp[2*(D+DN)+i+j*S],gsl_vector_get(con,i),net_growth_Temp[j*S+i],net_growth_Temp[D+j*S+i],net_growth_Temp[2*D+j*S+i]);
		}

	}
	
	
	for(j=0; j<Z; j++) 
	{
        //nurients
        for(m=(D+j*N); m<(D+j*N+1); m++)
            fprintf(file2,"%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%d\n",j,1000,0,2,2,iniB[m],meanB[m],CV[m],meanB_tot[m],CV_tot[m],Biomass_Temp[m],Biomass_Temp[D+DN+m],Biomass_Temp[2*(D+DN)+m],0,0,0,0);
       
        for(m=(D+j*N+1); m<(D+j*N+N); m++)
            fprintf(file2,"%d,%d,%d,%d,%d,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%.9g,%d,%d,%d,%d\n",j,2000,0,2,2,iniB[m],meanB[m],CV[m],meanB_tot[m],CV_tot[m],Biomass_Temp[m],Biomass_Temp[D+DN+m],Biomass_Temp[2*(D+DN)+m],0,0,0,0);
    }

    Write_params_to_file(params,paramslength); // write all simulation parameters from the 'params' array to a file.
    
	// close output files
    fclose(file1);
    fclose(file2);
    fclose(file3);

    return;
}

//***************************************************
//  Prepare Timeseries for a file
//***************************************************
/*
Preparation of file to save the timeseries with biomasses for each species on each patch at each evaluated time step.
*/
static void Prepare_timeseries_file(FILE *timeseries)
{
	
	strcat(timeseries_file,web); 
	strcat(timeseries_file,"_"); 
	strcat(timeseries_file,landscape); 
	strcat(timeseries_file,"_"); 
	strcat(timeseries_file,bash_remspp); 
	strcat(timeseries_file,"_"); 
	strcat(timeseries_file,bash_nutsupply); 
	strcat(timeseries_file,"_"); 
	strcat(timeseries_file,bash_invasion);
	strcat(timeseries_file,".csv"); 

    timeseries = fopen(timeseries_file,"w");
    fprintf(timeseries,"%s,","time");

    for(int i=0; i<D+DN; i++)
        fprintf(timeseries,"%d,",i);
   
    fprintf(timeseries,"\n");
    fclose(timeseries);

    return;
}

//***************************************************
// Write Timeseries to a file
//***************************************************
/*
Function to write down the biomasses for each time step and appends the values in a new row. If a timeseries file is created, mean biomasses are not calculated.
Input: biomass vector (B); current time step (t); file to store the timeseries.
*/
static void Write_timeseries_to_file(double B[], double t, FILE *timeseries)
{
 
    timeseries = fopen(timeseries_file,"a");
    fprintf(timeseries,"%g,",t);                                               // Print time t.
    
    for(int i=0; i<D+DN; i++)
        fprintf(timeseries,"%.8g,", B[i]);                                     // Print for each species and patch biomass densities resp. nutrient concentrations.
    fprintf(timeseries,"\n");

    fclose(timeseries);

    return;
}

//*****************************************************************
//  Write all fixed parameter values ('params' array to a file
//*****************************************************************
/*
Function to save all parameters stored in the 'params' array including handling times, respiration rates, assimilation efficiencies, ...
*/
static void Write_params_to_file(double params[],int paramslength)
{
    FILE *file1;
    char str[99];
    char params_out[99] = "params_";
    strcpy(params_out,output_dir);
    strcat(params_out,"params_");
    strcat(params_out,web);
	strcat(params_out,"_");
	strcat(params_out,landscape);
    strcat(params_out,"_");
    strcat(params_out,bash_remspp);
    strcat(params_out,"_");
    strcat(params_out,bash_nutsupply); 
    strcat(params_out,"_");
	strcat(params_out,bash_invasion);     
    strcat(params_out,".csv");      

    file1 = fopen(params_out,"w");

    fprintf(file1, "%s, %s\n", "variable","value");
    for (int i=0; i<S*S; i++)
        fprintf(file1, "%s, %g\n", "attack.rates",params[i]);
    for (int i=S*S; i<2*S*S; i++)
        fprintf(file1, "%s, %g\n", "handling.times",params[i]);
    for (int i=2*S*S; i<2*S*S+S_b; i++)
        fprintf(file1, "%s, %g\n", "plant.respiration",params[i]);
    for (int i=2*S*S+S_b; i<2*S*S+S; i++)
        fprintf(file1, "%s, %g\n", "animal.respiration",params[i]);
    for (int i=2*S*S+S; i<2*S*S+2*S; i++)
        fprintf(file1, "%s, %g\n", "interference",params[i]);
    for (int i=2*S*S+2*S; i<2*S*S+3*S; i++)
        fprintf(file1, "%s, %g\n", "plant.identifier",params[i]);
    for (int i=2*S*S+3*S; i<2*S*S+4*S; i++)
        fprintf(file1, "%s, %g\n", "hill.exponent",params[i]);
    for (int i=2*S*S+4*S; i<2*S*S+5*S; i++)
        fprintf(file1, "%s, %g\n", "bodymasses",params[i]);
    for (int i=2*S*S+5*S; i<2*S*S+5*S+S_b; i++)
        fprintf(file1, "%s, %g\n", "plant.assimilarion.efficiency",params[i]);
    for (int i=2*S*S+5*S+S_b; i<2*S*S+6*S; i++)
        fprintf(file1, "%s, %g\n", "animal.assimilarion.efficiency",params[i]);
    for (int i=2*S*S+6*S; i<2*S*S+6*S+N*S_b; i++)
        fprintf(file1, "%s, %g\n", "nutrient.uptake.half.sat.densities",params[i]);
    for (int i=2*S*S+6*S+N*S_b; i<2*S*S+6*S+N*S_b+N; i++)
        fprintf(file1, "%s, %g\n", "nutrient.supply.concentrations",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N; i<2*S*S+6*S+N*S_b+N+S_b; i++)
        fprintf(file1, "%s, %g\n", "max.nutrient.uptake.rates",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N+S_b; i<2*S*S+6*S+N*S_b+N+S_b+1; i++)
        fprintf(file1, "%s, %g\n", "content.of.first.nutrient.in.plant",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N+S_b+1; i<2*S*S+6*S+N*S_b+N+S_b+2; i++)
        fprintf(file1, "%s, %g\n", "content.of.second.nutrient.in.plant",params[i]);
    for (int i=2*S*S+6*S+N*S_b+N+S_b+2; i<2*S*S+6*S+N*S_b+N+S_b+2+S; i++)
        fprintf(file1, "%s, %g\n", "extinction.threshold",params[i]);
    for (int i=2*S*S+7*S+N*S_b+N+S_b+2; i<2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z; i++)
        fprintf(file1, "%s, %g\n", "SW.dispersal.matrix",params[i]);
    for (int i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z; i<2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z; i++)
        fprintf(file1, "%s, %g\n", "patch.location.XY.alternating",params[i]);
    int i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z;
    fprintf(file1, "%s, %g\n", "max.animal.emigration.rate",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+1;
    fprintf(file1, "%s, %g\n", "shape.of.animal.emigration.function",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+2;
    fprintf(file1, "%s, %g\n", "max.plant.emigration.rate",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+3;
    fprintf(file1, "%s, %g\n", "shape.of.plant.emigration.function",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+4;
    fprintf(file1, "%s, %g\n", "RGG.epsilon",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+5;
    fprintf(file1, "%s, %g\n", "RGG.dispersal.distance.D_0",params[i]);
    i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+6;
    fprintf(file1, "%s, %g\n", "RGG.dispersal.success.theta",params[i]);	
	for (int i=2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+7; i<2*S*S+7*S+N*S_b+N+S_b+2+S*Z*Z+2*Z+7+D; i++)
        fprintf(file1, "%s, %g\n", "net.growth.rate",params[i]);   
    for (int i=2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+8+D+DN; i<2*S*S+7*S+S_b*N+N+S_b+2+S*Z*Z+2*Z+8+D+DN+S; i++)
        fprintf(file1, "%s, %g\n", "D_Max",params[i]); 
  
	fclose(file1);
	
    return;
}

//***************************************************************************
// Read in input file to get the number of species in the web
//***************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: two columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0).
The function defines the global variables of species number (S), number of basal (S_b) and consumer species (S_c).
*/
static void getInputSpp()
{  
    char spp_in[99] = "bodymass_";
	strcpy(spp_in,directory);
    strcat(spp_in,"webs/bodymass_");
	strcat(spp_in,web); 
	strcat(spp_in,"_"); 
    strcat(spp_in,bash_remspp);
    strcat(spp_in,".csv");
    
    std::ifstream input(spp_in);                                       
    std::vector<double> spp;
    std::string line;

    if(!input.is_open())
    {
        exit(1);
    }

    double col1, col2; 

    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line, ',');  
        spp.push_back(col1);
    }

    S_b=0;
    S_c=0;
    
    for(int i=1; i<spp.size(); i++)										// count number of basal and consumber species.
    {
        if(spp[i]>0)
            S_b += 1;
        else
            S_c += 1;
    }
    
    input.close();  // close the input file.
    
    return;
}

//  ***************************************************************************
//  Read in input file containg species body masses to create web.
//  ***************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: two columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0).
Body masses are saved in a vector (mass).
*/
static void getInputMass(gsl_vector *mass)
{
    char mass_in[99] = "bodymass_";
	strcpy(mass_in,directory);
    strcat(mass_in,"webs/bodymass_");
	strcat(mass_in,web); 
	strcat(mass_in,"_"); 
    strcat(mass_in,bash_remspp);
    strcat(mass_in,".csv");
    
    std::ifstream input(mass_in);                                       
    std::vector<double> masses;
    std::vector<double> D_Maxes;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "getInputMass: There was a problem opening the bodymass file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line, ','))
    {
        std::istringstream iss(line); 
        iss >> col1 >> col2;
		getline(input, line);  
        masses.push_back(col1);											// extract body masses from input file.
    }

    for(int i=0; i<(masses.size()-1); i++)
        gsl_vector_set(mass,i,masses[i+1]);								// save body mass of species i in vector mass.
 
    input.close();  // close the input file
    
    return;
}

//  ***************************************************************************
//  Read in input file containg species body masses to create web.
//  ***************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: two columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0).
Dispesral distances are saved in a vector (D_Max).
*/

static void getInputDisp(gsl_vector *D_Max)
{
    char mass_in[99] = "bodymass_";
	strcpy(mass_in,directory);
    strcat(mass_in,"webs/bodymass_");
	strcat(mass_in,web); 
	strcat(mass_in,"_"); 
    strcat(mass_in,bash_remspp);
    strcat(mass_in,".csv");
    
    std::ifstream input(mass_in);                                       
    std::vector<double> D_Maxes;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "getInputDisp: There was a problem opening the dispersal file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line,',');
        getline(input, line,',');
        D_Maxes.push_back(col1);										    // extract body masses from simulation run without invasive species from input file.
    }
    for(int i=0; i<(D_Maxes.size()-1); i++)
        gsl_vector_set(D_Max,i,D_Maxes[i+1]);								// save dippersal distances of species i in D_Max.
    
    input.close();  // close the input file
    
    return;
}

//  ***************************************************************************
//  Read in input file containg species body masses to create web.
//  ***************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: three columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0); column 3: indicates if invasive species (1) or not (0). 
Body masses are saved in a vector (mass).
*/
static void getInputInvasionBodymass(gsl_vector *mass)
{

    char inv_in[99] = "inputBodymass_";
	strcpy(inv_in,invasion_dir);
    strcat(inv_in,"inputBodymass_");
	strcat(inv_in,web); 
	strcat(inv_in,"_"); 
    strcat(inv_in,landscape);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_remspp);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_nutsupply);
    strcat(inv_in,"_"); 
    strcat(inv_in,"0");														// import Biommasses_tend from simulation run without invasive species
    strcat(inv_in,".csv");
    
    std::ifstream input(inv_in);                                       
    std::vector<double> invMass;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "getInputInvasionBodymass: There was a problem opening the input file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line, ','))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line);
        invMass.push_back(col1);										// extract body masses from simulation run without invasive species from input file.
    }

	
	for(int i=0; i<(invMass.size()-1); i++)
   	  gsl_vector_set(mass,i,invMass[i+1]);							// save body mass of species i in vector mass.
  
	input.close();  // close the input file.

    return;
}


//  ***************************************************************************
//  Read in input file containg species body masses to create web.
//  ***************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: three columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0); column 3: indicates if invasive species (1) or not (0). 
Dipsersal distances are saved in a vector (D_Max).
*/
static void getInputInvasionDisp(gsl_vector *D_Max)
{

    char inv_in[99] = "inputBodymass_";
	strcpy(inv_in,invasion_dir);
    strcat(inv_in,"inputBodymass_");
	strcat(inv_in,web); 
	strcat(inv_in,"_"); 
    strcat(inv_in,landscape);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_remspp);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_nutsupply);
    strcat(inv_in,"_"); 
    strcat(inv_in,"0");														// import Biommasses_tend from simulation run without invasive species
    strcat(inv_in,".csv");
    
    std::ifstream input(inv_in);                                       
    std::vector<double> invD_Max;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "getInputInvasionBodymass: There was a problem opening the input file\n";
        exit(1);
    }

    double col1,col2;

      while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line,',');
        getline(input, line,',');
        getline(input, line,',');
        getline(input, line,',');
        invD_Max.push_back(col1);										// extract body masses from simulation run without invasive species from input file.
    }

	
    for(int i=0; i<(invD_Max.size()-1); i++)
   	  gsl_vector_set(D_Max,i,invD_Max[i+1]);							// save dispersal distances of species i in vector D_Max.
      
	input.close();  // close the input file.
    return;
}


//  *****************************************************************************************
//  Read in input file invasion/inputBodymass_*.csv to get which species is invasive spp
//  *****************************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: three columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0); column 3: indicates if invasive species (1) or not (0). 
Indicator for invasive spp is saved in a vector (InvSppvec).
*/
static void getInputInvasionInvSpp(gsl_vector *InvSppvec)
{

    char inv_in[99] = "inputBodymass_";
	strcpy(inv_in,invasion_dir); 
    strcat(inv_in,"inputBodymass_");
	strcat(inv_in,web); 
	strcat(inv_in,"_"); 
    strcat(inv_in,landscape);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_remspp);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_nutsupply);
    strcat(inv_in,"_"); 
    strcat(inv_in,"0");														
    strcat(inv_in,".csv");
    
    std::ifstream input(inv_in);                                       
    std::vector<double> invInvSpp;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "getInputInvasionInvSpp: There was a problem opening the input file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line,',');
        getline(input, line,',');
        invInvSpp.push_back(col1);										// extract body masses from simulation run without invasive species from input file.
    }

	for(int i=0; i<(invInvSpp.size()-1); i++)
   	  gsl_vector_set(InvSppvec,i,invInvSpp[i+1]);						// save whether species i is invasive species in vector InvSppvec.
   
    input.close();  // close the input file.

    return;
}

//  ***********************************************************************************
//  Read in input file containg biomasses from simulation without invasive species. 
//  ***********************************************************************************
/* 
Import bioamasses from tend from the simulation without invasive species to initialize biomassess after invasion.
Structure of input files: one colunn. Column 1: Biomass_tend.
Biomasses are saved in an array (B[]). 
*/
static void getInputInvasionBiomass(double B[])
{
    char inv_in[99] = "inputBiomass_";
	strcpy(inv_in,invasion_dir);
    strcat(inv_in,"inputBiomass_");
	strcat(inv_in,web); 
	strcat(inv_in,"_"); 
    strcat(inv_in,landscape);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_remspp);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_nutsupply);
    strcat(inv_in,"_"); 
    strcat(inv_in,"0");														// import Biommasses_tend from simulation run without invasive species
    strcat(inv_in,".csv");
    
    std::ifstream input(inv_in);                                       
    std::vector<double> invB;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "getInputInvasionBiomass: There was a problem opening the input file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        invB.push_back(col1);											// extract Biomassess_tend from simulation run without invasive species from input file.
    }

	for(int i=0; i<(invB.size()-1); i++)
	{
		B[i] = invB[i+1];											// save Biomassess_tend of species i in vector iniB
	}
	
    input.close();  // close the input file.
    
    return;
}

//  ***************************************************************************
//  Read in input file containg species body masses to create web.
//  ***************************************************************************
/*
Option to use a specific foodweb with a given number of species including their body masses.
Structure of input files: three columns. Column 1: body mass; column 2: inidicates if basal species (1) or consumer (0); column 3: indicates if invasive species (1) or not (0). 
Count of consumer and basal species. 
*/
static void getInputInvasionSpp()
{
	
    char inv_in[99] = "inputBodymass_";
	strcpy(inv_in,invasion_dir);
    strcat(inv_in,"inputBodymass_");
	strcat(inv_in,web); 
	strcat(inv_in,"_"); 
    strcat(inv_in,landscape);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_remspp);
  	strcat(inv_in,"_"); 
    strcat(inv_in,bash_nutsupply);
    strcat(inv_in,"_"); 
    strcat(inv_in,"0");														// import Biommasses_tend from simulation run without invasive species
    strcat(inv_in,".csv");
        
    std::ifstream input(inv_in);                                       
    std::vector<double> invSpp;
    std::string line;

    if(!input.is_open())
    {
        exit(1);
    }

    double col1,col2;
    
    while(std::getline(input, line))
    {
	    std::istringstream iss(line);
        iss >> col1 >> col2;
        std::getline(input, line,',');  
        invSpp.push_back(col1);
    }

    S_b=0;
    S_c=0;
    
    for(int i=1; i<invSpp.size(); i++)									// count number of basal and consumber species.
    {
        if(invSpp[i]>0)
            S_b += 1;
        else
            S_c += 1;
		
    }
        
    input.close();  // close the input file
    
    return;
}


//***************************************************************************
// Read in landscape file to get the number of patches Z in the landscape
// ***************************************************************************
/*
Option to use a specific landscape with Z patches.
Structure of input files: two columns. Column 1: x-coordinates of patches; column 2: y-coordinates of patches.
*/
static void getInputZ()
{
    char landsc_in[99] = "landscape_";
	strcpy(landsc_in,directory);
    strcat(landsc_in,"landscapes/landscape_");
    strcat(landsc_in,landscape);
    strcat(landsc_in,".csv");
      
    std::ifstream input(landsc_in);                               
    std::vector<double> x;
    std::vector<double> y;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the landscape file\n";
        exit(1);
    }

    double col1, col2; 

    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line, ',');  
        x.push_back(col1);
     }
    
    Z=0;

    for(int i=1; i<x.size(); i++)										// count number of patches.
      Z += 1; 
     	   
    input.close();  // close the input file

    return;
}

//***************************************************
// Read in X locations
//***************************************************
/*
Option to use a specific landscape with Z patches.
Structure of input files: two columns. Column 1: X locations of patches; column 2: Y locations of patches.
*/
static void getInputLocX(gsl_vector *Loc)
{
    
    char landsc_in[99] = "landscape_";
	strcpy(landsc_in,directory);
    strcat(landsc_in,"landscapes/landscape_");
    strcat(landsc_in,landscape);
    strcat(landsc_in,".csv");
    

    std::ifstream input(landsc_in);                                       
    std::vector<double> x;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the landscape file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line, ','))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line);  
        x.push_back(col1);												// extract X locations from the input file
    }
    
    for(int i=0; i<(x.size()-1); i++)
        gsl_vector_set(Loc,0+i*2,x[i+1]);								// save X in vector Loc (alternating XY values)

    input.close();  // close the input file

    return;
}

//*************************************************
// Read in y-coordinates
//*************************************************
/*
Option to use a specific landscape with Z patches.
Structure of input files: two columns. Column 1: X locations of patches; column 2: Y locations of patches.
*/
static void getInputLocY(gsl_vector *Loc)
{
	char landsc_in[99] = "landscape_";
	strcpy(landsc_in,directory);
    strcat(landsc_in,"landscapes/landscape_");
    strcat(landsc_in,landscape);
    strcat(landsc_in,".csv");
    

    std::ifstream input(landsc_in);                                       
    std::vector<double> y;
    std::string line;

    if(!input.is_open())
    {
        std::cerr << "There was a problem opening the landscape file\n";
        exit(1);
    }

    double col1,col2;

    while(getline(input,line))
    {
        std::istringstream iss(line);
        iss >> col1 >> col2;
        getline(input, line, ',');  
        y.push_back(col1);												// extract Y locations from the input file
     }

    for(int i=0; i<(y.size()-1); i++)
        gsl_vector_set(Loc,1+i*2,y[i+1]);								// save Y in vector Loc (alternating XY values)

    input.close();  // close the input file

    return;
}


//***************************************************************************************************
// Helper functions 
//***************************************************************************************************
/*
1. Function show_matrix prints a matrix with m rows and n cols in the terminal
*/
static void show_matrix(gsl_matrix *M, int rows, int cols)
{

    int i, j;
    for(i=0; i<(rows); i++)
    {
        for(j=0; j<(cols); j++)
        {

            if(j==cols-1)
                printf("%.3g \n",
                       (double)gsl_matrix_get(M, i, j));
            else
                printf("%.3g \t",
                       (double)gsl_matrix_get(M, i, j));
        }
    }
    printf("\n"); 
    
    return;
}

/*
2. Function show_vector prints a vector with of length m in the terminal
*/
static void show_vector(gsl_vector *A, int Num)
{
    int i,j;
    for(i = 0; i<Num; i++)
    {
        printf("%g \n",gsl_vector_get(A,i));
    }
    printf("\n"); 
    return;
}

/*
3. Calculate standard deviation of a gsl_vector
Input: gsl_vector, 3 integer numbers: start (START), end (P), sample size (div)
*/
double calc_sd(gsl_vector *con, int START, int P, int div)
{
    double sum = 0.0, mean = 0, sd = 0.0;
    int i;

    for(i = START; i < P; ++i)
        sum += gsl_vector_get(con,i);

    mean = sum/div;                                                             // calculate the arithmetic mean 

    for(i = START; i < P; ++i)
        sd += pow(gsl_vector_get(con,i) - mean, 2);

    return sqrt(sd/div);
}

/*
4. Calcuatlate the arithmetic mean of a gsl_vector 
Input: gsl_vector, 3 integer numbers: start (START), end (P), sample size (div)
*/
double calc_mean(gsl_vector *con, int START, int P, int div)
{
    double sum = 0.0;
    int i;

    for(i = START; i < P; ++i)
        sum += gsl_vector_get(con,i);

    return (sum/div);
}

/* 
5. Get a random number from a Gaussian distribution with specified parameters and within specified limits
Input: mean, variance, minimum value, maximum value, initalised gsl random number generator r
*/
double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff ,gsl_rng *r)
{
    double par = low_cutoff - 1;                                          // initialise parameter with a value outside the desired range

    while(par < low_cutoff || par > high_cutoff)
        par = gsl_ran_gaussian(r,sigma) + mean;

    return par;
}
