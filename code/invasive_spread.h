/* Header to "ivasive_spread.cpp" */ 

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <fstream>     		 	// for input/output from or to files
#include <iostream>			// for input/output on terminal
#include <sstream>
#include <vector>
#include <cstdlib>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>		// random number generator
#include <gsl/gsl_randist.h>		// random number distributions
#include <gsl/gsl_blas.h>		// linear algebra routines
#include <gsl/gsl_linalg.h>		// linear algebra
#include <gsl/gsl_sort_vector.h>	// vector operations
#include <gsl/gsl_odeiv.h>              // ODE solver

#include <cvode/cvode.h>        	/* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  	/* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       	/* prototype for CVDense */

static void getInputMass(gsl_vector *mass);
static void getInputDisp(gsl_vector *D_Max);
static void getInputSpp();
static void getInputZ(); 
static void getInputLocX(gsl_vector *Loc);
static void getInputLocY(gsl_vector *Loc);
static void getInputInvasionBodymass(gsl_vector *mass);
static void getInputInvasionDisp(gsl_vector *D_Max);
static void getInputInvasionInvSpp(gsl_vector *InvSppvec);
static void getInputInvasionBiomass(double B[]);
static void getInputInvasionSpp();

double *web_calc(gsl_rng *r);
static void pdef_structure(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass);

static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, double params[]);
double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff ,gsl_rng *r);
static void Generate_RGG_structure(gsl_vector *mass, gsl_vector *D_Max, gsl_vector *Loc, gsl_matrix *SW, gsl_vector *con, double params[], double RGG_params[]);
static void set_spatial_parameters(gsl_matrix *SW, gsl_vector *Loc, double params[]);
static void initialise_biomass_densities(double B[], double iniB[], double params[], gsl_rng *r);

static void solve_ode(double B[], double meanB[],  double meanB_tot[], double CV[],  double CV_tot[], double params[], gsl_matrix *CovM, double Biomass_Temp[], double net_growth_Temp[]);
static void dynamics(double B[], double Bdot[], void *params);
static void Extinction(double B[], int Num);
static void Extinct_Species(double B[], int Num, int Pat);
static int Cvode_rates(realtype t, N_Vector y, N_Vector ydot, void *params);
static void net_rates(gsl_vector *TBvec, gsl_vector *Ivec, gsl_vector *Dvec, gsl_vector *Xvec, gsl_vector *Emvec, gsl_vector *N_in, gsl_vector *N_out, gsl_vector *N_up, void *params);

static void calc_patch_distances(gsl_matrix *P_Dist, gsl_vector *Loc, double RGG_params[]);
static void calc_spatial_connectance(gsl_matrix *SW, gsl_vector *con, double RGG_params[]);

static void output(gsl_matrix *SW, gsl_vector *mass, gsl_vector *con, double iniB[], double meanB[], double meanB_tot[], double CV[], double CV_tot[], double params[], int paramslength, double Biomass_Temp[], double RGG_params[], double net_growth_Temp[], gsl_vector *InvSppvec, double B[]);
static void Prepare_timeseries_file(FILE *timeseries);
static void Write_timeseries_to_file(double B[], double t, FILE *timeseries);
static void Write_params_to_file(double params[], int paramslength);

double calc_sd(gsl_vector *con, int START, int P, int div);
double calc_mean(gsl_vector *con, int START, int P, int div);

static void show_matrix(gsl_matrix *A, int Num, int N2);
static void show_vector(gsl_vector *A, int Num);
