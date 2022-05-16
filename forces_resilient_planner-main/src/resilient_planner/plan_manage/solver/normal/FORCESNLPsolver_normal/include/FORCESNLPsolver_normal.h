/*
FORCESNLPsolver_normal : A fast customized optimization solver.

Copyright (C) 2013-2021 EMBOTECH AG [info@embotech.com]. All rights reserved.


This software is intended for simulation and testing purposes only. 
Use of this software for any commercial purpose is prohibited.

This program is distributed in the hope that it will be useful.
EMBOTECH makes NO WARRANTIES with respect to the use of the software 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. 

EMBOTECH shall not have any liability for any damage arising from the use
of the software.

This Agreement shall exclusively be governed by and interpreted in 
accordance with the laws of Switzerland, excluding its principles
of conflict of laws. The Courts of Zurich-City shall have exclusive 
jurisdiction in case of any dispute.

*/

/* Generated by FORCESPRO v4.4.0 on Monday, June 28, 2021 at 3:20:00 PM */
#ifndef FORCESNLPsolver_normal_H
#define FORCESNLPsolver_normal_H

#ifndef SOLVER_STDIO_H
#define SOLVER_STDIO_H
#include <stdio.h>
#endif


#ifndef SOLVER_STANDARD_TYPES
#define SOLVER_STANDARD_TYPES

typedef signed char solver_int8_signed;
typedef unsigned char solver_int8_unsigned;
typedef char solver_int8_default;
typedef signed short int solver_int16_signed;
typedef unsigned short int solver_int16_unsigned;
typedef short int solver_int16_default;
typedef signed int solver_int32_signed;
typedef unsigned int solver_int32_unsigned;
typedef int solver_int32_default;
typedef signed long long int solver_int64_signed;
typedef unsigned long long int solver_int64_unsigned;
typedef long long int solver_int64_default;

#endif


/* DATA TYPE ------------------------------------------------------------*/
typedef double FORCESNLPsolver_normal_float;
typedef double FORCESNLPsolver_normal_callback_float;

typedef double FORCESNLPsolver_normalinterface_float;

/* SOLVER SETTINGS ------------------------------------------------------*/

/* MISRA-C compliance */
#ifndef MISRA_C_FORCESNLPsolver_normal
#define MISRA_C_FORCESNLPsolver_normal (0)
#endif

/* restrict code */
#ifndef RESTRICT_CODE_FORCESNLPsolver_normal
#define RESTRICT_CODE_FORCESNLPsolver_normal (0)
#endif

/* print level */
#ifndef SET_PRINTLEVEL_FORCESNLPsolver_normal
#define SET_PRINTLEVEL_FORCESNLPsolver_normal    (1)
#endif

/* timing */
#ifndef SET_TIMING_FORCESNLPsolver_normal
#define SET_TIMING_FORCESNLPsolver_normal    (1)
#endif

/* Numeric Warnings */
/* #define PRINTNUMERICALWARNINGS */

/* maximum number of iterations  */
#define SET_MAXIT_FORCESNLPsolver_normal			(200)	

/* scaling factor of line search (FTB rule) */
#define SET_FLS_SCALE_FORCESNLPsolver_normal		(FORCESNLPsolver_normal_float)(0.99)      

/* maximum number of supported elements in the filter */
#define MAX_FILTER_SIZE_FORCESNLPsolver_normal	(200) 

/* maximum number of supported elements in the filter */
#define MAX_SOC_IT_FORCESNLPsolver_normal			(4) 

/* desired relative duality gap */
#define SET_ACC_RDGAP_FORCESNLPsolver_normal		(FORCESNLPsolver_normal_float)(0.0001)

/* desired maximum residual on equality constraints */
#define SET_ACC_RESEQ_FORCESNLPsolver_normal		(FORCESNLPsolver_normal_float)(1E-06)

/* desired maximum residual on inequality constraints */
#define SET_ACC_RESINEQ_FORCESNLPsolver_normal	(FORCESNLPsolver_normal_float)(1E-06)

/* desired maximum violation of complementarity */
#define SET_ACC_KKTCOMPL_FORCESNLPsolver_normal	(FORCESNLPsolver_normal_float)(1E-06)


/* SOLVER RETURN CODES----------------------------------------------------------*/
/* solver has converged within desired accuracy */
#define OPTIMAL_FORCESNLPsolver_normal      (1)

/* maximum number of iterations has been reached */
#define MAXITREACHED_FORCESNLPsolver_normal (0)

/* solver has stopped due to a timeout */
#define TIMEOUT_FORCESNLPsolver_normal   (2)

/* wrong number of inequalities error */
#define INVALID_NUM_INEQ_ERROR_FORCESNLPsolver_normal  (-4)

/* factorization error */
#define FACTORIZATION_ERROR_FORCESNLPsolver_normal   (-5)

/* NaN encountered in function evaluations */
#define BADFUNCEVAL_FORCESNLPsolver_normal  (-6)

/* no progress in method possible */
#define NOPROGRESS_FORCESNLPsolver_normal   (-7)

/* invalid values in parameters */
#define PARAM_VALUE_ERROR_FORCESNLPsolver_normal   (-11)

/* too small timeout given */
#define INVALID_TIMEOUT_FORCESNLPsolver_normal   (-12)

/* licensing error - solver not valid on this machine */
#define LICENSE_ERROR_FORCESNLPsolver_normal  (-100)

/* INTEGRATORS RETURN CODE ------------*/
/* Integrator ran successfully */
#define INTEGRATOR_SUCCESS (11)
/* Number of steps set by user exceeds maximum number of steps allowed */
#define INTEGRATOR_MAXSTEPS_EXCEEDED (12)





/* PARAMETERS -----------------------------------------------------------*/
/* fill this with data before calling the solver! */
typedef struct
{
    /* vector of size 9 */
    FORCESNLPsolver_normal_float xinit[9];

    /* vector of size 340 */
    FORCESNLPsolver_normal_float x0[340];

    /* vector of size 2600 */
    FORCESNLPsolver_normal_float all_parameters[2600];

    /* scalar */
    solver_int32_unsigned num_of_threads;


} FORCESNLPsolver_normal_params;  ///


/* OUTPUTS --------------------------------------------------------------*/
/* the desired variables are put here by the solver */
typedef struct
{
    /* vector of size 17 */
    FORCESNLPsolver_normal_float x01[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x02[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x03[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x04[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x05[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x06[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x07[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x08[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x09[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x10[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x11[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x12[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x13[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x14[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x15[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x16[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x17[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x18[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x19[17];

    /* vector of size 17 */
    FORCESNLPsolver_normal_float x20[17];


} FORCESNLPsolver_normal_output;  ///


/* SOLVER INFO ----------------------------------------------------------*/
/* diagnostic data from last interior point step */
typedef struct
{
    /* iteration number */
    solver_int32_default it;

	/* number of iterations needed to optimality (branch-and-bound) */
	solver_int32_default it2opt;
	
    /* inf-norm of equality constraint residuals */
    FORCESNLPsolver_normal_float res_eq;
	
    /* inf-norm of inequality constraint residuals */
    FORCESNLPsolver_normal_float res_ineq;

	/* norm of stationarity condition */
    FORCESNLPsolver_normal_float rsnorm;

	/* max of all complementarity violations */
    FORCESNLPsolver_normal_float rcompnorm;

    /* primal objective */
    FORCESNLPsolver_normal_float pobj;	
	
    /* dual objective */
    FORCESNLPsolver_normal_float dobj;	

    /* duality gap := pobj - dobj */
    FORCESNLPsolver_normal_float dgap;		
	
    /* relative duality gap := |dgap / pobj | */
    FORCESNLPsolver_normal_float rdgap;		

    /* duality measure */
    FORCESNLPsolver_normal_float mu;

	/* duality measure (after affine step) */
    FORCESNLPsolver_normal_float mu_aff;
	
    /* centering parameter */
    FORCESNLPsolver_normal_float sigma;
	
    /* number of backtracking line search steps (affine direction) */
    solver_int32_default lsit_aff;
    
    /* number of backtracking line search steps (combined direction) */
    solver_int32_default lsit_cc;
    
    /* step size (affine direction) */
    FORCESNLPsolver_normal_float step_aff;
    
    /* step size (combined direction) */
    FORCESNLPsolver_normal_float step_cc;    

	/* solvertime */
	FORCESNLPsolver_normal_float solvetime;   

	/* time spent in function evaluations */
	FORCESNLPsolver_normal_float fevalstime;  


} FORCESNLPsolver_normal_info;   ///









/* SOLVER FUNCTION DEFINITION -------------------------------------------*/
/* Time of Solver Generation: (UTC) Monday, June 28, 2021 3:20:01 PM */
/* User License expires on: (UTC) Friday, October 29, 2021 10:00:00 PM (approx.) (at the time of code generation) */
/* Solver Static License expires on: (UTC) Friday, October 29, 2021 10:00:00 PM (approx.) */
/* Solver Generation Request Id: 24aaf738-cd2a-4f0e-a19e-3998c0f9fce2 */
/* examine exitflag before using the result! */
#ifdef __cplusplus
extern "C" {
#endif		

typedef void (*FORCESNLPsolver_normal_extfunc)(FORCESNLPsolver_normal_float* x, FORCESNLPsolver_normal_float* y, FORCESNLPsolver_normal_float* lambda, FORCESNLPsolver_normal_float* params, FORCESNLPsolver_normal_float* pobj, FORCESNLPsolver_normal_float* g, FORCESNLPsolver_normal_float* c, FORCESNLPsolver_normal_float* Jeq, FORCESNLPsolver_normal_float* h, FORCESNLPsolver_normal_float* Jineq, FORCESNLPsolver_normal_float* H, solver_int32_default stage, solver_int32_default iterations, solver_int32_default threadID);

extern solver_int32_default FORCESNLPsolver_normal_solve(FORCESNLPsolver_normal_params *params, FORCESNLPsolver_normal_output *output, FORCESNLPsolver_normal_info *info, FILE *fs, FORCESNLPsolver_normal_extfunc evalextfunctions_FORCESNLPsolver_normal);	

;











#ifdef __cplusplus
}
#endif

#endif
