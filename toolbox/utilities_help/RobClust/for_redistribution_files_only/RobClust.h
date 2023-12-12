/*
 * MATLAB Compiler: 8.4 (R2022a)
 * Date: Wed Sep 21 15:28:00 2022
 * Arguments:
 * "-B""macro_default""-W""lib:RobClust,version=1.0""-T""link:lib""-d""C:\FSDA\u
 * tilities_help\RobClust\for_testing""-v""C:\Users\joshua\DropBoxAldo\Dropbox\_
 * _ALDO\Consulenza\JRC\RobClust.m"
 */

#ifndef RobClust_h
#define RobClust_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_RobClust_C_API 
#define LIB_RobClust_C_API /* No special import/export declaration */
#endif

/* GENERAL LIBRARY FUNCTIONS -- START */

extern LIB_RobClust_C_API 
bool MW_CALL_CONV RobClustInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_RobClust_C_API 
bool MW_CALL_CONV RobClustInitialize(void);

extern LIB_RobClust_C_API 
void MW_CALL_CONV RobClustTerminate(void);

extern LIB_RobClust_C_API 
void MW_CALL_CONV RobClustPrintStackTrace(void);

/* GENERAL LIBRARY FUNCTIONS -- END */

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_RobClust_C_API 
bool MW_CALL_CONV mlxRobClust(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[]);

/* C INTERFACE -- MLX WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

/* C INTERFACE -- MLF WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- START */

extern LIB_RobClust_C_API bool MW_CALL_CONV mlfRobClust(int nargout, mxArray** out, mxArray* yy, mxArray* XX, mxArray* nXX, mxArray* pXX, mxArray* kmax, mxArray* whichClustFun, mxArray* plotClust, mxArray* plotsAutoSel);

#ifdef __cplusplus
}
#endif
/* C INTERFACE -- MLF WRAPPERS FOR USER-DEFINED MATLAB FUNCTIONS -- END */

#endif
