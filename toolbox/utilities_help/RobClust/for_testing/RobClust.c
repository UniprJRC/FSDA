/*
 * MATLAB Compiler: 8.4 (R2022a)
 * Date: Wed Sep 21 15:28:00 2022
 * Arguments:
 * "-B""macro_default""-W""lib:RobClust,version=1.0""-T""link:lib""-d""C:\FSDA\u
 * tilities_help\RobClust\for_testing""-v""C:\Users\joshua\DropBoxAldo\Dropbox\_
 * _ALDO\Consulenza\JRC\RobClust.m"
 */

#define EXPORTING_RobClust 1
#include "RobClust.h"

static HMCRINSTANCE _mcr_inst = NULL; /* don't use nullptr; this may be either C or C++ */

#if defined( _MSC_VER) || defined(__LCC__) || defined(__MINGW64__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#define NOMINMAX
#include <windows.h>
#undef interface

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        if (GetModuleFileName(hInstance, path_to_dll, _MAX_PATH) == 0)
            return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultPrintHandler(const char *s)
{
    return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern C block */
#endif

#ifdef __cplusplus
extern "C" { // sbcheck:ok:extern_c
#endif

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0;
    size_t len = 0;
    len = strlen(s);
    written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
    if (len > 0 && s[ len-1 ] != '\n')
        written += mclWrite(2 /* stderr */, "\n", sizeof(char));
    return written;
}

#ifdef __cplusplus
} /* End extern C block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_RobClust_C_API
#define LIB_RobClust_C_API /* No special import/export declaration */
#endif

LIB_RobClust_C_API 
bool MW_CALL_CONV RobClustInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler)
{
    int bResult = 0;
    if (_mcr_inst)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!GetModuleFileName(GetModuleHandle("RobClust"), path_to_dll, _MAX_PATH))
        return false;
    {
        mclCtfStream ctfStream = 
            mclGetEmbeddedCtfStream(path_to_dll);
        if (ctfStream) {
            bResult = mclInitializeComponentInstanceEmbedded(&_mcr_inst,
                                                             error_handler, 
                                                             print_handler,
                                                             ctfStream);
            mclDestroyStream(ctfStream);
        } else {
            bResult = 0;
        }
    }  
    if (!bResult)
    return false;
    return true;
}

LIB_RobClust_C_API 
bool MW_CALL_CONV RobClustInitialize(void)
{
    return RobClustInitializeWithHandlers(mclDefaultErrorHandler, mclDefaultPrintHandler);
}

LIB_RobClust_C_API 
void MW_CALL_CONV RobClustTerminate(void)
{
    if (_mcr_inst)
        mclTerminateInstance(&_mcr_inst);
}

LIB_RobClust_C_API 
void MW_CALL_CONV RobClustPrintStackTrace(void) 
{
    char** stackTrace;
    int stackDepth = mclGetStackTrace(&stackTrace);
    int i;
    for(i=0; i<stackDepth; i++)
    {
        mclWrite(2 /* stderr */, stackTrace[i], sizeof(char)*strlen(stackTrace[i]));
        mclWrite(2 /* stderr */, "\n", sizeof(char)*strlen("\n"));
    }
    mclFreeStackTrace(&stackTrace, stackDepth);
}


LIB_RobClust_C_API 
bool MW_CALL_CONV mlxRobClust(int nlhs, mxArray *plhs[], int nrhs, mxArray *prhs[])
{
    return mclFeval(_mcr_inst, "RobClust", nlhs, plhs, nrhs, prhs);
}

LIB_RobClust_C_API 
bool MW_CALL_CONV mlfRobClust(int nargout, mxArray** out, mxArray* yy, mxArray* XX, 
                              mxArray* nXX, mxArray* pXX, mxArray* kmax, mxArray* 
                              whichClustFun, mxArray* plotClust, mxArray* plotsAutoSel)
{
    return mclMlfFeval(_mcr_inst, "RobClust", nargout, 1, 8, out, yy, XX, nXX, pXX, kmax, whichClustFun, plotClust, plotsAutoSel);
}

