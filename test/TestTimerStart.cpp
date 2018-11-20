
#include "Timer.h"
#include "def.h"

//#include "matrix.h"
#include <iostream>
using namespace std;
using namespace ROPTLIB;

int main(){
    unsigned long t0;
    double st0;
    
    t0 = getTickCount();
    st0 = static_cast<double>(getTickCount() - t0) / CLK_PS;
    
    #ifdef _WIN64
        #ifdef _DEBUG
	        _CrtDumpMemoryLeaks();
        #endif
    #endif
	return 0;
}

#ifdef MATLAB_MEX_FILE


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{        
    unsigned long t0;
    double st0;
    
    t0 = getTickCount();
    st0 = static_cast<double>(getTickCount() - t0) / CLK_PS;
}


#endif
