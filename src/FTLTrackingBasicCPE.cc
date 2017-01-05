#include "RecoTracker/FastTimeMatching/interface/FTLTrackingBasicCPE.h"

FTLTrackingBasicCPE::FTLTrackingBasicCPE(const FastTimeGeometry* geom, float fixedError, float fixedErrorBH)  : 
    geom_(geom),     
    fixedError_(fixedError), fixedError2_(fixedError*fixedError),
    fixedErrorBH_(fixedErrorBH), fixedErrorBH2_(fixedErrorBH*fixedErrorBH) 
{
}

