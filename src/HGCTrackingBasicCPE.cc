#include "RecoParticleFlow/HGCTracking/interface/HGCTrackingBasicCPE.h"

HGCTrackingBasicCPE::HGCTrackingBasicCPE(const CaloGeometry* geom, float fixedError, float fixedErrorBH)  : 
    geom_(geom), 
    geomEE_(static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(HGCalDetId(HGCEE,1,1,0,0,0)))),
    geomFH_(static_cast<const HGCalGeometry*>(geom->getSubdetectorGeometry(HGCalDetId(HGCHEF,1,1,0,0,0)))),
    geomBH_(geom->getSubdetectorGeometry(DetId::Hcal, 2)),
    fixedError_(fixedError), fixedError2_(fixedError*fixedError),
    fixedErrorBH_(fixedErrorBH), fixedErrorBH2_(fixedErrorBH*fixedErrorBH) 
{
}

