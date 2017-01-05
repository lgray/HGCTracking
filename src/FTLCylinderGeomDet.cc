#include "RecoTracker/FastTimeMatching/interface/FTLCylinderGeomDet.h"
#include "DataFormats/GeometrySurface/interface/MediumProperties.h"
#include "DataFormats/GeometrySurface/interface/BoundPlane.h"

#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

FTLBarrelGeomDet::FTLBarrelGeomDet(int subdet, int zside, int layer, float z, float rmin, float rmax, float zmin, float zmax, float radlen, float xi) :
  GeomDet( Plane::build(Plane::PositionType(0,0,z), Plane::RotationType(), RectangularPlaneBounds(0.5*(rmin + rmax),0.5*(rmin + rmax), zmax-zmin).clone()).get() ),
  subdet_(subdet), zside_(zside), layer_(layer) 
{
  setDetId(FastTimeDetId(subdet,0,0,zside));
  
  if (radlen > 0) {
    (const_cast<Plane &>(surface())).setMediumProperties(MediumProperties(radlen,xi));
  }
}

