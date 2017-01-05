#include "RecoTracker/FastTimeMatching/interface/FTLTracker.h"
#include "RecoTracker/FastTimeMatching/interface/ftldebug.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "FWCore/Utilities/interface/Exception.h"


FTLTracker::FTLTracker(const FastTimeGeometry *geom) :
    geom_(geom)
{
  if (ftltracking::g_debuglevel > 0) std::cout << "Making the FTLTracker" << std::endl;
  //makeBarrel(FastTimeDetId::FastTimeBarrel, 1); // type , number of layers
  makeDisks(FastTimeDetId::FastTimeEndcap, 1);

  auto ptrSort = [](const FTLDiskGeomDet *a, const FTLDiskGeomDet *b) -> bool { return (*a) < (*b); };
  std::sort(disksPos_.begin(), disksPos_.end(), ptrSort);
  std::sort(disksNeg_.begin(), disksNeg_.end(), ptrSort);
}

void FTLTracker::makeDisks(int subdet, int disks)
{    
  std::vector<float>  rmax(disks, 0), rmin(disks, 9e9);
  std::vector<double> zsumPos(disks), zsumNeg(disks);
  std::vector<int>    countPos(disks), countNeg(disks);
  const std::vector<DetId> & ids = geom_->getValidDetIds();
  if (ftltracking::g_debuglevel > 0) std::cout << "on subdet " << subdet << " I got a total of " << ids.size() << " det ids " << std::endl;
  for (auto & i : ids) {
    const GlobalPoint & pos = geom_->getPosition(i); 
    float z = pos.z();
    float rho = pos.perp();
    int side = z > 0 ? +1 : -1;
    int layer = 0;
    if (subdet == 5) {
      if (std::abs(z) < 400 || std::abs(z) > 600 || (layer == 4 && rho > 240)) continue; // bad, not BH
    }
    (side > 0 ? zsumPos : zsumNeg)[layer] += z;
    (side > 0 ? countPos : countNeg)[layer]++;
    if (rho > rmax[layer]) rmax[layer] = rho;
    if (rho < rmin[layer]) rmin[layer] = rho;
  }
  for (int i = 0; i < disks; ++i) {
    float radlen=0.65f, xi=radlen * 0.42e-3; // see DataFormats/GeometrySurface/interface/MediumProperties.h    
    if (countPos[i]) {
      printf("Positive disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i]);
      addDisk(new FTLDiskGeomDet(subdet, +1, i+1, zsumPos[i]/countPos[i], rmin[i], rmax[i], radlen, xi));
    }
    if (countNeg[i]) {
      printf("Negative disk %2d at z = %+7.2f   %6.1f <= rho <= %6.1f\n", i+1, zsumNeg[i]/countPos[i], rmin[i], rmax[i]);
      addDisk(new FTLDiskGeomDet(subdet, -1, i+1, zsumNeg[i]/countNeg[i], rmin[i], rmax[i], radlen, xi));
    }
  }
}


const FTLDiskGeomDet * FTLTracker::nextDisk(const FTLDiskGeomDet *from, PropagationDirection direction) const 
{
    const std::vector<FTLDiskGeomDet *> & vec = (from->zside() > 0 ? disksPos_ : disksNeg_);
    auto it = std::find(vec.begin(), vec.end(), from);
    if (it == vec.end()) throw cms::Exception("LogicError", "nextDisk called with invalid starting disk");
    if (direction == alongMomentum) {
        if (*it == vec.back()) return nullptr;
        return *(++it);
    } else {
        if (it == vec.begin()) return nullptr;
        return *(--it);
    }
}

const FTLDiskGeomDet * FTLTracker::idToDet(DetId id) const 
{
    // FIXME: can be made less stupid
    for (const FTLDiskGeomDet * disk : disksPos_) {
        if (disk->geographicalId() == id) return disk;
    }
    for (const FTLDiskGeomDet * disk : disksNeg_) {
        if (disk->geographicalId() == id) return disk;
    }
    return nullptr;
}



