#ifndef RecoTracker_FastTimeMatching_ftldebug_h
#define RecoTracker_FastTimeMatching_ftldebug_h

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"

namespace ftltracking { 
    extern int g_debuglevel;
}

typedef std::unordered_multimap<uint32_t,std::pair<const CaloParticle *,float> > CaloTruthRevMap;

#endif
