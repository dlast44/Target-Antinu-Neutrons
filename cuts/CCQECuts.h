//File: CCQECuts
//Brief: Cuts for CCQE selection derived from Tejin's cuts.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

#ifndef David_CCQECuts_H
#define David_CCQECuts_H

#include "PlotUtils/Cut.h"

namespace MyCCQECuts{

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class PMuRange: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor                                                          
    PMuRange(const std::string& name, const double pMin, const double pMax):
    PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(pMin), fMax(pMax)
    {
    }

  private:
  // THE cut function                                                     
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      // Call a CVUniverse member function to make the cut                
      return univ.GetMuonP() >= fMin && univ.GetMuonP() <= fMax;
    }
    
    const double fMin;
    const double fMax;

  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsAntiNu: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    IsAntiNu(): PlotUtils::Cut<UNIVERSE, EVENT>("Helicity == 2") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNuHelicity() == 2;
    }
  };

  //May want to at some point make this a maximum so that it could
  //be tested or a sideband...
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class IsSingleTrack: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    IsSingleTrack(): PlotUtils::Cut<UNIVERSE, EVENT>("Single Track") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNTracks() == 1;
    }
  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class NoMichels: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    NoMichels(): PlotUtils::Cut<UNIVERSE, EVENT>("No Improved Michels") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      return univ.GetNImprovedMichel() <= 0;
    }
  };

  // Maybe break this apart/ use the EVENT class to handle this cut. Currently implementing as before.
  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class AllEMBlobsCuts: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    AllEMBlobsCuts(): PlotUtils::Cut<UNIVERSE, EVENT>("EM Blobs Cuts") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      std::vector<double> EMBlobInfo = univ.GetEMNBlobsTotalEnergyTotalNHits();
      return EMBlobInfo.at(0) < 2 && EMBlobInfo.at(1) >= 10.0*EMBlobInfo.at(2);
    }

  };

  template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class RecoilCut: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    RecoilCut(): PlotUtils::Cut<UNIVERSE, EVENT>("Recoil Cut") {}

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
    {
      double Q2GeV = univ.GetQ2QEPickledGeV();
      double recoilEGeV = univ.GetDANRecoilEnergyGeV();
      if (Q2GeV < 0.0 || recoilEGeV < 0.0) return false;
      else if (Q2GeV < 0.3) return (recoilEGeV < (0.04+0.43*Q2GeV));
      else if (Q2GeV < 1.4) return (recoilEGeV < (0.08+0.3*Q2GeV));
      else return (recoilEGeV < 0.5);
    }

  };

}
#endif //David_CCQECuts_H
