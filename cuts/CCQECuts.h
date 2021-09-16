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

  //Not sure how I wish to handle this... might use Event... might call function
  //3 times... Only worry is if either is a reasonable sideband that I wouldn't want
  //the two criteria to be appleid together...  
    /*  
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
  class MaxEMBlobs: public PlotUtils::Cut<UNIVERSE, EVENT>
  {
  public:
    // Constructor
    MaxEMBlobs(const std::string& name, const double nMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMax(nMax) 
    {
    }

  private:
    bool checkCut(const UNIVERSE& univ, EVENT& evt) const override
    {
      return evt.GetNEMBlobs() >= fMax;
    }

    const double fMax;
  };
*/

}
#endif //David_CCQECuts_H
