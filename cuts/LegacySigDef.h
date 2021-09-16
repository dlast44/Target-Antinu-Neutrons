//This is what we'll do for reco
/*template <class UNIVERSE>
using q3Signal = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetTrueq3>;*/ //TODO: Except I didn't actually write this for SignalConstraints :(

//PlotUtils includes
#include "PlotUtils/Cut.h"

//Package includes
#include "event/CVUniverse.h"

template <class UNIVERSE>
class Q3Limit: public PlotUtils::SignalConstraint<UNIVERSE>
{
  public:
    Q3Limit(const double q3Max): PlotUtils::SignalConstraint<UNIVERSE>("q3 < " + std::to_string(q3Max) + " GeV"),
                                 fQ3Max(q3Max)
    {
    }

  private:
    double fQ3Max; //Maximum q3 allowed in GeV/c

    bool checkConstraint(const UNIVERSE& univ) const //override
    {
      double trueq3 = univ.Getq3True()/pow(10,3);
      return trueq3 < fQ3Max;
    }
};
