//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "event/CVUniverse.h"

class NeutronVariables: public Study
{
  private:
    std::vector<NeutronVariable*> fVars;
  public:
    NeutronVariables(): Study()
    {
      std::vector<double> myBlobEBins;
      const double myBlobEBinWidth = 3.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myBlobEBins.push_back(myBlobEBinWidth * whichBin);
      fVars = {
	new NeutronVariable("leabBlobE","E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE,&NeutronCandidates::NeutCand::GetTotalE),
      };
    }

    void SaveOrDraw(TDirectory& outDir)
    {
      /*
       outDir.cd();
       m_VarToGENIELabel->visit([](HIST& wrapper)
                                {
                                  wrapper.SyncCVHistos();
                                  wrapper.hist->Write();
                                });
      */
       //TODO: You could do plotting here
    }

  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight)
    {
      return;
      //(*m_VarToGENIELabel)[univ.GetInteractionType()].FillUniverse(&univ, fReco(univ, evt), weight);
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }
};
