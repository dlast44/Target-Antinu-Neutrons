//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/Variable.h"
#include "event/CVUniverse.h"

class NeutronVariables: public Study
{
  private:
  //std::vector<NeutronVariable*> fAllVars;
    std::vector<NeutronVariable*> fLeadVars;
  public:
    NeutronVariables(std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		     std::vector<CVUniverse*>& data_error_bands): Study()
    {
      std::vector<double> myBlobEBins;
      const double myBlobEBinWidth = 3.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myBlobEBins.push_back(myBlobEBinWidth * whichBin);
      fLeadVars = {
	new NeutronVariable("leadBlobE","E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),//&NeutronCandidates::NeutCand::GetTotalE),
	//new Variable("leadBlobE","E [MeV]", myBlobEBins,&CVUniverse::GetNNeutBlobs),//&NeutronCandidates::NeutCand::GetTotalE),
      };
      for(auto& var: fLeadVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fLeadVars) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& var : fLeadVars) var->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& var : fLeadVars) var->WriteData(outFile);
    }

  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      NeutronCandidates::NeutCand leadCand = evt.GetLeadingNeutCand();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      if (evt.IsSignal()){
	for (auto& var : fLeadVars){
	  (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	}
      }
      else{
	for (auto& var : fLeadVars){
	  (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	}
      }
      return;
    }

    //All of your plots happen here so far.
    void fillSelectedSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight)
    {
      return;
      //(*m_VarToGENIELabel)[univ.GetInteractionType()].FillUniverse(&univ, fReco(univ, evt), weight);
    }

    //Do nothing for now...  Good place for efficiency denominators in the future.
    void fillTruthSignal(const CVUniverse& univ, const NeutronEvent& evt, const double weight) { return; }
};
