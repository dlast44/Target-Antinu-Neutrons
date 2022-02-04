//studies includes
#include "studies/Study.h"

//David's includes
#include "event/NeutCands.h"
#include "util/NeutronVariable.h"
#include "util/EventVariable.h"
#include "util/Variable.h"
#include "util/GetBackgroundID.h"
#include "event/CVUniverse.h"

class RecoilSB: public Study
{
  private:
    std::vector<Variable*> fVars;

  public:
    RecoilSB(std::vector<Variable*> vars,
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		     std::vector<CVUniverse*>& data_error_bands): Study()
    {
      for (auto& var : vars){
	fVars.push_back(new Variable((var->GetName()+"_PreRecoilCut").c_str(), var->GetAxisLabel(), var->GetBinVec(),var->GetRecoFunc(),var->GetTrueFunc()));
      }

      for(auto& var: fVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fVars) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& var : fVars) var->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& var : fVars) var->WriteData(outFile);
    }
    
  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      std::bitset<64> SBStat = evt.GetSideBandStat();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = evt.GetLeadingNeutCand().GetPDGBin();
      
      if (evt.IsMC()){
	if (evt.IsSignal()){
	  for (auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    var->selectedSignalReco->FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }
	}
	
	else{
	  int bkgd_ID = -1;	  
	  bkgd_ID = util::GetBackgroundID(univ);
	  for(auto& var: fVars){
	    var->selectedMCReco->FillUniverse(&univ, var->GetRecoValue(univ), weight); //"Fake data" for closure
	    (*var->m_backgroundHists)[bkgd_ID].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(univ), weight);
	  }
	}
      }
      
      else{
	for (auto& var : fVars){
	  var->dataHist->FillUniverse(&univ, var->GetRecoValue(univ), 1);
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
