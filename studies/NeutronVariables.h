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
    std::vector<NeutronVariable*> fLeadVars;
    std::vector<NeutronVariable*> fAllVars;
    std::vector<NeutronVariable*> fTgtVars;
    std::vector<NeutronVariable*> fTrackVars;
    double fBound;

  public:
    NeutronVariables(double tgtBoundary, 
		     std::map<std::string, std::vector<CVUniverse*>>& mc_error_bands,
		     std::map<std::string, std::vector<CVUniverse*>>& truth_error_bands,
		     std::vector<CVUniverse*>& data_error_bands): Study(), fBound(tgtBoundary)
    {
      std::vector<double> myBlobEBins;
      const double myBlobEBinWidth = 3.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myBlobEBins.push_back(myBlobEBinWidth * whichBin);

      std::vector<double> mySelfAngleBins;
      const double mySelfAngleBinWidth = 3.2/180.;
      for (int whichBin=0; whichBin < 181; ++whichBin) mySelfAngleBins.push_back(mySelfAngleBinWidth * whichBin);

      std::vector<double> myPDGBins;
      const double myPDGBinWidth = 1.;
      for (int whichBin=0; whichBin < 11; ++whichBin) myPDGBins.push_back(myPDGBinWidth * whichBin);

      std::vector<double> myLenBins;
      const double myLenBinWidth = 10.;
      for (int whichBin=0; whichBin < 51; ++whichBin) myLenBins.push_back(myLenBinWidth * whichBin);

      std::vector<double> myDEDXBins;
      const double myDEDXBinWidth = 2.;
      for (int whichBin=0; whichBin < 26; ++whichBin) myDEDXBins.push_back(myDEDXBinWidth * whichBin);

      std::vector<double> myVtxDistBins;
      const double myVtxDistBinWidth = 10.;
      for (int whichBin=0; whichBin < 101; ++whichBin) myVtxDistBins.push_back(myVtxDistBinWidth * whichBin);

      std::vector<double> myZPosBins;
      const double myZPosBinWidth = 10.;
      for (int whichBin=0; whichBin < ((9300-fBound)/10)+1; ++whichBin) myZPosBins.push_back(myZPosBinWidth * whichBin + fBound);

      fLeadVars = {
	new NeutronVariable("leadBlob_blobE","E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable("leadBlob_SelfAngle","Angle [rad]", mySelfAngleBins,&NeutronCandidates::NeutCand::GetAngleToFP),
	new NeutronVariable("leadBlob_primary_parent","", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable("leadBlob_length","len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable("leadBlob_avg_dEdx","dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable("leadBlob_dist","dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable("leadBlob_Zdist","Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
	new NeutronVariable("leadBlob_ZPos","Z [mm]", myZPosBins,&NeutronCandidates::NeutCand::GetZPos),
      };
      fAllVars = {
	new NeutronVariable("All_blobE","E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable("All_primary_parent","", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable("All_length","len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable("All_avg_dEdx","dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable("All_dist","dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable("All_Zdist","Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
      };
      fTgtVars = {
	new NeutronVariable("target_blobE","E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable("target_primary_parent","", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable("target_length","len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable("target_avg_dEdx","dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable("target_dist","dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable("target_Zdist","Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
      };
      fTrackVars = {
	new NeutronVariable("tracker_blobE","E [MeV]", myBlobEBins,&NeutronCandidates::NeutCand::GetTotalE),
	new NeutronVariable("tracker_primary_parent","", myPDGBins,&NeutronCandidates::NeutCand::GetPDGBin),
	new NeutronVariable("tracker_length","len. [mm]", myLenBins,&NeutronCandidates::NeutCand::GetLength),
	new NeutronVariable("tracker_avg_dEdx","dE/dx [MeV/mm]", myDEDXBins,&NeutronCandidates::NeutCand::GetDEDX),
	new NeutronVariable("tracker_dist","dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxDist),
	new NeutronVariable("tracker_Zdist","Z dist. [mm]", myVtxDistBins,&NeutronCandidates::NeutCand::GetVtxZDist),
      };

      for(auto& var: fLeadVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fLeadVars) var->InitializeDATAHists(data_error_bands);
      for(auto& var: fAllVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fAllVars) var->InitializeDATAHists(data_error_bands);
      for(auto& var: fTgtVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fTgtVars) var->InitializeDATAHists(data_error_bands);
      for(auto& var: fTrackVars) var->InitializeMCHists(mc_error_bands, truth_error_bands);
      for(auto& var: fTrackVars) var->InitializeDATAHists(data_error_bands);
    }

    void SaveOrDraw(TDirectory& outDir){
      return;
    }

    void SaveOrDrawMC(TFile& outFile)
    {
      for (auto& var : fLeadVars) var->WriteMC(outFile);
      for (auto& var : fAllVars) var->WriteMC(outFile);
      for (auto& var : fTgtVars) var->WriteMC(outFile);
      for (auto& var : fTrackVars) var->WriteMC(outFile);
    }

    void SaveOrDrawData(TFile& outFile)
    {
      for (auto& var : fLeadVars) var->WriteData(outFile);
      for (auto& var : fAllVars) var->WriteData(outFile);
      for (auto& var : fTgtVars) var->WriteData(outFile);
      for (auto& var : fTrackVars) var->WriteData(outFile);
    }

  private:
    //Overriding base class functions
    //Do nothing for now...  Good place for data comparisons in the future. 
    void fillSelected(const CVUniverse& univ, const NeutronEvent& evt, const double weight) {
      NeutronCandidates::NeutCand leadCand = evt.GetLeadingNeutCand();
      NeutronCandidates::NeutCands neutCands = evt.GetNeutCands();
      int intType = evt.GetIntType();
      int tgtType = evt.GetTgtZ();
      int leadBlobType = leadCand.GetPDGBin();

      if (evt.IsSignal()){
	for (auto& var : fLeadVars){
	  (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	}
	for (auto& cand: neutCands.GetCandidates()){
	  for (auto& var : fAllVars){
	    (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	  }
	  if (cand.second.GetBegPos().Z() <= fBound){
	    for (auto& var : fTgtVars){
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    }
	  }
	  else {
	    for (auto& var : fTrackVars){
	      (*var->m_SigIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_SigLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    }
	  }
	}
      }

      else {
	for (auto& var : fLeadVars){
	  (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	  (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(leadCand), weight);
	}
	for (auto& cand: neutCands.GetCandidates()){
	  for (auto& var : fAllVars){
	    (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	  }
	  if (cand.second.GetBegPos().Z() <= fBound){
	    for (auto& var : fTgtVars){
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    }
	  }
	  else {
	    for (auto& var : fTrackVars){
	      (*var->m_BkgIntTypeHists)[intType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgTargetTypeHists)[tgtType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	      (*var->m_BkgLeadBlobTypeHists)[leadBlobType].FillUniverse(&univ, var->GetRecoValue(cand.second), weight);
	    }
	  }
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
