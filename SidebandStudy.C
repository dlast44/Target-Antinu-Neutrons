//File: SideBandStudy.C
//Info: Quick macro to find sideband in recoil v. Q2 space with the most consistency with background distributions in pTmu in the signal region.
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//C++ includes
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <unordered_map>
#include <bitset>
#include <time.h>
#include <sys/stat.h>

//ROOT includes
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"
#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TMath.h"
#include "TColor.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"

using namespace std;
using namespace PlotUtils;

double Chi2(TH1D* h1, TH1D* h2){
  h1->Scale(1.0/h1->Integral(0,-1));
  h2->Scale(1.0/h2->Integral(0,-1));
  double chi2 = -999.0;
  if (h1->GetNbinsX() != h2->GetNbinsX() ) return chi2;
  chi2=0.0;
  int nBins = h1->GetNbinsX();
  for (int iBin=1; iBin <= nBins; ++iBin){
    double con1 = h1->GetBinContent(iBin);
    double con2 = h2->GetBinContent(iBin);
    double diff = con1-con2;
    double sum = con1+con2;
    if (sum > 1e-10) chi2 += (diff*diff)/sum;
  }
  return chi2;
}

TH2D* GetTH2D(TString name, TFile* file){
  cout << "Getting: " << name << endl;
  MnvH2D* m = (MnvH2D*)file->Get(name);
  TH2D h1 = m->GetCVHistoWithStatError();
  TH2D* h = (TH2D*)h1.Clone();
  return h;
}

TH1D* GetTH1D(TString name, TFile* file){
  cout << "Getting: " << name << endl;
  MnvH1D* m = (MnvH1D*)file->Get(name);
  TH1D h1 = m->GetCVHistoWithStatError();
  TH1D* h = (TH1D*)h1.Clone();
  return h;
}

void SidebandStudy() {

  //Define these to be the correct files...
  TString fileName="/minerva/data/users/dlast/TargetNeutronsAna/NewStructure/Recoil2D_Exploration_2022_03_14/MyANA_CCQENu_tracker_submitted_2022_03_14_13_53_MnvTune_v1_RecoilExploration_muonPT_wNeutronCuts_neutKE_10_noSYST/runEventLoopMC_SkippedSyst_MnvTuneV1_FVregion_Tracker_wNeutCuts_neutKE_10.000000.root";

  TFile* inFile = new TFile(fileName,"READ");

  map<int,int> recoilCut;
  recoilCut[0] = 8;
  recoilCut[1] = 8;
  recoilCut[2] = 8;
  recoilCut[3] = 8;
  recoilCut[4] = 8;
  recoilCut[5] = 8;
  recoilCut[6] = 8;
  recoilCut[7] = 8;
  recoilCut[8] = 10;
  recoilCut[9] = 11;
  recoilCut[10] = 14;
  recoilCut[11] = 17;
  recoilCut[12] = 20;
  recoilCut[13] = 23;
  recoilCut[14] = 26;

  map<int,TH1D*> Proj_QE;
  map<int,TH1D*> Proj_RES;
  map<int,TH1D*> Proj_DIS;
  map<int,TH1D*> Proj_2p2h;
  map<int,TH1D*> Proj_OtherInt;

  TH2D* h_QE_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_bkg_IntType_QE",inFile);
  TH2D* h_RES_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_bkg_IntType_RES",inFile);
  TH2D* h_DIS_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_bkg_IntType_DIS",inFile);
  TH2D* h_2p2h_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_bkg_IntType_2p2h",inFile);
  TH2D* h_OtherInt_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_bkg_IntType_Other",inFile);

  TH1D* h_QE = GetTH1D("pTmu_bkg_IntType_QE",inFile);
  TH1D* h_RES = GetTH1D("pTmu_bkg_IntType_RES",inFile);
  TH1D* h_DIS = GetTH1D("pTmu_bkg_IntType_DIS",inFile);
  TH1D* h_2p2h = GetTH1D("pTmu_bkg_IntType_2p2h",inFile);
  TH1D* h_OtherInt = GetTH1D("pTmu_bkg_IntType_Other",inFile);

  map<int,TH1D*> Proj_PiC;
  map<int,TH1D*> Proj_Pi0;
  map<int,TH1D*> Proj_NPi;
  map<int,TH1D*> Proj_Other;

  TH2D* h_PiC_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_by_BKG_Label_1chargePi",inFile);
  TH2D* h_Pi0_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_by_BKG_Label_1neutPi",inFile);
  TH2D* h_NPi_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_by_BKG_Label_NPi",inFile);
  TH2D* h_Other_Pre = GetTH2D("TwoD_recQ2Bin_PreRecoilCut_pTmu_PreRecoilCut_by_BKG_Label_Other",inFile);

  TH1D* h_PiC = GetTH1D("pTmu_background_1chargePi",inFile);
  TH1D* h_Pi0 = GetTH1D("pTmu_background_1neutPi",inFile);
  TH1D* h_NPi = GetTH1D("pTmu_background_NPi",inFile);
  TH1D* h_Other = GetTH1D("pTmu_background_Other",inFile);

  map<int,TH1D*> Proj_BKG;

  TH2D* h_BKG_Pre = (TH2D*)h_PiC_Pre->Clone();
  h_BKG_Pre->Add(h_Pi0_Pre);
  h_BKG_Pre->Add(h_NPi_Pre);
  h_BKG_Pre->Add(h_Other_Pre);

  TH1D* h_BKG = (TH1D*)h_PiC->Clone();
  h_BKG->Add(h_Pi0);
  h_BKG->Add(h_NPi);
  h_BKG->Add(h_Other);

  TH1D* h_Tot_Pre = GetTH1D("recQ2Bin_PreRecoilCut_data",inFile);

  TH1D* h_BKG_x = new TH1D("h_BKG_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_QE_x = new TH1D("h_QE_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_RES_x = new TH1D("h_RES_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_DIS_x = new TH1D("h_DIS_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_2p2h_x = new TH1D("h_2p2h_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_OtherInt_x = new TH1D("h_OtherInt_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_PiC_x = new TH1D("h_PiC_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_Pi0_x = new TH1D("h_Pi0_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_NPi_x = new TH1D("h_NPi_x",";Distance from cut line in bins;#chi^{2}",25,0,25);
  TH1D* h_Other_x = new TH1D("h_Other_x",";Distance from cut line in bins;#chi^{2}",25,0,25);

  for (int A=7; A<8; ++A){
    for(int B=A+1; B<25; ++B){
      TH1D* h_BKG_SB;
      TH1D* h_PiC_SB;
      TH1D* h_Pi0_SB;
      TH1D* h_NPi_SB;
      TH1D* h_Other_SB;
      TH1D* h_QE_SB;
      TH1D* h_RES_SB;
      TH1D* h_DIS_SB;
      TH1D* h_2p2h_SB; 
      TH1D* h_OtherInt_SB;

      double selected = 0.0;

      for(int iQ2=0; iQ2<15; ++iQ2){
	int lowBin = iQ2*50+recoilCut[iQ2]+A;
	int hiBin = min(iQ2*50+recoilCut[iQ2]+B,(iQ2+1)*50);

	TString name = "Junk_"+to_string(iQ2);

	TString name_BKG = name+"BKG";

	TString name_QE = name+"QE";
	TString name_RES = name+"RES";
	TString name_DIS = name+"DIS";
	TString name_2p2h = name+"2p2h";
	TString name_OtherInt = name+"OtherInt";

	TString name_PiC = name+"PiC";
	TString name_Pi0 = name+"Pi0";
	TString name_NPi = name+"NPi";
	TString name_Other = name+"Other";

	Proj_QE[iQ2]=h_QE_Pre->ProjectionY(name_QE, lowBin, hiBin);
	Proj_RES[iQ2]=h_RES_Pre->ProjectionY(name_RES, lowBin, hiBin);
	Proj_DIS[iQ2]=h_DIS_Pre->ProjectionY(name_DIS, lowBin, hiBin);
	Proj_2p2h[iQ2]=h_2p2h_Pre->ProjectionY(name_2p2h, lowBin, hiBin);
	Proj_OtherInt[iQ2]=h_OtherInt_Pre->ProjectionY(name_OtherInt, lowBin, hiBin);

	Proj_PiC[iQ2]=h_PiC_Pre->ProjectionY(name_PiC, lowBin, hiBin);
	Proj_Pi0[iQ2]=h_Pi0_Pre->ProjectionY(name_Pi0, lowBin, hiBin);
	Proj_NPi[iQ2]=h_NPi_Pre->ProjectionY(name_NPi, lowBin, hiBin);
	Proj_Other[iQ2]=h_Other_Pre->ProjectionY(name_Other, lowBin, hiBin);

	Proj_BKG[iQ2]=h_BKG_Pre->ProjectionY(name_BKG, lowBin, hiBin);

	selected+=h_Tot_Pre->Integral(lowBin, hiBin);

	if (iQ2==0){
	  TH1D* h = (TH1D*)(Proj_QE[iQ2]->Clone());
	  h_QE_SB = (TH1D*)(Proj_QE[iQ2]->Clone());
	  h_RES_SB = (TH1D*)(Proj_RES[iQ2]->Clone());
	  h_DIS_SB = (TH1D*)(Proj_DIS[iQ2]->Clone());
	  h_2p2h_SB = (TH1D*)(Proj_2p2h[iQ2]->Clone());
	  h_OtherInt_SB = (TH1D*)(Proj_OtherInt[iQ2]->Clone());

	  h_PiC_SB = (TH1D*)(Proj_PiC[iQ2]->Clone());
	  h_Pi0_SB = (TH1D*)(Proj_Pi0[iQ2]->Clone());
	  h_NPi_SB = (TH1D*)(Proj_NPi[iQ2]->Clone());
	  h_Other_SB = (TH1D*)(Proj_Other[iQ2]->Clone());

	  h_BKG_SB = (TH1D*)(Proj_BKG[iQ2]->Clone());
	}

	else{
	  h_QE_SB->Add(Proj_QE[iQ2]);
	  h_RES_SB->Add(Proj_RES[iQ2]);
	  h_DIS_SB->Add(Proj_DIS[iQ2]);
	  h_2p2h_SB->Add(Proj_2p2h[iQ2]);
	  h_OtherInt_SB->Add(Proj_OtherInt[iQ2]);

	  h_PiC_SB->Add(Proj_PiC[iQ2]);
	  h_Pi0_SB->Add(Proj_Pi0[iQ2]);
	  h_NPi_SB->Add(Proj_NPi[iQ2]);
	  h_Other_SB->Add(Proj_Other[iQ2]);

	  h_BKG_SB->Add(Proj_BKG[iQ2]);
	}	
      }

      cout << "" << endl;
      cout << "Distance past cut A: " << A << " Distance past cut B: " << B << endl;
      cout << "Total BKG % in Sideband: " << 100.0*h_BKG_SB->Integral()/selected << endl;
      cout << "Total stats in Sideband: " << selected << endl;
      cout << "Total BKG stats in Sideband: " << h_BKG_SB->Integral() << endl;
      cout << "Chi2 BKG: " << Chi2(h_BKG_SB,h_BKG) << endl;
      h_BKG_x->Fill(B,Chi2(h_BKG_SB,h_BKG));
      cout << "Chi2 QE: " << Chi2(h_QE_SB,h_QE) << endl;
      h_QE_x->Fill(B,Chi2(h_QE_SB,h_QE));
      cout << "Chi2 RES: " << Chi2(h_RES_SB,h_RES) << endl;
      h_RES_x->Fill(B,Chi2(h_RES_SB,h_RES));
      cout << "Chi2 DIS: " << Chi2(h_DIS_SB,h_DIS) << endl;
      h_DIS_x->Fill(B,Chi2(h_DIS_SB,h_DIS));
      cout << "Chi2 2p2h: " << Chi2(h_2p2h_SB,h_2p2h) << endl;
      h_2p2h_x->Fill(B,Chi2(h_2p2h_SB,h_2p2h));
      cout << "Chi2 Other Int: " << Chi2(h_OtherInt_SB,h_OtherInt) << endl;
      h_OtherInt_x->Fill(B,Chi2(h_OtherInt_SB,h_OtherInt));
      cout << "Chi2 PiC: " << Chi2(h_PiC_SB,h_PiC) << endl;
      h_PiC_x->Fill(B,Chi2(h_PiC_SB,h_PiC));
      cout << "Chi2 Pi0: " << Chi2(h_Pi0_SB,h_Pi0) << endl;
      h_Pi0_x->Fill(B,Chi2(h_Pi0_SB,h_Pi0));
      cout << "Chi2 NPi: " << Chi2(h_NPi_SB,h_NPi) << endl;
      h_NPi_x->Fill(B,Chi2(h_NPi_SB,h_NPi));
      cout << "Chi2 Other: " << Chi2(h_Other_SB,h_Other) << endl;
      h_Other_x->Fill(B,Chi2(h_Other_SB,h_Other));
    }
  }

  gStyle->SetOptStat(0);
  TCanvas* c1 = new TCanvas("c1","",800,800);
  TLegend* leg = new TLegend(0.8,0.8,0.95,0.95);
  TLegend* leg2 = new TLegend(0.8,0.8,0.95,0.95);
  c1->cd();
  h_BKG_x->Draw("hist");
  h_BKG_x->GetYaxis()->SetRangeUser(0,0.5);
  h_BKG_x->SetTitle("All Backgrounds #chi^{2} Scan");
  c1->Print("../SideBandStudy7/SideBandStudy_BKG_Only.pdf");
  c1->Print("../SideBandStudy7/SideBandStudy_BKG_Only.png");
  h_QE_x->SetTitle("Backgrounds By Interaction Type #chi^{2} Scan");
  h_QE_x->GetYaxis()->SetRangeUser(0,0.5);
  h_QE_x->Draw("hist");
  h_2p2h_x->SetLineColor(kBlack);
  h_2p2h_x->Draw("same,hist");
  h_OtherInt_x->SetLineColor(kOrange);
  h_OtherInt_x->Draw("same,hist");
  h_RES_x->SetLineColor(kRed);
  h_RES_x->Draw("same,hist");
  h_DIS_x->SetLineColor(kGreen);
  h_DIS_x->Draw("same,hist");
  leg->AddEntry(h_QE_x,"QE");
  leg->AddEntry(h_RES_x,"RES");
  leg->AddEntry(h_DIS_x,"DIS");
  leg->AddEntry(h_2p2h_x,"2p2h");
  leg->AddEntry(h_OtherInt_x,"Other");
  leg->Draw();
  c1->Print("../SideBandStudy7/SideBandStudy_IntTypes.pdf");
  c1->Print("../SideBandStudy7/SideBandStudy_IntTypes.png");
  h_Other_x->SetLineColor(kBlack);
  h_Other_x->GetYaxis()->SetRangeUser(0,0.5);
  h_Other_x->SetTitle("Backgrounds By FS #chi^{2} Scan");
  h_Other_x->Draw("hist");
  h_PiC_x->Draw("same,hist");
  h_Pi0_x->SetLineColor(kRed);
  h_Pi0_x->Draw("same,hist");
  h_NPi_x->SetLineColor(kGreen);
  h_NPi_x->Draw("same,hist");
  leg2->AddEntry(h_PiC_x,"single #pi^{#pm}");
  leg2->AddEntry(h_Pi0_x,"single #pi^{0}");
  leg2->AddEntry(h_NPi_x,"N #pi");
  leg2->AddEntry(h_Other_x,"Other");
  leg2->Draw();
  c1->Print("../SideBandStudy7/SideBandStudy_BKGTypes.pdf");
  c1->Print("../SideBandStudy7/SideBandStudy_BKGTypes.png");

  return;
}
