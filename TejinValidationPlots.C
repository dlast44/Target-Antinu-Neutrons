//File: signalBKGStack.cxx
//Info: This is a script to run a loop over all MC breakdown plots in a single histos file and save nice plots from them. Primarily used for validation against older plots.
//
//Usage: signalBKGStack <histos_file> <output_directory> <plot_label>
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

void TejinValidationPlots() {

  TString dir = "/minerva/data/users/dlast/TargetNeutronsAna/NewStructure/TejinValidation/";

  TFile* mcFile = new TFile(dir+"runEventLoopMC.root","READ");
  TFile* dataFile = new TFile(dir+"runEventLoopData.root","READ");

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  c1->cd();

  TParameter<double>* mc_POT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* data_POT = (TParameter<double>*)dataFile->Get("POTUsed");

  double scale = data_POT->GetVal()/mc_POT->GetVal();

  MnvH1D* h_data = (MnvH1D*)dataFile->Get("pmu_data");
  h_data->SetLineColor(kBlack);
  h_data->SetFillColor(kBlack);

  /*
  MnvH1D* h_QE_H = (MnvH1D*)mcFile->Get("pmu_sig_IntType_QE_H");
  h_QE_H->Scale(scale);
  h_QE_H->SetLineColor(kRed);
  h_QE_H->SetFillColor(kRed);
  */

  MnvH1D* h_QE_Other = (MnvH1D*)mcFile->Get("pmu_sig_IntType_QE_Other");
  h_QE_Other->Scale(scale);
  h_QE_Other->SetLineColor(kGreen);
  h_QE_Other->SetFillColor(kGreen);

  MnvH1D* h_2p2h = (MnvH1D*)mcFile->Get("pmu_sig_IntType_2p2h");
  h_2p2h->Scale(scale);
  h_2p2h->SetLineColor(kMagenta-4);
  h_2p2h->SetFillColor(kMagenta-4);

  MnvH1D* h_RES = (MnvH1D*)mcFile->Get("pmu_sig_IntType_RES");
  h_RES->Scale(scale);
  h_RES->SetLineColor(kBlue+1);
  h_RES->SetFillColor(kBlue+1);

  MnvH1D* h_DIS = (MnvH1D*)mcFile->Get("pmu_sig_IntType_DIS");
  h_DIS->Scale(scale);
  h_DIS->SetLineColor(kOrange+2);
  h_DIS->SetFillColor(kOrange+2);

  MnvH1D* h_Other = (MnvH1D*)mcFile->Get("pmu_sig_IntType_Other");
  h_Other->Scale(scale);
  h_Other->SetLineColor(kBlack);
  h_Other->SetFillColor(kBlack);

  MnvH1D* h_back_1PiC = (MnvH1D*)mcFile->Get("pmu_background_1chargePi");
  h_back_1PiC->Scale(scale);
  h_back_1PiC->SetLineColor(kBlue-7);
  h_back_1PiC->SetFillColor(kBlue-7);
  h_back_1PiC->SetFillStyle(3144);

  MnvH1D* h_back_1Pi0 = (MnvH1D*)mcFile->Get("pmu_background_1neutPi");
  h_back_1Pi0->Scale(scale);
  h_back_1Pi0->SetLineColor(kGreen-3);
  h_back_1Pi0->SetFillColor(kGreen-3);
  h_back_1Pi0->SetFillStyle(3144);

  MnvH1D* h_back_NPi = (MnvH1D*)mcFile->Get("pmu_background_NPi");
  h_back_NPi->SetLineColor(kRed-4);
  h_back_NPi->SetFillColor(kRed-4);
  h_back_NPi->SetFillStyle(3144);

  MnvH1D* h_back_Other = (MnvH1D*)mcFile->Get("pmu_background_Other");
  h_back_Other->SetLineColor(kBlack);
  h_back_Other->SetFillColor(kBlack);
  h_back_Other->SetFillStyle(3144);

  THStack* h = new THStack();
  h->Add(h_back_Other);
  h->Add(h_back_NPi);
  h->Add(h_back_1Pi0);
  h->Add(h_back_1PiC);

  h->Add(h_Other);
  h->Add(h_DIS);
  h->Add(h_RES);
  h->Add(h_2p2h);
  h->Add(h_QE_Other);
  //h->Add(h_QE_H);

  h->Draw("hist");
  h_data->Draw("same");
  //c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  leg->AddEntry(h_data,"DATA");

  //leg->AddEntry(h_QE_H,"QE-Like && QE H");
  leg->AddEntry(h_QE_Other,"QE-Like && QE Oth");
  leg->AddEntry(h_2p2h,"QE-Like && 2p2h");
  leg->AddEntry(h_RES,"QE-Like && RES");
  leg->AddEntry(h_DIS,"QE-Like && DIS");
  leg->AddEntry(h_Other,"QE-Like && OTH");

  leg->AddEntry(h_back_1PiC,"single #pi+/- in FS");
  leg->AddEntry(h_back_1Pi0,"single #pi0 in FS");
  leg->AddEntry(h_back_1PiC,"N#pi in FS");
  leg->AddEntry(h_back_1PiC,"other");

  leg->Draw();

  c1->Print(dir+"Plots/pmu_stacked_test.pdf");
  c1->Print(dir+"Plots/pmu_stacked_test.png");
}
