//File: BKGFitting.cxx
//Info: This script is intended to fit recoil plots using TFractionFitter
//The following is lifted/translated from git repo MinervaExpt/CCQENu/make_hists/nhv/bkgfitting/FractionFitTest.py
//
//Usage: BKGFitting <mc_file> <data_file> <outdir> optional: <lowFitBinNum> <hiFitBinNum> TODO: Save the information beyond just printing it out
//Author: David Last dlast@sas.upenn.edu/lastd44@gmail.com

//TODO: Same SegFault Business From My Plotting Code... I'm assuming I just need to delete things carefully that I'm not yet.

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
#include "TParameter.h"
#include "TFractionFitter.h"

//PlotUtils includes??? Trying anything at this point...
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"

#ifndef NCINTEX
#include "Cintex/Cintex.h"
#endif

using namespace std;
using namespace PlotUtils;

//TCanvas* DrawRatioFromTH1Ds maybe?

TCanvas* DrawRatioFromMnvH1Ds(MnvH1D* h_data, MnvH1D* h_sig, MnvH1D* h_bkg){

  MnvH1D* mcSum = (MnvH1D*)h_sig->Clone();
  mcSum->Add(h_bkg);
  
  h_sig->SetLineColor(TColor::GetColor("#999933"));
  h_sig->SetFillColor(TColor::GetColor("#999933"));
  h_bkg->SetLineColor(TColor::GetColor("#882255"));
  h_bkg->SetFillColor(TColor::GetColor("#882255"));

  THStack* h = new THStack();
  h->Add((TH1D*)h_bkg->GetCVHistoWithError().Clone());
  h->Add((TH1D*)h_sig->GetCVHistoWithError().Clone());

  TH1D* dataHist = (TH1D*)h_data->GetCVHistoWithError().Clone();

  TCanvas* c1 = new TCanvas("c1","c1",1200,800);
  TPad* top = new TPad("Overlay","Overlay",0,0.078+0.2,1,1);
  TPad* bottom = new TPad("Ratio","Ratio",0,0,1,0.078+0.2);
  top->Draw();
  bottom->Draw();
  top->cd();

  double bottomArea = bottom->GetWNDC()*bottom->GetHNDC();
  double topArea = top->GetWNDC()*top->GetHNDC();

  double areaScale = topArea/bottomArea;

  cout << "areaScale: " << areaScale << endl;

  h->Draw("hist");
  h->SetMaximum((dataHist->GetMaximum())*1.05);
  c1->Update();
  dataHist->Draw("same");
  c1->Update();

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
 
  leg->AddEntry(dataHist,"DATA");
  leg->AddEntry(h_sig,"Signal");
  leg->AddEntry(h_bkg,"BKG");

  leg->Draw();
  c1->Update();

  bottom->cd();
  bottom->SetTopMargin(0.05);
  bottom->SetBottomMargin(0.3);

  MnvH1D* ratio = (MnvH1D*)h_data->Clone();
  ratio->Divide(ratio,mcSum);

  TH1D* mcRatio = new TH1D(mcSum->GetTotalError(false, true, false));
  for (int iBin=1; iBin <= mcRatio->GetXaxis()->GetNbins(); ++iBin){
    mcRatio->SetBinError(iBin, max(mcRatio->GetBinContent(iBin),1.0e-9));
    mcRatio->SetBinContent(iBin,1);
  }

  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(3);
  ratio->SetTitle("");
  //ratio->SetTitleSize(0);                                                     
  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetTitleSize(0.05*areaScale);
  ratio->GetYaxis()->SetTitleOffset(0.75/areaScale);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*areaScale);

  ratio->GetXaxis()->SetLabelSize(ratio->GetXaxis()->GetLabelSize()*areaScale);
  ratio->GetXaxis()->SetTitleSize(0.04*areaScale);
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->Draw();

  mcRatio->SetLineColor(kRed);
  mcRatio->SetLineWidth(3);
  mcRatio->SetFillColorAlpha(kPink + 1, 0.4);
  mcRatio->Draw("E2 SAME");

  TH1D* straightLine = (TH1D*)mcRatio->Clone();
  straightLine->SetFillStyle(0);
  straightLine->Draw("HIST SAME");

  c1->Update();
  return c1;
}

bool PathExists(string path){
  struct stat buffer;
  return (stat (path.c_str(), &buffer) == 0);
}

int main(int argc, char* argv[]) {

  gStyle->SetOptStat(0);

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc < 4 || argc > 6) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName = string(argv[1]);
  string DATAfileName = string(argv[2]);
  string outDir = string(argv[3]);

  //int lowBin = 11;//For 200 MeV for the neutron sample.
  int lowBin = 6;//For 200 MeV for the neutron sample.
  int hiBin = 25;
  double binWidth = 20.0;//hard-coded from the recoil variable for now.
  if (argc > 4) lowBin = atoi(argv[4]);
  if (argc > 5) hiBin = atoi(argv[5]);

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = MCfileName;
  size_t pos=0;

  if (PathExists(outDir)){
    cout << "Thank you for choosing a path for output files that exists." << endl;
  }
  else{
    cout << "Output directory doesn't exist. Exiting" << endl;
    return 3;
  }


  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "MC Input need be .root file." << endl;
    return 4;
  }

  cout << "Input MC file name parsed to: " << fileNameStub << endl;

  rootExt = ".root";
  slash = "/";
  token = "";
  fileNameStub = DATAfileName;
  pos=0;

  //cout << sigNameStub << endl;
  while ((pos = fileNameStub.find(slash)) != string::npos){
    //cout << sigNameStub << endl;
    token = fileNameStub.substr(0, pos);
    //cout << token << endl;
    fileNameStub.erase(0, pos+slash.length());
  }
  //cout << sigNameStub << endl;
  if ((pos=fileNameStub.find(rootExt)) == string::npos){
    cout << "DATA Input need be .root file." << endl;
    return 5;
  }

  cout << "Input Data file name parsed to: " << fileNameStub << endl;

  //cout << "Setting up MnvPlotter" << endl;
  //MnvPlotter* plotter = new MnvPlotter(kCCQEAntiNuStyle);

  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  double POTscale = dataPOT->GetVal()/mcPOT->GetVal();
  //cout << "POT scale factor: " << scale << endl;
  
  TString varName = "recoilE_PreRecoilCut";

  MnvH1D* dataHist = (MnvH1D*)(dataFile->Get(varName+"_data"))->Clone();
  MnvH1D* sigHist = (MnvH1D*)(mcFile->Get(varName+"_selected_signal_reco"))->Clone();
  MnvH1D* chargePiHist = (MnvH1D*)(mcFile->Get(varName+"_background_1chargePi"))->Clone();
  MnvH1D* neutPiHist = (MnvH1D*)(mcFile->Get(varName+"_background_1neutPi"))->Clone();
  MnvH1D* NPiHist = (MnvH1D*)(mcFile->Get(varName+"_background_NPi"))->Clone();
  MnvH1D* otherHist = (MnvH1D*)(mcFile->Get(varName+"_background_Other"))->Clone();

  MnvH1D* bkgTotHist = chargePiHist->Clone();
  bkgTotHist->Add(neutPiHist);
  bkgTotHist->Add(NPiHist);
  bkgTotHist->Add(otherHist);

  sigHist->Scale(POTscale);
  bkgTotHist->Scale(POTscale);

  MnvH1D* mcTotHist = bkgTotHist->Clone();
  mcTotHist->Add(sigHist);

  TH1D* hData = (TH1D*)dataHist->GetCVHistoWithStatError().Clone();
  TH1D* hMC = (TH1D*)mcTotHist->GetCVHistoWithStatError().Clone();
  TH1D* hSig = (TH1D*)sigHist->GetCVHistoWithStatError().Clone();
  TH1D* hBKG = (TH1D*)bkgTotHist->GetCVHistoWithStatError().Clone();

  double dataInt = hData->Integral(lowBin,hiBin);
  double mcInt = hMC->Integral(lowBin,hiBin);
  double sigInt = hSig->Integral(lowBin,hiBin);
  double bkgInt = hBKG->Integral(lowBin,hiBin);

  double sigFrac = sigInt/mcInt;
  double bkgFrac = bkgInt/mcInt;

  cout << "Initial Sig Frac." << sigFrac << endl;

  double scale = dataInt/mcInt;
  hMC->Scale(scale);
  hSig->Scale(scale);
  hBKG->Scale(scale);
  
  TObjArray* mcList = new TObjArray(2);
  mcList->Add(hSig);
  mcList->Add(hBKG);

  TFractionFitter* fit = new TFractionFitter(hData,mcList,"V");
  /* Unclear how to translate into c++ code... very confused since I see this used elsewhere as well... outdated root thing maybe?
  TVirtualFitter* vFit = fit->GetFitter();
  vFit->Config().ParSettings(0).Set("sig", sigFrac, binWidth, 0.0, 1.0);
  vFit->Config().ParSettings(1).Set("bkg", bkgFrac, binWidth, 0.0, 1.0);
  */
  //This constraint requires one to area normalize the totalMChist to the dataHist.
  fit->Constrain(0.0,0.0,1.0);
  fit->Constrain(1.0,0.0,1.0);
  fit->SetRangeX(lowBin,hiBin);

  int status = fit->Fit();
  double scale0, scale1, err0, err1;

  if (status==0){
    cout << "Fit Successful!" << endl;
    fit->GetResult(0,scale0,err0);
    fit->GetResult(1,scale1,err1);
    cout << "Par 0: " << scale0 << " with error: " << err0 << endl;
    cout << "Par 1: " << scale1 << " with error: " << err1 << endl;
  }
  else{
    cout << "FIT FAILED. Exiting..." << endl;
    return 6;
  }

  //scale factor by newly derived signal fraction.
  scale0=scale0/sigFrac;
  scale1=scale1/bkgFrac;
  err0=err0/sigFrac;
  err1=err1/bkgFrac;

  //Parameters to maybe be used later to scale histograms that aren't the one in question. Necessary because of the initial scaling of the histo to the area normalization in the fit region.
  double scale0_full=scale0*scale;
  double scale1_full=scale1*scale;
  double err0_full=err0*scale;
  double err1_full=err1*scale;

  TCanvas* c1 = DrawRatioFromMnvH1Ds(dataHist, sigHist, bkgTotHist);
  TPad* top = (TPad*)c1->GetPrimitive("Overlay");
  c1->Print((TString)outDir+"recoilE_preFit_POTScale.pdf");
  c1->Print((TString)outDir+"recoilE_preFit_POTScale.png");
  top->SetLogy();
  c1->Update();
  c1->Print((TString)outDir+"recoilE_preFit_POTScale_log.pdf");
  c1->Print((TString)outDir+"recoilE_preFit_POTScale_log.png");  
  delete c1;

  //Scaling to the area normalizaion
  sigHist->Scale(scale);
  bkgTotHist->Scale(scale);

  c1 = DrawRatioFromMnvH1Ds(dataHist, sigHist, bkgTotHist);
  top = (TPad*)c1->GetPrimitive("Overlay");
  c1->Print((TString)outDir+"recoilE_preFit_areaScale.pdf");
  c1->Print((TString)outDir+"recoilE_preFit_areaScale.png");
  top->SetLogy();
  c1->Update();
  c1->Print((TString)outDir+"recoilE_preFit_areaScale_log.pdf");
  c1->Print((TString)outDir+"recoilE_preFit_areaScale_log.png");  
  delete c1;

  sigHist->Scale(1.0/scale);
  bkgTotHist->Scale(1.0/scale);

  sigHist->Scale(scale0);
  bkgTotHist->Scale(scale1);

  c1 = DrawRatioFromMnvH1Ds(dataHist, sigHist, bkgTotHist);
  top = (TPad*)c1->GetPrimitive("Overlay");
  c1->Print((TString)outDir+"recoilE_postFit_fitScaleONLY.pdf");
  c1->Print((TString)outDir+"recoilE_postFit_fitScaleONLY.png");
  top->SetLogy();
  c1->Update();
  c1->Print((TString)outDir+"recoilE_postFit_fitScaleONLY_log.pdf");
  c1->Print((TString)outDir+"recoilE_postFit_fitScaleONLY_log.png");  
  delete c1;

  sigHist->Scale(1.0/scale0);
  bkgTotHist->Scale(1.0/scale1);

  sigHist->Scale(scale0_full);
  bkgTotHist->Scale(scale1_full);

  c1 = DrawRatioFromMnvH1Ds(dataHist, sigHist, bkgTotHist);
  top = (TPad*)c1->GetPrimitive("Overlay");
  c1->Print((TString)outDir+"recoilE_postFit_fitAreaScale.pdf");
  c1->Print((TString)outDir+"recoilE_postFit_fitAreaScale.png");
  top->SetLogy();
  c1->Update();
  c1->Print((TString)outDir+"recoilE_postFit_fitAreaScale_log.pdf");
  c1->Print((TString)outDir+"recoilE_postFit_fitAreaScale_log.png");  
  delete c1;

  sigHist->Scale(1.0/scale0_full);
  bkgTotHist->Scale(1.0/scale1_full);

  c1 = DrawRatioFromMnvH1Ds(dataHist, sigHist, bkgTotHist);
  top = (TPad*)c1->GetPrimitive("Overlay");
  c1->Print((TString)outDir+"recoilE_checkScaling.pdf");
  c1->Print((TString)outDir+"recoilE_checkScaling.png");
  top->SetLogy();
  c1->Update();
  c1->Print((TString)outDir+"recoilE_checkScaling_log.pdf");
  c1->Print((TString)outDir+"recoilE_checkScaling_log.png");  
  delete c1;

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
