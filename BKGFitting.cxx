//File: BKGFitting.cxx
//Info: This script is intended to fit recoil plots using TFractionFitter
//The following is lifted/translated from git repo MinervaExpt/CCQENu/make_hists/nhv/bkgfitting/FractionFitTest.py
//
//Usage: BKGFitting <mc_file> <data_file> TODO: Save the information beyond just printing it out
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

int main(int argc, char* argv[]) {

  gStyle->SetOptStat(0);

  #ifndef NCINTEX
  ROOT::Cintex::Cintex::Enable();
  #endif

  //Pass an input file name to this script now
  if (argc != 3) {
    cout << "Check usage..." << endl;
    return 2;
  }

  string MCfileName=string(argv[1]);
  string DATAfileName=string(argv[2]);

  string rootExt = ".root";
  string slash = "/";
  string token;
  string fileNameStub = MCfileName;
  size_t pos=0;

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

  //double scale = dataPOT->GetVal()/mcPOT->GetVal();
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

  int lowBin = 6;
  int hiBin = 25;
  double binWidth = 20.0;//hard-coded from the recoil variable for now.

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
  double val0, val1, err0, err1;

  if (status==0){
    cout << "Fit Successful!" << endl;
    fit->GetResult(0,val0,err0);
    fit->GetResult(1,val1,err1);
    cout << "Par 0: " << val0 << " with error: " << err0 << endl;
    cout << "Par 1: " << val1 << " with error: " << err1 << endl;
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
