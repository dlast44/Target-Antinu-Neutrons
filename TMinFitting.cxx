//File: TMinFitting.cxx
//Info: This script is intended to fit for scale factors (and scale factors only) in whatever distribution.
//
//Usage: TMinFitting <mc_file> <data_file>
//
//IGNORE FOR THE TIME BEING <outdir> <do fits in bins of muon momentum (only 0 means no)> optional: <lowFitBinNum> <hiFitBinNum> TODO: Save the information beyond just printing it out
//
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
#include "Math/IFunction.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "TMinuitMinimizer.h"

//Analysis includes
#include "fits/ScaleFactors.h"

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

  string MCfileName = string(argv[1]);
  string DATAfileName = string(argv[2]);
  
  /*
  int lowBin1 = 6;
  int hiBin1 = 25;
  */

  int lowBin = 25;
  int hiBin = 50;

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

  TFile* mcFile = new TFile(MCfileName.c_str(),"READ");
  TFile* dataFile = new TFile(DATAfileName.c_str(),"READ");

  TParameter<double>* mcPOT = (TParameter<double>*)mcFile->Get("POTUsed");
  TParameter<double>* dataPOT = (TParameter<double>*)dataFile->Get("POTUsed");

  std::vector<TString> tags = {"_PreRecoilCut"};

  double POTscale = dataPOT->GetVal()/mcPOT->GetVal();
  //cout << "POT scale factor: " << scale << endl;
  
  TString varName = "recoilE";

  for (int iTag=0; iTag < tags.size(); ++iTag){

    TString tag = tags.at(iTag);
    TString name = varName+tag;

    cout << "Testing Fitting for: " << name << endl;

    MnvH1D* dataHist = (MnvH1D*)(dataFile->Get(name+"_data"))->Clone();
    MnvH1D* sigHist = (MnvH1D*)(mcFile->Get(name+"_selected_signal_reco"))->Clone();
    sigHist->Scale(POTscale);

    MnvH1D* chargePiHist = (MnvH1D*)(mcFile->Get(name+"_background_1chargePi"))->Clone();
    chargePiHist->Scale(POTscale);
    MnvH1D* neutPiHist = (MnvH1D*)(mcFile->Get(name+"_background_1neutPi"))->Clone();
    neutPiHist->Scale(POTscale);
    MnvH1D* NPiHist = (MnvH1D*)(mcFile->Get(name+"_background_NPi"))->Clone();
    NPiHist->Scale(POTscale);
    MnvH1D* otherHist = (MnvH1D*)(mcFile->Get(name+"_background_Other"))->Clone();
    otherHist->Scale(POTscale);

    MnvH1D* QEHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_QE"))->Clone();
    QEHist->Scale(POTscale);
    MnvH1D* RESHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_RES"))->Clone();
    RESHist->Scale(POTscale);
    MnvH1D* DISHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_DIS"))->Clone();
    DISHist->Scale(POTscale);
    MnvH1D* MECHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_2p2h"))->Clone();
    MECHist->Scale(POTscale);
    MnvH1D* OtherIntTypeHist = (MnvH1D*)(mcFile->Get(name+"_bkg_IntType_Other"))->Clone();
    OtherIntTypeHist->Scale(POTscale);

    MnvH1D* bkgTotHist = chargePiHist->Clone();
    bkgTotHist->Add(neutPiHist);
    bkgTotHist->Add(NPiHist);
    bkgTotHist->Add(otherHist);

    MnvH1D* bkg1PiHist = chargePiHist->Clone();
    bkg1PiHist->Add(neutPiHist);

    TH1D* dataHistCV = (TH1D*)dataHist->GetCVHistoWithStatError().Clone();

    vector<TH1D*> fitHistsCV1 = {(TH1D*)bkg1PiHist->GetCVHistoWithStatError().Clone()};
    vector<TH1D*> unfitHistsCV1 = {(TH1D*)sigHist->GetCVHistoWithStatError().Clone()};
    unfitHistsCV1.push_back((TH1D*)NPiHist->GetCVHistoWithStatError().Clone());
    unfitHistsCV1.push_back((TH1D*)otherHist->GetCVHistoWithStatError().Clone());

    vector<TH1D*> fitHistsCV2 = {(TH1D*)bkg1PiHist->GetCVHistoWithStatError().Clone()};
    fitHistsCV2.push_back((TH1D*)NPiHist->GetCVHistoWithStatError().Clone());
    vector<TH1D*> unfitHistsCV2 = {(TH1D*)sigHist->GetCVHistoWithStatError().Clone()};
    unfitHistsCV2.push_back((TH1D*)otherHist->GetCVHistoWithStatError().Clone());

    vector<TH1D*> fitHistsCV3 = {(TH1D*)neutPiHist->GetCVHistoWithStatError().Clone()};
    fitHistsCV3.push_back((TH1D*)chargePiHist->GetCVHistoWithStatError().Clone());
    fitHistsCV3.push_back((TH1D*)NPiHist->GetCVHistoWithStatError().Clone());
    vector<TH1D*> unfitHistsCV3 = {(TH1D*)sigHist->GetCVHistoWithStatError().Clone()};
    unfitHistsCV3.push_back((TH1D*)otherHist->GetCVHistoWithStatError().Clone());

    fit::ScaleFactors func1(fitHistsCV1,unfitHistsCV1,dataHistCV,lowBin,hiBin);
    fit::ScaleFactors func2(fitHistsCV2,unfitHistsCV2,dataHistCV,lowBin,hiBin);
    fit::ScaleFactors func3(fitHistsCV3,unfitHistsCV3,dataHistCV,lowBin,hiBin);
    
    auto* mini1 = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
    auto* mini2 = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);
    auto* mini3 = new ROOT::Minuit2::Minuit2Minimizer(ROOT::Minuit2::kMigrad);

    int nextPar = 0;
    for(unsigned int i=0; i < fitHistsCV1.size(); ++i){
      string var = "par"+to_string(i);
      mini1->SetVariable(nextPar,var,1.0,1.0);
      nextPar++;
    }

    if (nextPar != func1.NDim()){
      cout << "The number of parameters was unexpected for some reason for fitHists1." << endl;
      return 6;
    }

    nextPar = 0;
    for(unsigned int i=0; i < fitHistsCV2.size(); ++i){
      string var = "par"+to_string(i);
      mini2->SetVariable(nextPar,var,1.0,1.0);
      nextPar++;
    }

    if (nextPar != func2.NDim()){
      cout << "The number of parameters was unexpected for some reason for fitHists2." << endl;
      return 6;
    }

    nextPar = 0;
    for(unsigned int i=0; i < fitHistsCV3.size(); ++i){
      string var = "par"+to_string(i);
      mini3->SetVariable(nextPar,var,1.0,1.0);
      nextPar++;
    }

    if (nextPar != func3.NDim()){
      cout << "The number of parameters was unexpected for some reason for fitHists3." << endl;
      return 6;
    }

    mini1->SetFunction(func1);
    mini2->SetFunction(func2);
    mini3->SetFunction(func3);

    cout << "Fitting 1" << endl;
    if(!mini1->Minimize()){
      cout << "Printing Results." << endl;
      mini1->PrintResults();
      cout << "FIT FAILED" << endl;
      //return 7;
    }
    else{
      cout << "Printing Results." << endl;
      mini1->PrintResults();
      //cout << mini->X() << endl;
      //cout << mini->Errors() << endl;
      cout << "FIT SUCCEEDED" << endl;
    }

    cout << "Fitting 2" << endl;
    if(!mini2->Minimize()){
      cout << "Printing Results." << endl;
      mini2->PrintResults();
      cout << "FIT FAILED" << endl;
      //return 7;
    }
    else{
      cout << "Printing Results." << endl;
      mini2->PrintResults();
      //cout << mini->X() << endl;
      //cout << mini->Errors() << endl;
      cout << "FIT SUCCEEDED" << endl;
    }

    cout << "Fitting 3" << endl;
    if(!mini3->Minimize()){
      cout << "Printing Results." << endl;
      mini3->PrintResults();
      cout << "FIT FAILED" << endl;
      //return 7;
    }
    else{
      cout << "Printing Results." << endl;
      mini3->PrintResults();
      //cout << mini->X() << endl;
      //cout << mini->Errors() << endl;
      cout << "FIT SUCCEEDED" << endl;
    }
  }

  cout << "Closing Files... Does this solve the issue of seg fault." << endl;
  mcFile->Close();
  dataFile->Close();

  cout << "HEY YOU DID IT!!!" << endl;
  return 0;
}
